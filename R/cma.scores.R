cma.scores <- function(cma.alter = NULL,
                       cma.cov,
                       cma.samp,
                       scores = c("CaMP", "logLRT"),
                       cma.data = NULL, ##use if provided, if not, create it within the object
                       coverage = NULL, ##use if provided, if not, create it within the object
                       passenger.rates = t(data.frame(0.55*rep(1.0e-6,25))),
                       ##compute.poisson.BF=FALSE,
                       ##compute.binomial.posterior=FALSE,
                       allow.separate.rates = TRUE,
                       filter.above=0, # mutations per Mb above threshold
                       filter.below=0, # mutations per Mb below threshold
                       filter.threshold=0, # minimal size for applying above filters
                       filter.mutations=0,
                       aa=1e-10, # hyperparameters for binomial posterior
                       bb=1e-10,
                       priorH0=1-300/13020, # hyperparameters for poisson BF
                       prior.a0=100,
                       prior.a1=5,
                       prior.fold=10
                       ){

  if(is.null(coverage))
    {
      ##make sure everything is a data frame
      cma.alter <- as.data.frame(cma.alter)
      cma.cov <- as.data.frame(cma.cov)
      cma.samp <- as.data.frame(cma.samp)

      ##now keep only the mutations in cma.alter
      cma.alter$Type <- as.character(cma.alter$Type)
      cma.alter <- cma.alter[cma.alter$Type=="Mut", ]
      
      ##first create coverage object
      coverage <- make.cov.obj(cma.cov, cma.samp)
    }
  
  ##get number of genes
  number.genes <- nrow(coverage)

  ##get average coverage of genes in discovery screen
  GeneTotalSize <-
    coverage[,"CoverageAtRiskDiscovery.ins.del"]/
    coverage[,"NumberTumorsAnalyzedDiscovery"]
  
  if(is.null(cma.data))
    {
      ##now create mutations object
      MutPass <- make.mut.obj(cma.cov, coverage, cma.alter, pass.rates = passenger.rates)
      cma.data <- MutPass$cma.data
      passenger.rates <- MutPass$passenger.rates
    }

  ##create empty data frame of scores if no data is available
  if (nrow(cma.data)==0){
    out <- data.frame( matrix(0,0,8) )
    colnames(out) <- c("CaMP","neglogPg","logLRT",
                           "logitBinomialPosteriorDriver",
                       "PoissonlogBF","PoissonPosterior","Poissonlmlik0",
                           "Poissonlmlik1")
    return(out)
  }
  
  if( allow.separate.rates & nrow(passenger.rates)==2 ) {
    separate.rates <- TRUE
  } else{
    separate.rates <- FALSE
  }
  
  if ( sum(is.na(cma.data[,"NumberTumorsAnalyzedValidation"]))>0 ) {
    return("I got you finally: NA in NumberTumorsAnalyzedValidation")
  }
#  cma.data <- cma.data[!is.na(cma.data[,"NumberTumorsAnalyzedValidation"]),] #CCC
  
  ColMutDisc <- grep("MutationsDiscovery",colnames(cma.data))
  ColMutVali <- grep("MutationsValidation",colnames(cma.data))
  ColCovDisc <- grep("CoverageAtRiskDiscovery",colnames(cma.data))[1:8]
  ColCovVali <- grep("CoverageAtRiskValidation",colnames(cma.data))[1:8]
  ColCovDiscIndel <- grep("CoverageAtRiskDiscovery.ins.del",colnames(cma.data))
  ColCovValiIndel <- grep("CoverageAtRiskValidation.ins.del",colnames(cma.data))
  
  nT <- length(ColMutDisc)
  nS <- cma.data$NumberTumorsAnalyzed
  nSd <- cma.data$NumberTumorsAnalyzedDiscovery
  nSv <- cma.data$NumberTumorsAnalyzedValidation
  nGG <- number.genes
  xxD <- data.matrix( cma.data[,ColMutDisc] ) # mutations discovery
  xxV <- data.matrix( cma.data[,ColMutVali] ) # mutations validation
  xx <- xxD + xxV
  nG <- nrow(cma.data)
  
  if (nG>1) {
    t1 <- cma.data[,ColCovDisc]
    t2 <- matrix( rep(data.matrix(t1),3), nrow(t1), 3*ncol(t1) )[,c(8*(0:2)+1,8*(0:2)+2,8*(0:2)+3,8*(0:2)+4,8*(0:2)+5,8*(0:2)+6,8*(0:2)+7,8*(0:2)+8)]
    rownames(t2) <- rownames(t1)
    tnsD <- round ( cbind(t2,cma.data[,ColCovDiscIndel]) )
    t1 <- cma.data[,ColCovVali]
    t2 <- matrix( rep(data.matrix(t1),3), nrow(t1), 3*ncol(t1) )[,c(8*(0:2)+1,8*(0:2)+2,8*(0:2)+3,8*(0:2)+4,8*(0:2)+5,8*(0:2)+6,8*(0:2)+7,8*(0:2)+8)]
    rownames(t2) <- rownames(t1)
    tnsV <- round ( cbind(t2,cma.data[,ColCovValiIndel]) )
    tns <- tnsD + tnsV
  }
  if (nG==1) {
    t1 <- cma.data[,ColCovDisc]
    t2 <- as.double( rep(t1,3)[c(8*(0:2)+1,8*(0:2)+2,8*(0:2)+3,8*(0:2)+4,8*(0:2)+5,8*(0:2)+6,8*(0:2)+7,8*(0:2)+8)]  )
    tnsD <- t ( as.matrix( round ( c(t2,cma.data[,ColCovValiIndel]) ) ) )
    t1 <- cma.data[,ColCovVali]
    t2 <- as.double( rep(t1,3)[c(8*(0:2)+1,8*(0:2)+2,8*(0:2)+3,8*(0:2)+4,8*(0:2)+5,8*(0:2)+6,8*(0:2)+7,8*(0:2)+8)]  )
    tnsV <- t ( as.matrix( round ( c(t2,cma.data[,ColCovValiIndel]) ) ) )
    tns <- tnsD + tnsV
  }

  br <- data.matrix( passenger.rates )
  if ( !separate.rates ) br <- colMeans(br)
  if ( separate.rates ) {
    rrD <- matrix(as.numeric(br[1,]),nG,nT,byrow=TRUE)
    rrV <- matrix(as.numeric(br[2,]),nG,nT,byrow=TRUE)
  } else{
    rr <- rrD <- rrV <- matrix(as.numeric(br),nG,nT,byrow=TRUE)
  }

  rownames(tnsD) <- rownames(xxD)
  rownames(tnsV) <- rownames(xxV)
  colnames(tnsD) <- colnames(xxD)
  colnames(tnsV) <- colnames(xxV)
  
  rrD.mle <- xxD / tnsD
  rrV.mle <- xxV / tnsV
  rrD.mle[is.nan(rrD.mle)] = 0
  rrV.mle[is.nan(rrV.mle)] = 0
  
  if ( separate.rates ) {
    ppgtD <- matrix(dbinom(xxD,tnsD,rrD),nG,nT) # point prob by gene and type
    ppgtV <- matrix(dbinom(xxV,tnsV,rrV),nG,nT)
    ppgt <- ppgtD * ppgtV
  } else{
    ppgt <- matrix(dbinom(xxD+xxV,tnsD+tnsV,rr),nG,nT)
  }

  multinomial=FALSE
  if (multinomial){
    tnsDns <- round( cma.data[,ColCovDisc] )
    tnsDindels <- round( cma.data[,ColCovDiscIndel] )
    tnsVns <- round( cma.data[,ColCovVali] )
    tnsVindels <- round( cma.data[,ColCovValiIndel] )
    pp <- matrix(NA,nG,9)
    for (gg in 1:nG) {
      for (ii in 1:8) {
        yy <- xx[gg,(1:3)+3*(ii-1)]
        yy <- c(tnsDns[gg,ii]+tnsVns[gg,ii]-sum(yy),yy)
        pro <- c(br[(1:3)+3*(ii-1)])
        pro <- c(1-sum(pro),pro)
        pp[gg,ii] <- dmultinom(yy, prob=pro)
      }
      pp[gg,9] <- dbinom(xx[gg,25],tnsDindels[gg]+tnsVindels[gg],br[25])
    }
    lppg.mult <- rowSums(log10(pp))
  }

  ppgtD.mle <- matrix(dbinom(xxD,tnsD,rrD.mle),nG,nT) # point prob by gene and type
  ppgtV.mle <- matrix(dbinom(xxV,tnsV,rrV.mle),nG,nT)
  ppgt.mle <- ppgtD.mle * ppgtV.mle

  lppg <- rowSums(log10(ppgt)) # log point prob by gene

  out <- data.frame( matrix(NA,nG,8) )
  rownames(out) <- rownames(cma.data)
  colnames(out) <- c("CaMP","neglogPg","logLRT",
                     "logitBinomialPosteriorDriver",
                     "PoissonlogBF","PoissonPosterior","Poissonlmlik0",
                     "Poissonlmlik1")

  out[,"CaMP"] <- - log10( nGG / rank(lppg) ) - lppg
  out[,"neglogPg"] <- - lppg

  out[,"logLRT"] <- rowSums(log10(ppgt.mle)) - lppg
  camp.pg <- 10^lppg

  # p-values
  if (nrow(passenger.rates)==2) {
    br <- colMeans(passenger.rates,2)
    rr <- matrix(as.numeric(br),nG,nT,byrow=TRUE)
  }

  if ("logitBinomialPosteriorDriver" %in% scores) {
    aa = matrix( aa, nG, nT)
    bb = matrix( bb, nG, nT)
    post <- matrix( pbeta( rr, aa + xx, bb + tns - xx ), nG, nT )
    post <- 1 -  apply( post, 1, prod )
    out[,"logitBinomialPosteriorDriver"] <- log( post / (1-post) )
  }

  if (length(grep("Poisson", scores)) > 0) {
    a0 = matrix( prior.a0, nG, nT)
    a1 = matrix( prior.a1, nG, nT)
    b0 = matrix( rep(prior.a0-1,nT)/br, nG, nT, byrow=TRUE)
    b1 = matrix( rep(prior.a1-1,nT)/(prior.fold*br), nG, nT, byrow=TRUE)

    lmm0 <- a0 * log( b0 ) +
                                        #    xx * log ( tns ) +
      lgamma ( xx + a0 ) -
        lgamma ( a0 ) -
          lgamma ( xx + 1 ) -
            ( xx + a0 ) * log ( b0 + tns )
    lm0 <- rowSums(lmm0)
    lmm1 <- a1 * log( b1 ) +
                                        #    xx * log ( tns ) +
      lgamma ( xx + a1 ) -
        lgamma ( a1 ) -
          lgamma ( xx + 1 ) -
            ( xx + a1 ) * log ( b1 + tns )
    lm1 <- rowSums(lmm1)

    post <- priorH0 * exp(lm0) / ( priorH0 * exp(lm0) + (1-priorH0) * exp(lm1) )
    lbf <- ( lm1 - lm0 ) / log(10)

    out[,c("PoissonlogBF","PoissonPosterior","Poissonlmlik0","Poissonlmlik1")] <- cbind(lbf,post,lm0,lm1)

  }
  
  if (nrow(xx)>0) {
    pass <- rep(TRUE,nG)
    if (filter.above + filter.below > 0){
      NucsSeq <- rowSums(cma.data[,c(ColCovDiscIndel,ColCovValiIndel)])
      below <- GeneTotalSize < filter.threshold
      above <- !below
      if ( sum(below==TRUE) > 1 ) {
        pass[below==TRUE]  <- rowSums(xx[below==TRUE,]) / NucsSeq[below==TRUE] > filter.below / 10^6
      }
      if ( sum(below==TRUE) ==1 ) {
        pass[below==TRUE]  <- sum( xx[below==TRUE,] ) / NucsSeq[below==TRUE] > filter.below / 10^6
      }
      if ( sum(above==TRUE) >1 ) {
        pass[above==TRUE]  <- rowSums(xx[above==TRUE,]) / NucsSeq[above==TRUE] > filter.above / 10^6
      }
      if ( sum(above==TRUE) ==1 ) {
        pass[above==TRUE] <- sum( xx[above==TRUE,] ) / NucsSeq[above==TRUE] > filter.above / 10^6
      }
      out[!pass,"neglogPg"] <- 0
      out[!pass,"CaMP"] <- -Inf
      out[!pass,"logLRT"] <- 0
      out[!pass,"logitBinomialPosteriorDriver"] <- -Inf
      out[!pass,c("PoissonlogBF","PoissonPosterior","Poissonlmlik0","Poissonlmlik1")] <-
        c(1,priorH0,NA,NA)
    }
  }
  
  ##only return the scores which are selected
  out <- out[, scores, drop = FALSE]
  
  return(out)
}
