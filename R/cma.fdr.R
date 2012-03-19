cma.fdr <- function(cma.alter,
                    cma.cov,
                    cma.samp,
                    scores = c("CaMP", "logLRT"),
                    passenger.rates = t(data.frame(.55*rep(1.0e-6,25))),
                    allgenes=TRUE,
                    estimate.p0=FALSE,
                    p0.step=1,
                    p0=1,
                    eliminate.noval=FALSE,
                    filter.threshold=0, 
                    filter.above=0,
                    filter.below=0,
                    filter.mutations=0, 
                    aa=1e-10, 
                    bb=1e-10,
                    priorH0=1-500/13020, 
                    prior.a0=100,
                    prior.a1=5,
                    prior.fold=10,
                    M=2,
                    DiscOnly=FALSE,
                    PrevSamp="Sjoeblom06",
                    KnownCANGenes=NULL,
                    showFigure=FALSE,
                    cutoffFdr=0.1){
  
  ##if DiscOnly is TRUE, then only consider discovery samples
  if(DiscOnly)
    {
      ##remove all data from prevalence samples
      cma.samp$NrSamp[cma.samp$Screen == "Prev"] <- 0
      cma.cov$Coverage[cma.cov$Screen == "Prev"] <- 0
      cma.alter <- cma.alter[cma.alter$Screen == "Disc", ]
    }
  
  ##make sure everything is a data frame
  cma.alter <- as.data.frame(cma.alter)
  cma.cov <- as.data.frame(cma.cov)
  cma.samp <- as.data.frame(cma.samp)

  ##now keep only the mutations in cma.alter
  cma.alter$Type <- as.character(cma.alter$Type)
  cma.alter <- cma.alter[cma.alter$Type=="Mut", ]
    
  ##first create coverage object
  coverage <- make.cov.obj(cma.cov, cma.samp)
  ##now create mutations object
  MutPass <- make.mut.obj(cma.cov, coverage, cma.alter, pass.rates = passenger.rates)
  cma.data <- MutPass$cma.data
  passenger.rates <- MutPass$passenger.rates

  ##get number of genes
  number.genes <- nrow(coverage)
  
  if (estimate.p0) { allgenes==TRUE; eliminate.noval=FALSE }
  if (allgenes) {
    cma.data <- add.zeroes.cma.data(cma.alter = cma.alter,
                                    cma.cov = cma.cov,
                                    cma.samp = cma.samp)
    filter.above=0
    filter.below=0
    eliminate.noval==FALSE
   }

  # observed statistics
  obs.score <-
    cma.scores(cma.data = cma.data,
               coverage = coverage,
               scores = scores,
               passenger.rates = passenger.rates,
               ##compute.poisson.BF=FALSE,
               ##compute.binomial.posterior=FALSE,
               filter.above=filter.above,
               filter.below=filter.below,
               filter.threshold=filter.threshold,
               aa=aa, 
               bb=bb,
               priorH0=priorH0, 
               prior.a0=prior.a0,
               prior.a1=prior.a1,
               ##prior.fold=prior.fold)[,c("CaMP",scores)]
               prior.fold=prior.fold)[,scores,drop=FALSE]
  
  ColMutDisc <- grep("MutationsDiscovery",colnames(cma.data))
  ColMutVali <- grep("MutationsValidation",colnames(cma.data))

  pass <- rep(TRUE,nrow(obs.score))
  if ( eliminate.noval ) {
    DiscoveryMutations <- apply(cma.data[,ColMutDisc],1,sum)
    ValidationMutations <- apply(cma.data[,ColMutVali],1,sum)
    pass <- DiscoveryMutations > 0 & ValidationMutations > 0
  }
  if (filter.mutations>0) {
    pass <- pass & ( apply(cma.data[,c(ColMutDisc,ColMutVali)],1,sum) >= filter.mutations )
  }

  obs.score <- obs.score[pass,,drop=FALSE]
  
  # null distribution
  null.score <- data.frame()
  for (m in 1:M){
    null <- cma.simulator(cma.cov=cma.cov,
                          cma.samp=cma.samp,
                          passenger.rates,
                          eliminate.noval=eliminate.noval,
                          filter.mutations=filter.mutations,
                          PrevSamp=PrevSamp,
                          KnownCANGenes=KnownCANGenes)

    ms <- cma.scores(cma.data=null,
                     coverage=coverage,
                     scores = scores,
                     passenger.rates = passenger.rates,
                     ##compute.poisson.BF=FALSE,
                     ##compute.binomial.posterior=FALSE,
                     filter.above=filter.above,
                     filter.below=filter.below,
                     filter.threshold=filter.threshold,
                     aa=aa, 
                     bb=bb,
                     priorH0=priorH0, 
                     prior.a0=prior.a0,
                     prior.a1=prior.a1,
                     prior.fold=prior.fold
                     )[,scores,drop=FALSE]
    
    null.score <- rbind(null.score,ms) 
  }
  
  if (nrow(null.score)==0) stop("no null genes made it through to be validated; increase M")
  out <- list()

  ##browser()
  
  for (ss in 1:length(scores) ){
    ##obs <- obs.score[,ss+1]; names(obs) <- rownames(obs.score)
    
    obs <- obs.score[,ss]
    names(obs) <- rownames(obs.score)
    
    obs[obs == Inf] = 10^6
    
    tmp <- ebfdr(obs.score=obs,
                 null.score=null.score[,ss],
                 ##null.score=null.score,
                 mass=M,
                 estimate.p0 = estimate.p0,
                 p0.step = p0.step,
                 p0=p0,
                 showFigure=showFigure,
                 cutoffFdr=cutoffFdr,
                 score.name=scores[ss])
    ##tmp <- cbind( obs.score[rownames(tmp),"CaMP"], tmp )
    ##colnames(tmp)[1] <- "CaMP"

    if(ss<length(scores) & showFigure)
      {
        readline(paste("Score: ",scores[ss],". Hit return for next score.\n",
                       sep=""))
      }
    if(ss==length(scores) & showFigure)
      {
        cat(paste("Score: ",scores[ss],".",sep=""))
      }
    
    out[[ss]] <- tmp
  }
  names(out) <- scores
  
  return(out)
}
