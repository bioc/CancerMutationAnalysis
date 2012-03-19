cma.simulator <-
  function(cma.cov,
           cma.samp,
           passenger.rates = t(data.frame(0.55*rep(1.0e-6,25))),
           multinomial=FALSE,
           simulation.type="rates",
           eliminate.noval=TRUE,
           filter.mutations=0, # restrict output to at least filter.mutations in total
           random.passenger.rates=FALSE,
           lognormal.sigma=.58,
           null.genome=TRUE,
           ncangenes=NULL,
           random.ncangenes=FALSE,
           driver.rates=NULL,
           driver.prevalences=NULL,
           PrevSamp="Sjoeblom06",
           KnownCANGenes=NULL
           ){

    ##make sure everything is a data frame
    cma.cov <- as.data.frame(cma.cov)
    cma.samp <- as.data.frame(cma.samp)
    
    ##create coverage object
    coverage <- make.cov.obj(cma.cov, cma.samp)
    
    ColCovDisc <- grep("CoverageAtRiskDiscovery",colnames(coverage))[1:8]
    ColCovVali <- grep("CoverageAtRiskValidation",colnames(coverage))[1:8]
    ColCovDiscIndel <- grep("CoverageAtRiskDiscovery.ins.del",colnames(coverage))
    ColCovValiIndel <- grep("CoverageAtRiskValidation.ins.del",colnames(coverage))
    nT <- length(ColCovDisc)*3+1
    nGG <- sum(cma.samp$Screen == "Disc")
    gene.names <- as.character(cma.samp$Gene[cma.samp$Screen == "Disc"])
    
    nSd <- coverage$NumberTumorsAnalyzedDiscovery
    nSv <- coverage$NumberTumorsAnalyzedValidation
    ## for the validation, if the gene was validated, the actual sample is used,
    ##if not, the sample sizes from the validation phase are randomly imputed
    if (sum(nSv==0)==0) {
      nSv[nSv==0] <- sample( nSv[nSv>0], size=sum(nSv==0), replace=TRUE)
    }
    nS <- nSv + nSd
    
    ## coverage
    
    t1 <- coverage[,ColCovDisc]
    t2 <- matrix( rep(data.matrix(t1),3), nrow(t1), 3*ncol(t1) )[,c(8*(0:2)+1,8*(0:2)+2,8*(0:2)+3,8*(0:2)+4,8*(0:2)+5,8*(0:2)+6,8*(0:2)+7,8*(0:2)+8)]
    tnsD <- round ( cbind(t2,coverage[,ColCovDiscIndel]) )
    
    t1 <- coverage[,ColCovVali]
    t2 <- matrix( rep(data.matrix(t1),3), nrow(t1), 3*ncol(t1) )[,c(8*(0:2)+1,8*(0:2)+2,8*(0:2)+3,8*(0:2)+4,8*(0:2)+5,8*(0:2)+6,8*(0:2)+7,8*(0:2)+8)]
    tnsV <- round ( cbind(t2,coverage[,ColCovValiIndel]) )
    tns <- tnsD + tnsV
    
    tnsDsum <- rowSums(tnsD)
    tnsVsum <- rowSums(tnsV)
    scaleDV <- ifelse( is.finite( mean(tnsVsum[tnsVsum>0] / tnsDsum[tnsVsum>0] ) ),
                      mean(tnsVsum[tnsVsum>0] / tnsDsum[tnsVsum>0]) , 0)
    tnsV[tnsVsum==0,] <- round (scaleDV * tnsD[tnsVsum==0,] )
    tns <- tnsD + tnsV
    
    ## passenger rates matrix
    br <- data.matrix( passenger.rates )
    rate.multiplier <- rate.multiplier.D <- rate.multiplier.V <-
      matrix ( 1, nGG, nT )
    if (random.passenger.rates) {
      rate.multiplier <- rate.multiplier.D <- rate.multiplier.V <-
        exp( rnorm(nGG,-(lognormal.sigma^2)/2,lognormal.sigma) )
      rate.multiplier <- rate.multiplier.D <- rate.multiplier.V <-
        matrix ( rate.multiplier, nGG, nT )
    }
    
    ## spike-ins
    if (null.genome==FALSE){
      if (random.ncangenes) ncangenes <- max(1,rbinom(1,nGG,ncangenes/nGG))
      ii.can <- sample(1:nGG,ncangenes)
      if (simulation.type=="prevalence"){
        yy <- rbinom(ncangenes,nS,driver.prevalences)
        driver.rates <- matrix( yy / tns[ii.can,nT], length(yy), nT)
      }
      ii <- sample(1:nrow(driver.rates),ncangenes,replace=TRUE)
      driver.rates <- driver.rates[ii,]
      if ( nrow(br)==1 ) {
        driver.rates <- as.matrix(driver.rates)
        driver.multiplier <- driver.rates / matrix(as.numeric(br),ncangenes,nT,byrow=TRUE)
        rate.multiplier[ii.can,] <- driver.multiplier
        rate.multiplier[rate.multiplier==0] <- 1
      }
      if ( nrow(br)==2 ) {
        driver.multiplier <- driver.rates / matrix(as.numeric(br[1,]),ncangenes,nT,byrow=TRUE)
        rate.multiplier.D[ii.can,] <- driver.multiplier
        rate.multiplier.D[rate.multiplier.D==0] <- 1
        driver.multiplier <- driver.rates / matrix(as.numeric(br[2,]),ncangenes,nT,byrow=TRUE)
        rate.multiplier.V[ii.can,] <- driver.multiplier
        rate.multiplier.V[rate.multiplier.V==0] <- 1
      }
      gene.truth <- rep(FALSE,nGG); gene.truth[ii.can] <- TRUE
    }

    if ( nrow(br)==1 ) {
      rrD <- rrV <- rate.multiplier * matrix(as.numeric(br),nGG,nT,byrow=TRUE)
    } else {
      if ( nrow(br)==2 ) {
        rrD <- rate.multiplier * matrix(as.numeric(br[1,]),nGG,nT,byrow=TRUE)
        rrV <- rate.multiplier * matrix(as.numeric(br[2,]),nGG,nT,byrow=TRUE)
      }
      else {
        print("Rate object does not have 1 or 2 rows")
        stop()
      }
    }
    xxD <- matrix(rbinom(nGG*nT,data.matrix(tnsD),rrD),nGG,nT)
    xxD[is.nan(xxD)] <- max(xxD[!is.nan(xxD)])
    
    xxV <- matrix(rbinom(nGG*nT,data.matrix(tnsV),rrV),nGG,nT)
    if( sum( is.nan(xxV) ) == nGG*nT )
      {
        xxV[is.nan(xxV)] <- 0
      } else
    {
      xxV[is.nan(xxV)] <- max(xxV[!is.nan(xxV)])
    }
    
    ##add rownames here, so we can keep track of genes
    rownames(xxD) <- rownames(xxV) <- rownames(t1)
      
    tot.disc <- rowSums(xxD)
    
    if (multinomial){
      tnsDns <- round( coverage[,ColCovDisc] )
      tnsDindels <- round( coverage[,ColCovDiscIndel] )
      tnsVns <- round( coverage[,ColCovVali] )
      tnsVindels <- round( coverage[,ColCovValiIndel] )
      tnsDnssum <- rowSums(tnsDns)
      tnsVnssum <- rowSums(tnsVns)
      tnsVns[tnsVnssum==0,] <- round ( mean(tnsVnssum[tnsVnssum>0] /
               tnsDnssum[tnsVnssum>0]) * tnsDns[tnsVnssum==0,] )
      tnsVindels[tnsVindels==0] <- round ( mean(tnsVindels[tnsVindels>0] /
                   tnsDindels[tnsVindels>0]) * tnsDindels[tnsVindels==0] )
      
      for (gg in 1:nGG) {
        for (ii in 1:8) {
          proD <- c(rrD[gg,(1:3)+3*(ii-1)])
          proD <- c(1-sum(proD),proD)
          xxD[gg,(1:3)+3*(ii-1)] <- rmultinom(1, tnsDns[gg,ii], prob=proD)[-1]
        }
        xxD[gg,nT] <- rbinom(1,tnsDindels[gg],rrD[nT])
      }
      
      tot.disc <- rowSums(xxD)
      gg.disc <- (1:nGG)[tot.disc>0]
      for (gg in gg.disc) {
        for (ii in 1:8) {
          proV <- c(rrV[gg,(1:3)+3*(ii-1)])
          proV <- c(1-sum(proV),proV)
          xxV[gg,(1:3)+3*(ii-1)] <- rmultinom(1, tnsVns[gg,ii], prob=proV)[-1]
        }
        xxV[gg,nT] <- rbinom(1,tnsVindels[gg],rrV[nT])
      }
    }

    ##remove all the stuff that doesn't conform to the rules, depending on what genes are analyzed in the prevalence samples
    if(PrevSamp=="Sjoeblom06")
      {
        xxV[tot.disc==0,] <- 0
        tnsV[tot.disc==0,] <- 0
        nSv[tot.disc==0] <- 0
      }
    if(PrevSamp=="Parsons11")
      {
        ##remove all the stuff that doesn't conform to the rules
        GenesInVali <- (tot.disc >= 2) | (tot.disc >= 1 & names(tot.disc) %in% KnownCANGenes)
        ##add in validation samples for the genes making it through validation
        nSv[GenesInVali] <- max(nSv)
        ##zero out everything that's not in the validation stage
        xxV[!GenesInVali,] <- 0
        tnsV[!GenesInVali,] <- 0
        nSv[!GenesInVali] <- 0
      }
        
    tot.vali <- rowSums(xxV)
    
    xx <- xxD + xxV
    makesit <- rowSums(xx) >= filter.mutations
    if (eliminate.noval) makesit <- makesit & tot.disc>0 & tot.vali>0
    
    tums <- cbind( nSd, nSv )
    cove <- cbind( tnsD[,c(0:7*3+1,25)], tnsV[,c(0:7*3+1,25)] )
    
    contexts <- c("C.in.CpG",
                  "G.in.CpG",
                  "G.in.GpA",
                  "C.in.TpC",
                  "A",
                  "C.not.in.CpG.or.TpC",
                  "G.not.in.CpG.or.GpA",
                  "T",
                  "ins.del")
    mut.contexts <- c("C.in.CpG.to.G",
                      "C.in.CpG.to.A",
                      "C.in.CpG.to.T",
                      "G.in.CpG.to.C",
                      "G.in.CpG.to.A",
                      "G.in.CpG.to.T",
                      "G.in.GpA.to.C",
                      "G.in.GpA.to.A",
                      "G.in.GpA.to.T",
                      "C.in.TpC.to.G",
                      "C.in.TpC.to.A",
                      "C.in.TpC.to.T",
                      "A.to.C",
                      "A.to.G",
                      "A.to.T",
                      "C.not.in.CpG.or.TpC.to.G",
                      "C.not.in.CpG.or.TpC.to.A",
                      "C.not.in.CpG.or.TpC.to.T",
                      "G.not.in.CpG.or.GpA.to.C",
                      "G.not.in.CpG.or.GpA.to.A",
                      "G.not.in.CpG.or.GpA.to.T",
                      "T.to.C",
                      "T.to.G",
                      "T.to.A",
                      "ins.del")
    colsDisc <- paste("MutationsDiscovery", mut.contexts, sep=".")
    colsVal <- paste("MutationsValidation", mut.contexts, sep=".")
    colsTums <- c("NumberTumorsAnalyzedDiscovery",
                  "NumberTumorsAnalyzedValidation")
    colsCoveDisc <- paste("CoverageAtRiskDiscovery", contexts, sep=".")
    colsCoveVali <- paste("CoverageAtRiskValidation", contexts, sep=".")
    columns <- c(colsDisc,
                 colsVal,
                 colsCoveDisc,
                 colsCoveVali,
                 colsTums)
    out <- data.frame(cbind(xxD,xxV,cove,tums))
    colnames(out) <- columns
    rownames(out) <- rownames(coverage)
    
    if (null.genome==FALSE){
      out <- cbind(out,gene.truth)
      colnames(out) <- c(columns,"TrueDriver")
    }
    
    out <- out[makesit==TRUE,]
    
    return(out)
  }
