add.zeroes.cma.data <- function(cma.alter,cma.cov,cma.samp,
                                passenger.rates =
                                t(data.frame(0.55*rep(1.0e-6,25)))) {
  
  ##make sure everything is a data frame
  cma.alter <- as.data.frame(cma.alter)
  cma.cov <- as.data.frame(cma.cov)
  cma.samp <- as.data.frame(cma.samp)
  
  ##first create coverage object
  coverage <- make.cov.obj(cma.cov, cma.samp)
  ##now create mutations object
  MutPass <- make.mut.obj(cma.cov, coverage, cma.alter, pass.rates = passenger.rates)
  cma.data <- MutPass$cma.data
  passenger.rates <- MutPass$passenger.rates

  ColMut <- grep("Mutations",colnames(cma.data))
  ColCov1 <- grep("Coverage",colnames(cma.data))
  ColCov2 <- grep("Coverage",colnames(coverage))

  zeroes <- ! ( rownames(coverage) %in% rownames(cma.data) )
  nm <- nrow(cma.data)
  nn <- sum(zeroes)
  cma.data[(nm+1):(nm+nn),ColMut] <- matrix(0,nn,50)
  rownames( cma.data )[(nm+1):(nm+nn)] <- rownames( coverage )[zeroes]
  cma.data[(nm+1):(nm+nn),"NumberTumorsAnalyzedDiscovery"] <-
    coverage[zeroes,"NumberTumorsAnalyzedDiscovery"]
  cma.data[(nm+1):(nm+nn),"NumberTumorsAnalyzedValidation"] <- rep(0,nn)
  cma.data[(nm+1):(nm+nn),ColCov1] <- coverage[zeroes,ColCov2]

  return(cma.data)
}
