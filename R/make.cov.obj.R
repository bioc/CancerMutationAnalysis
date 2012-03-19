make.cov.obj <- function(cma.cov, cma.samp)
  {
    ##first create column names
    tempCols <- as.character(cma.cov$WtNuc[(1:9)*2-1])
    ##label indel column
    tempCols[tempCols == ""] <- "ins.del"
    tempCols2 <- as.character(cma.cov$Context[(1:9)*2-1])
    tempCols2[tempCols2=="All"] <- ""
    tempCols2 <- paste("in",tempCols2,sep=".")
    tempCols <- paste(tempCols, tempCols2, sep=".")
    tempCols[tempCols == "A.in."] <- "A"
    tempCols[tempCols == "T.in."] <- "T"
    tempCols[tempCols == "C.in."] <- "C.not.in.CpG.or.TpC"
    tempCols[tempCols == "G.in."] <- "G.not.in.CpG.or.GpA"
    tempCols[tempCols == "ins.del.in."] <- "ins.del"
    tempCols <- c(paste("CoverageAtRiskDiscovery",tempCols,sep="."),
                  paste("CoverageAtRiskValidation",tempCols,sep="."))
    tempCols <- c(tempCols, c("NumberTumorsAnalyzedDiscovery",
                              "NumberTumorsAnalyzedValidation"))
    
    ##check if Gene is in the same order in cma.samp and in cma.cov
    if(!identical(unique(cma.cov$Gene), unique(cma.samp$Gene)))
      {
        stop("Check that the `Gene` component of the GeneCov and GeneSamp objects is the same")
      }
    
    tempCov <- cbind(matrix(cma.cov$Coverage[cma.cov$Screen=="Disc"],
                            nrow = length(unique(cma.cov$Gene)),
                            ncol = 9, byrow = TRUE),
                     matrix(cma.cov$Coverage[cma.cov$Screen=="Prev"],
                            nrow = length(unique(cma.cov$Gene)),
                            ncol = 9, byrow = TRUE),
                     cma.samp$NrSamp[cma.samp$Screen=="Disc"],
                     cma.samp$NrSamp[cma.samp$Screen=="Prev"])
    rownames(tempCov) <- unique(cma.cov$Gene)
    colnames(tempCov) <- tempCols
    tempCov <- as.data.frame(tempCov)
    Cov <- tempCov
        
    Cov
  }
