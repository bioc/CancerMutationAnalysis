make.mut.obj <- function(cma.cov, Cov, cma.alter, pass.rates)
  {
    ##get rid of genes named "NULL"
    cma.alter <- cma.alter[cma.alter$Gene != "NULL",]

    ##create mutations object
    tempCols <- as.character(cma.cov$WTNuc[(1:9)*2-1])
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
    tempCols <- tempCols[tempCols != "ins.del.in."]
    tempCols <- paste(rep(tempCols, each=3), "to",
                      c("G","A","T",
                        "C","A","T",
                        "C","A","T",
                        "G","A","T",
                        "C","G","T",
                        "G","A","T",
                        "C","A","T",
                        "C","G","A"),sep=".")
    tempCols <- c(tempCols, "ins.del")
    tempCols <- c(paste("MutationsDiscovery",tempCols,sep="."),
                  paste("MutationsValidation",tempCols,sep="."))
    tempCols <- c(tempCols, colnames(Cov))

    cma.data <- matrix(0, nrow =
                       length(unique(cma.alter$Gene[cma.alter$Type=="Mut"])),
                       ncol = length(tempCols))
    colnames(cma.data) <- tempCols
    rownames(cma.data) <-
      as.character(sort(unique(cma.alter$Gene[cma.alter$Type=="Mut"])))
    cma.data[,grep("Coverage",colnames(cma.data))] <-
      as.matrix(Cov[rownames(cma.data),grep("Coverage",colnames(Cov))])
    cma.data[,grep("NumberTumors",colnames(cma.data))] <-
      as.matrix(Cov[rownames(cma.data),grep("NumberTumors",colnames(Cov))])
    ##now fill in the mutations
    ##first create the mutational types from the Alter object
    JustMuts <- cma.alter[cma.alter$Type=="Mut",]
    MutTypes <- apply(JustMuts[,c("WTNuc","Context")], 1,
                      paste, sep="", collapse=".in.")
    MutTypes <- paste(MutTypes, JustMuts[,"MutNuc"], sep=".to.")
    MutTypes[grep("All",MutTypes)] <- "ins.del"
    MutTypes[grep("G.in..to.",MutTypes)] <-
      gsub("G.in..to.", "G.not.in.CpG.or.GpA.to.",
           MutTypes[grep("G.in..to.",MutTypes)])
    MutTypes[grep("C.in..to.",MutTypes)] <-
      gsub("C.in..to.", "C.not.in.CpG.or.TpC.to.",
           MutTypes[grep("C.in..to.",MutTypes)])
    MutTypes[grep("A.in..to.",MutTypes)] <-
    gsub("A.in..to.", "A.to.",
         MutTypes[grep("A.in..to.",MutTypes)])
    MutTypes[grep("T.in..to.",MutTypes)] <-
      gsub("T.in..to.", "T.to.",
         MutTypes[grep("T.in..to.",MutTypes)])
    ##now add whether it is discovery or validation
    tempScreen <- as.character(JustMuts$Screen)
    tempScreen[tempScreen=="Disc"] <- "Discovery"
    tempScreen[tempScreen=="Prev"] <- "Validation"
    MutTypes <- paste("Mutations",tempScreen,".",MutTypes,sep="")
    ##now add the gene-specific mutations to the Mutations object
    for(gene in rownames(cma.data))
      {
        ##get the mutations for that specific genes
        nn <- table(MutTypes[which(JustMuts$Gene == gene)])
        cma.data[gene,names(nn)] <- nn
    }
    cma.data <- as.data.frame(cma.data)

    ##create column names based just on nucleotide contexts
    contexts <- c("C.in.CpG", "G.in.CpG", "G.in.GpA", "C.in.TpC", "A",
                  "C.not.in.CpG.or.TpC", "G.not.in.CpG.or.GpA", "T",
                  "ins.del")
    MutTypes.contexts <- c()
    for(c in contexts)
      {
        ##get nucleotide
        nucl <- strsplit(c, "\\.")[[1]][1]
        ##make mutation types
        if(nucl != "ins")
          {
            for(m in c("C","G","A","T"))
              {
                if(nucl != m)
                  {
                    MutTypes.contexts <-
                      c(MutTypes.contexts,
                        paste(c,"to",m,sep="."))
                  }
              }
          }
        else
          {
            MutTypes.contexts <-
              c(MutTypes.contexts,
                c)
          }
      }
    colsDisc.contexts <- paste("MutationsDiscovery.", MutTypes.contexts, sep = "")
    colsVal.contexts <- paste("MutationsValidation.", MutTypes.contexts, sep = "")
    colsCovDisc.contexts <- paste("CoverageAtRiskDiscovery.", contexts, sep = "")
    colsCovVal.contexts <- paste("CoverageAtRiskValidation.", contexts, sep = "")
    colnames.contexts <- c(colsDisc.contexts,colsVal.contexts,
                           colsCovDisc.contexts,colsCovVal.contexts,
                           "NumberTumorsAnalyzedDiscovery",
                           "NumberTumorsAnalyzedValidation")
    if(!identical(colnames(cma.data),colnames.contexts))
      {
        stop("There's a problem with the mutations object - check the colnames!")
    }

  ##make sure the pass.rates object has the proper names
  names(pass.rates) <-
    c("C.in.CpG.to.G", "C.in.CpG.to.A", "C.in.CpG.to.T",
      "G.in.CpG.to.C", "G.in.CpG.to.A", "G.in.CpG.to.T",
      "G.in.GpA.to.C", "G.in.GpA.to.A", "G.in.GpA.to.T",
      "C.in.TpC.to.G", "C.in.TpC.to.A", "C.in.TpC.to.T",
      "A.to.C", "A.to.G", "A.to.T",
      "C.not.in.CpG.or.TpC.to.G", "C.not.in.CpG.or.TpC.to.A",
      "C.not.in.CpG.or.TpC.to.T",
      "G.not.in.CpG.or.GpA.to.C", "G.not.in.CpG.or.GpA.to.A",
      "G.not.in.CpG.or.GpA.to.T",
      "T.to.C", "T.to.G", "T.to.A",
      "ins.del")

  if(!identical(names(pass.rates),MutTypes.contexts))
    {
      stop("There's a problem with the rates object - check the names!")
    }

    return(list(cma.data = cma.data,
                passenger.rates = pass.rates))
  }
