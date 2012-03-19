cma.set.stat <-
    function(cma.alter,
             cma.cov,
             cma.samp,
             GeneSets,
             ID2name=NULL,
             Scores,
             ##passenger rates
             passenger.rates = t(data.frame(0.55*rep(1.0e-6,25))),
             BH = TRUE,
             ##T/F whether to do gene method
             gene.method = FALSE,
             ##T/F whether to do subject method
             perm.null.method = TRUE,
             ##T/F whether to do heterogeneous subject method
             perm.null.het.method = FALSE,
             ##T/F whether to do passenger null method
             pass.null.method = FALSE,
             ##T/F whether to do heterogeneous passenger null method
             pass.null.het.method = FALSE,
             score = "logLRT",
             verbose = TRUE)
{
  ##make sure everything is a data frame
  cma.alter <- as.data.frame(cma.alter)
  cma.cov <- as.data.frame(cma.cov)
  cma.samp <- as.data.frame(cma.samp)

  ##remove all data from prevalence samples (since analysis only works for discovery samples)
  cma.samp$NrSamp[cma.samp$Screen == "Prev"] <- 0
  cma.cov$Coverage[cma.cov$Screen == "Prev"] <- 0
  cma.alter <- cma.alter[cma.alter$Screen == "Disc", ]
  
  ##first create coverage object
  Coverage <- make.cov.obj(cma.cov, cma.samp)
        
  ##now create mutations object
  MutPass <- make.mut.obj(cma.cov, Coverage, cma.alter, pass.rates = passenger.rates)
  cma.data <- MutPass$cma.data
  passenger.rates <- MutPass$passenger.rates

  ##create EventsBySample object
  ##first add in "missing" context
  cma.alter$WTNuc <- as.character(cma.alter$WTNuc)
  cma.alter$Context <- as.character(cma.alter$Context)
  cma.alter$MutNuc <- as.character(cma.alter$MutNuc)

  cma.alter$Context <- paste("in", cma.alter$Context, "to", sep=".")
  cma.alter$Context[cma.alter$WTNuc == "C" & cma.alter$Context == "in..to"] <-
    "not.in.CpG.or.TpC.to"
  cma.alter$Context[cma.alter$WTNuc == "G" & cma.alter$Context == "in..to"] <-
    "not.in.CpG.or.GpA.to"
  cma.alter$Context[cma.alter$WTNuc == "A" & cma.alter$Context == "in..to"] <-
    "to"
  cma.alter$Context[cma.alter$WTNuc == "T" & cma.alter$Context == "in..to"] <-
    "to"
  EventsBySample <- data.frame(Event = cma.alter$Type,
                               Sample = cma.alter$Sample,
                               Symbol = cma.alter$Gene,
                               MutationClass = paste(cma.alter$WTNuc,
                                 cma.alter$Context, cma.alter$MutNuc, sep = "."))
  EventsBySample$MutationClass <- as.character(EventsBySample$MutationClass)
  EventsBySample$MutationClass[EventsBySample$MutationClass == ".in.All.to.ins.del"] <-
    "ins.del"
  EventsBySample$MutationClass[EventsBySample$MutationClass == ".in..to."] <-
    NA
  EventsBySample$Sample <- as.character(EventsBySample$Sample)

  ##consider only point mutations
  EventsBySample <-
    EventsBySample[EventsBySample$Event == "Mut",]
    
  ##transform the GeneSets object into a list, if it is not already a list
  GeneSets <- as.list(GeneSets)
  
  ##get the unique samples that have at least one mutation
  GoodSamples <- unique(EventsBySample[,"Sample"])
  
  ##get the gene symbols
  GeneNames <- rownames(Coverage)
  GeneSymb <- shorten.gene.names(GeneNames)
  names(GeneSymb) <- GeneNames
  
  ##use the ID2name if it is supplied
  IDs <- unique(unlist(GeneSets))
  if(!is.null(ID2name))
    {
      GeneSets <- lapply(GeneSets,
                         function(set, i2n) {i2n[set]},
                         ID2name)
    }
  ##now eliminate all the NAs
  GeneSets <- lapply(GeneSets, function(x) {x[!is.na(x)]})

  ##store this version of the gene-sets to use for the permutation null methods
  GeneSetsSymbs <- GeneSets
  
  ##now change GeneSets to have transcripts, not genes
  ##function to get transcripts for a vector of genes
  ##genes is a vector of genes
  ##transc2genes is a "look-up" vector which has the transcripts as names and the genes as transcripts
  get.trans.from.genes <- function(genes, transc2genes)
    {
      names(transc2genes[which(transc2genes %in% genes)])
    }

  GeneSets <- lapply(GeneSets, get.trans.from.genes, GeneSymb)
  
  EventsBySampleBySet <-
    get.EventsBySampleBySet.from.EventsBySample(EventsBySample,
                                                GeneSets,
                                                GoodSamples)
  
  ##get number of events in each sample
  nr.events.sample <-
    sapply(split(EventsBySample$Sample,
                 EventsBySample$Sample),
           length)
  nr.events.sample <- nr.events.sample[GoodSamples]
  
  ##change the gene-set annotations to have transcripts, not genes
  ##gs.temp <- GeneSets
  ##for(set in names(GeneSets))
  ##  {
  ##    ##get all the entries in GeneNames corresponding
  ##    ##to the entries in GeneSymb
  ##    GeneSets[[set]] <- GeneNames[GeneSymb %in% GeneSets[[set]]]
  ##  }
  
  ##gene-oriented method
  if(gene.method)
    {
      message("Gene method!")
      message(date())
      
      ##GeneScores <- data.frame(matrix(0,length(unique(GeneSymb)),1))
      ##rownames(GeneScores) <- unique(GeneSymb)
      ##colnames(GeneScores) <- score
      
      GeneScores <- data.frame(matrix(0,length(GeneNames),1))
      rownames(GeneScores) <- GeneNames
      colnames(GeneScores) <- score
      
      ##transcripts2genes <- shorten.gene.names(rownames(Scores))
      ##names(transcripts2genes) <- rownames(Scores)
      ##get only the first transcript if there are multiple transcripts per gene
      ##Scores.short <- rep(0, length(unique(transcripts2genes)))
      ##names(Scores.short) <- unique(transcripts2genes)
      ##for(g in names(Scores.short))
      ##{
      ##    t <- names(transcripts2genes[transcripts2genes == g])
      ##    Scores.short[g] <- Scores[t[1],score]
      ##}
      
      ##message(head(Scores.short))
      ##message(head(GeneScores))
      
      GeneScores[rownames(Scores),] <- Scores[,score, drop = FALSE]
      
      ## put non-zero LRT scores in GeneScores data frame
      ##GeneScores[match(names(Scores.short),
      ##                 rownames(GeneScores)),score] <-
      ##                     Scores.short
      
      ##get p-values using the GSEA approach
      ##GeneScores.logLRT <- GeneScores[,score, drop = FALSE]
      
      p.values.gene <- sapply(GeneSets,runGSEA,
                              data=GeneScores,
                              datacol=1,
                              IDinRowNames=TRUE, absolute=TRUE,
                              type= "f", alternative="mixed",
                              ranks.only=TRUE, nsim=1)
      
      q.values.gene <- p2q(p.values.gene, BH = BH)
    }
  else
    {
      p.values.gene <- NA
      q.values.gene <- NA
    }
  
  ##now switch GeneSets object back to having just the genes, not the transcripts
  ##GeneSets <- gs.temp
  
  ##get p-values using the subject approach
  if(perm.null.method)
    {
      message("Permutation null w/o heterogeneity")
      message(date())
      
      gs <- GeneSymb
      names(gs) <- NULL
      
      p.values.subject <-
        perm.null.p.values(EventsBySampleBySet =
                           EventsBySampleBySet,
                           GoodSamples = GoodSamples,
                           GeneSymb = gs,
                           nr.events.sample =  nr.events.sample,
                           GeneSets = GeneSetsSymbs)

      q.values.subject <- p2q(p.values.subject, BH = BH)
    }
  else
    {
      p.values.subject <- NA
      q.values.subject <- NA
    }
  
  ##get p-values using the subject "heterogeneous" approach
  if(perm.null.het.method)
    {
      message("Permutation null w/ heterogeneity")
      message(date())
      
      ##dyn.load("get.cdf.het.berns.so")

      gs <- GeneSymb
      names(gs) <- NULL
      
      p.values.subject.het <-
        perm.null.het.p.values(EventsBySampleBySet =
                               EventsBySampleBySet,
                               GoodSamples = GoodSamples,
                               nr.events.sample =  nr.events.sample,
                               GeneSets = GeneSetsSymbs,
                               GeneSymb = gs)

      q.values.subject.het <- p2q(p.values.subject.het, BH = BH)
    }
  else
    {
      p.values.subject.het <- NA
      q.values.subject.het <- NA
    }
  
  ##passenger null method
  if(pass.null.method)
    {
      message("Passenger null w/o heterogeneity")
      message(date())
      
      p.values.pass.null <-
        pass.null.p.values(GeneSets = GeneSets,
                           GoodSamples = GoodSamples,
                           GeneSymb = GeneSymb,
                           nr.events.sample = nr.events.sample,
                           passenger.rates = passenger.rates,
                           Coverage = Coverage,
                           EventsBySample = EventsBySample,
                           EventsBySampleBySet = EventsBySampleBySet)

      q.values.pass.null <- p2q(p.values.pass.null, BH = BH)
    }
  else
    {
      p.values.pass.null <- NA
      q.values.pass.null <- NA
    }
  
  
  ##heterogeneous passenger null method
  if(pass.null.het.method)
    {
      message("Passenger null w/ heterogeneity")
      message(date())
      
      ##dyn.load("get.cdf.het.berns.so")
      
      p.values.pass.null.het <-
        pass.null.het.p.values(GeneSets = GeneSets,
                               GoodSamples = GoodSamples,
                               GeneSymb = GeneSymb,
                               nr.events.sample = nr.events.sample,
                               passenger.rates = passenger.rates,
                               Coverage = Coverage,
                               EventsBySample = EventsBySample,
                               EventsBySampleBySet = EventsBySampleBySet)

      q.values.pass.null.het <- p2q(p.values.pass.null.het, BH = BH)
    }
  else
    {
      p.values.pass.null.het <- NA
      q.values.pass.null.het <- NA
    }
  
  message(date())
  
  results <- data.frame(p.values.gene = p.values.gene,
                        q.values.gene = q.values.gene,
                        p.values.perm.null = p.values.subject,
                        q.values.perm.null = q.values.subject,
                        p.values.perm.null.het = p.values.subject.het,
                        q.values.perm.null.het = q.values.subject.het, 
                        p.values.pass.null = p.values.pass.null,
                        q.values.pass.null = q.values.pass.null, 
                        p.values.pass.null.het = p.values.pass.null.het,
                        q.values.pass.null.het = q.values.pass.null.het)
  
  rownames(results) <- names(GeneSets)

  ##only return results for methods which were implemented
  results <- results[,!is.na(colSums(results))]
  
  results
}

