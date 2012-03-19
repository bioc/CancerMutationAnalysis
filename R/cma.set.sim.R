##do a bunch of simulations to get p-values for the various methods
##pass.null = if TRUE, use passenger null instead of permutation null
##return.sim.data = if TRUE, return actual simulated objects as well
cma.set.sim <- function(cma.alter,
                        cma.cov,
                        cma.samp,
                        GeneSets,
                        passenger.rates = 
                        t(data.frame(0.55*rep(1.0e-6,25))),
                        ID2name = NULL,
                        BH = TRUE,
                        nr.iter,
                        pass.null = FALSE,
                        perc.samples = NULL,
                        spiked.set.sizes = NULL,
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
                        show.iter = TRUE,
                        KnownMountains = c("EGFR","SMAD4","KRAS",
                          "TP53","CDKN2A","MYC","MYCN","PTEN","RB1"),
                        exclude.mountains=TRUE,
                        verbose = TRUE)
{
  ##message a random number
  ##message(paste("Random number:", runif(1)))
  
  ##make sure everything is a data frame
  cma.alter <- as.data.frame(cma.alter)
  cma.cov <- as.data.frame(cma.cov)
  cma.samp <- as.data.frame(cma.samp)
  
  ##first create coverage object
  Coverage <- make.cov.obj(cma.cov, cma.samp)
  
  ##now create mutations object
  MutPass <- make.mut.obj(cma.cov, Coverage, cma.alter, pass.rates = passenger.rates)
  Mutations <- MutPass$cma.data
  passenger.rates <- MutPass$passenger.rates
  
  ##create EventsBySample object
  ##first add in "missing" context
  cma.alter$WTNuc <- as.character(cma.alter$WTNuc)
  cma.alter$Context <- as.character(cma.alter$Context)
  cma.alter$MuTNuc <- as.character(cma.alter$MutNuc)

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
  
  ##EventsBySampleBySet <-
  ##    get.EventsBySampleBySet.from.EventsBySample(EventsBySample,
  ##                                                GeneSets,
  ##                                                GoodSamples)
  
  ##get simulated data
  sim.data <- sim.data(cma.cov = cma.cov,
                       cma.samp = cma.samp,
                       GeneSets = GeneSets,
                       genes = GeneSymb,
                       passenger.rates = passenger.rates,
                       nr.iterations = nr.iter,
                       EventsBySample = EventsBySample,
                       ##EventsBySampleBySet = EventsBySampleBySet,
                       GoodSamples = GoodSamples,
                       pass.null = pass.null,
                       perc.samples = perc.samples,
                       spiked.set.sizes = spiked.set.sizes,
                       KnownMountains = KnownMountains,
                       exclude.mountains = exclude.mountains,
                       verbose = verbose)
  
  ##all.genes <-
  ##  shorten.gene.names(unique(as.character(sim.data$cma.cov[[1]]$Gene)))
  
  ##get mutated spiked-in genes
  ##spiked.sets <- setdiff(names(sim.data$GeneSets),
  ##                       names(GeneSets))
  
  ##message("all altered genes")
  ##message(sort(rownames(sim.data$Mutations[[1]])))
  ##for(i in 1:nr.iter)
  ##{
  ##    spiked.genes.trans <- setdiff(rownames(sim.data$Coverage[[i]]),
  ##                                  rownames(Coverage))
  ##    spiked.genes.trans <-
  ##        spiked.genes.trans[grep("gene.25.25.", spiked.genes.trans)]
  
  ##    spiked.genes.mut <- shorten.gene.names(spiked.genes.trans)
  
  ##message(paste("Mutations: Iteration #", i, sep = ""))
  ##get total number of mutations and coverage of spiked-in genes from gene.set.25.25
  ##    mut.spiked.genes <-
  ##        cbind(rowSums(sim.data$Mutations[[i]][spiked.genes.trans, 1:25]),
  ##              sim.data$Mutations[[i]][spiked.genes.trans, 68])
  ##    rownames(mut.spiked.genes) <- spiked.genes.mut
  ##    colnames(mut.spiked.genes) <- c("Mutations", "Coverage")
  ##message them
  ##message(mut.spiked.genes)
  
  ##message("Samples in which these genes are altered")
  ##message(sim.data$EventsBySample[[i]][sim.data$EventsBySample[[i]]$Symbol %in%
  ##                                   spiked.genes.mut,])
  ##}
  
  ##get spiked-in genes
  ##spiked.genes <- setdiff(rownames(sim.data$Coverage[[1]]),
  ##                        rownames(Coverage))
  
  results.sim <- list()
  for(i in 1:nr.iter)
    {
      message("")
      
      if(show.iter)
        {
          message(paste("Implement methods: Iteration #", i, sep = ""))
        }
      
      results.sim[[i]] <-
        cma.set.stat(cma.alter = 
                     sim.data$cma.alter[[i]],
                     cma.cov =
                     sim.data$cma.cov[[i]],
                     cma.samp =
                     sim.data$cma.samp[[i]],
                     Scores =
                     sim.data$Scores[[i]],
                     GeneSets =
                     sim.data$GeneSets,
                     passenger.rates = passenger.rates,
                     ID2name = NULL,
                     BH = BH,
                     gene.method = gene.method,
                     perm.null.method = perm.null.method,
                     perm.null.het.method = perm.null.het.method,
                     pass.null.method = pass.null.method,
                     pass.null.het.method = pass.null.het.method)
    }
  
  ##message ranks of spiked-in genes for basic and driver methods
  ##for(i in 1:nr.iter)
  ##{
  ##    message(paste("Iteration #", i, sep = ""))
  ##    ranks.basic <- rank(results.sim[[i]]$p.values.subject, ties.method = "min")
  ##    ranks.driver <- rank(results.sim[[i]]$p.values.drivers, ties.method = "min")
  ##    ranks.driver.pass <- rank(results.sim[[i]]$p.values.drivers.pass, ties.method = "min")
  
  ##    names(ranks.basic) <- names(ranks.driver) <- names(ranks.driver.pass) <-
  ##        rownames(results.sim[[i]])
  ##    message("ranks for basic method!")
  ##    message(ranks.basic[spiked.sets])
  ##    message("ranks for driver method!")
  ##    message(ranks.driver[spiked.sets])
  ##    message("ranks for driver method - pass!")
  ##    message(ranks.driver.pass[spiked.sets])
  ##}
  
  return.stuffs <- list(sim.data = sim.data,
                        results.sim = results.sim)
  
  return.obj <- new("SetMethodsSims")
  if(pass.null)
    {
      return.obj@null.dist <- "Passenger null"
    }
  else
    {
      return.obj@null.dist <- "Permutation null"
    }
  return.obj@perc.samples <- perc.samples
  return.obj@spiked.set.sizes <- spiked.set.sizes
  return.obj@GeneSets <- return.stuffs$sim.data$GeneSets
  return.obj@cma.alter <- return.stuffs$sim.data$cma.alter
  return.obj@cma.cov <- return.stuffs$sim.data$cma.cov
  return.obj@cma.samp <- return.stuffs$sim.data$cma.samp
  return.obj@Scores <- return.stuffs$sim.data$Scores
  return.obj@results <- return.stuffs$results.sim

  return(return.obj)
  
}
