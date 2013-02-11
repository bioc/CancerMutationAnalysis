##return objects of the form EventsBySampleBrain etc (actually a list
##of such objects, of length nr.iterations),
##as well an "augmented" GeneSets object, with the fake spiked-in sets
##with their fake genes
##nr.iterations = nr of iterations
##GoodSamples = see above
##perc.samples = percentages of samples in which the gene-sets are spiked in
##spiked.set.sizes = sizes of sets spiked in
##mutations.only = flag on whether should just use mutational event
##(if FALSE, also use amplifications and deletions)
##pass.null = if TRUE, use passenger null instead of permutation null
sim.data <- function(cma.cov,
                     cma.samp,
                     GeneSets,
                     genes,
                     passenger.rates,
                     nr.iterations = 100,
                     EventsBySample,
                     EventsBySampleBySet,
                     GoodSamples,
                     pass.null = FALSE,
                     perc.samples = c(25, 50, 75, 90),
                     spiked.set.sizes = c(20, 100, 500),
                     KnownMountains = c("EGFR","SMAD4","KRAS",
                     "TP53","CDKN2A","MYC","MYCN","PTEN","RB1"),
                     exclude.mountains=TRUE,
                     verbose = TRUE)
{
  cma.cov <- as.data.frame(cma.cov)
  cma.samp <- as.data.frame(cma.samp)

  ##first create coverage object
  Coverage <- make.cov.obj(cma.cov, cma.samp)

  ##get the gene symbols
  GeneNames <- rownames(Coverage)
  GeneSymb <- shorten.gene.names(GeneNames)

  EventsBySample <- EventsBySample[EventsBySample$Event=="Mut",]

  ##remove mutations in mountains
  if (exclude.mountains) {
    EventsBySample <-
      EventsBySample[!(EventsBySample[, "Symbol"] %in% KnownMountains), ]
  }

  ##create vectors of names of spiked sets and names of spiked genes ("fake sets" and "fake genes"),
  fake.sets <- c()
  fake.genes <- c()
  ##also create vector of names of spiked transcripts (simply add "@@00" after gene name)
  fake.transcripts <- c()
  if(length(spiked.set.sizes) > 0)
    {
      for(size in spiked.set.sizes)
        {
          for(perc in perc.samples)
            {
              ##create name for spiked set
              fake.set <- paste("gene.set", size, perc, sep = ".")
              ##add to vector of names of spiked-in sets
              fake.sets <- c(fake.sets,
                             fake.set)
              ##create genes inside spiked set and add the set to GeneSets object
              fake.genes.set <- paste("gene", size, perc, "nr", 1:size, sep = ".")
              ##add the genes to the fake.genes object
              fake.genes <- c(fake.genes, fake.genes.set)
              GeneSets[[fake.set]] <- fake.genes.set
            }
        }
      fake.transcripts <- paste(fake.genes, "@@00", sep="")
    }

  ##get number of events in each sample
  nr.events.sample <- rep(0, length(GoodSamples))
  names(nr.events.sample) <- GoodSamples

  for(sample in GoodSamples)
    {
      nr.events.sample[sample] <-
        nrow(EventsBySample[EventsBySample[, "Sample"] == sample, ])
    }

  ##create objects for simulated data
  sim.Coverage <- list()
  sim.Scores <- list()
  sim.cma.alter <- list()
  sim.cma.cov <- list()
  sim.cma.samp <- list()

  ##create a vector of correspondence between transcript name and gene-name
  transcripts <- c(rownames(Coverage), fake.transcripts)
  transcripts2genes <- shorten.gene.names(transcripts)
  names(transcripts2genes) <- transcripts

  ##create a vector of correspondence between mutation types and contexts
  mutations2contexts <- c(rep(c("C.in.CpG","G.in.CpG","G.in.GpA","C.in.TpC",
                                "A","C.not.in.CpG.or.TpC","G.not.in.CpG.or.GpA",
                                "T"), each = 3), "ins.del")
  names(mutations2contexts) <- names(passenger.rates)

  ##create a vector of correspondence between mutation types preceded by
  ##"MutationsDiscovery." and mutation types
  mutations2add <- names(passenger.rates)

  names(mutations2add) <- paste("MutationsDiscovery.",
                                names(passenger.rates),
                                sep="")

  for(i in 1:nr.iterations)
    {
      message(paste("Currently simulating data: Iteration #", i, sep=""))
      message(date())

      sim.EventsBySample.i <- c()

      ##if there are no spiked-in genes, keep same Coverage object
      sim.Coverage.i <- Coverage

      ##if there are spiked-in genes, then these objects will be augmented
      ##by some randomly selected coverages and sizes for the spike-ins
      if(length(spiked.set.sizes) > 0)
        {
          ##randomly sample some genes
          sampled.genes <- sample(rownames(Coverage),
                                  length(fake.transcripts))

          sim.Coverage.i[fake.transcripts,] <-
            Coverage[sampled.genes,]
        }

      sim.Coverage[[i]] <- sim.Coverage.i

      ##create object with expected number of mutations, given the coverages and
      ##background rates
      sim.exp.nr.mut.i <-
          as.matrix(sim.Coverage.i[, c(rep(1:8, each = 3), 9)]) *
              matrix(rep(unlist(passenger.rates), nrow(sim.Coverage.i)),
                     ncol = 25, byrow = TRUE)
      rownames(sim.exp.nr.mut.i) <- rownames(sim.Coverage.i)
      colnames(sim.exp.nr.mut.i) <- names(passenger.rates)

      if(pass.null)
      {
          ##get sample-specific constants we're multiplying the passenger null by
          sample.constants <- matrix(rep(nr.events.sample/mean(nr.events.sample),
                                         each=length(passenger.rates)),ncol=length(GoodSamples))
          colnames(sample.constants) <- GoodSamples

          ##make matrix of sample-specific passenger rates
          ##each column represents a sample
          ##(so columns differ just by a multiplicative constant)
          passenger.rates.mat <- matrix(rep(as.numeric(passenger.rates),
                                            length(GoodSamples)),
                                        ncol=length(GoodSamples))
          colnames(passenger.rates.mat) <- GoodSamples
          rownames(passenger.rates.mat) <- names(passenger.rates)
          passenger.rates.mat <- passenger.rates.mat*sample.constants

          ##distribute events over genes for non-spiked-in gene-sets
          for(sample in names(nr.events.sample))
          {
            ##get mutations
            ##message(sample)

            cma.cov.sample <- cma.cov
            cma.cov.sample$Coverage <- cma.cov$Coverage/length(GoodSamples)
            cma.samp.sample <- cma.samp
            cma.samp.sample$NrSamp[cma.samp.sample$NrSamp != 0] <- 1

            Mut <- cma.simulator(cma.cov = cma.cov.sample,
                                 cma.samp = cma.samp.sample,
                                 passenger.rates =
                                 t(data.frame(passenger.rates.mat[,sample])),
                                 eliminate.noval = FALSE)

            genes.with.events <- c()
            contexts <- c()

            ##get the rownames with mutations
            Mut.nonzero <- Mut[rowSums(Mut[,1:25]) > 0,]

            ##get the number of events for each gene
            nr.events.for.gene <- rowSums(Mut.nonzero[,1:25])

            ##write out the genes that have events
            genes.with.events <- rep(rownames(Mut.nonzero),
                                     nr.events.for.gene)

            ##transform Mut.nonzero into a vector to make it easier to manipulate
            Mut.nonzero.as.vect <-
              as.vector(t(as.matrix(Mut.nonzero[,1:25])))
            ##get the types of events for each gene
            types <- rep(rep(colnames(Mut.nonzero[,1:25]),
                             times = nrow(Mut.nonzero)),
                         Mut.nonzero.as.vect)

            ##get the number of events of each type
            events.per.type <-
              Mut.nonzero.as.vect[Mut.nonzero.as.vect > 0]

            contexts <- types

            ##add stuff to sim.EventsBySample.i
            sim.EventsBySample.i <- rbind(sim.EventsBySample.i,
                                          cbind("Mut", sample,
                                                shorten.gene.names(genes.with.events),
                                                mutations2add[contexts]))

            ##see if sample is altered for each of the spiked-in gene-sets
            ##considered
            unifs <- matrix(runif(length(perc.samples)*length(spiked.set.sizes), 0, 1),
                            ncol = length(perc.samples), nrow = length(spiked.set.sizes))

              ##make matrix of prob. that a gene-set is altered in this sample
              prob.set.alt <-
                  matrix(rep(perc.samples, length(spiked.set.sizes)),
                         ncol = length(perc.samples),
                         byrow = TRUE)/100

              ##get the entries where unifs < prob.set.alt -> these will correspond to the altered spiked-in sets
              alt.spiked <- which(prob.set.alt > unifs, arr.ind = TRUE)
              alt.spiked[,"col"] <- perc.samples[alt.spiked[,"col"]]
              alt.spiked[,"row"] <- spiked.set.sizes[alt.spiked[,"row"]]
              ##get the actual names of the spiked-in sets
              alt.spiked.sets <- sapply(as.data.frame(t(alt.spiked)),
                                        function(vect) {paste("gene.set.",
                                                              paste(vect,collapse="."),
                                                              sep="")})

              for(set in alt.spiked.sets)
              {
                  ##if sample is altered, randomly choose a gene
                  ##for it to be altered in - weigh by total expected number
                  ##of mutations

                  ##get genes in this spiked-in set
                  genes.in.spiked.set <-
                      GeneSets[[set]]

                  ##get corresponding transcripts
                  transc.in.spiked.set <-
                      paste(genes.in.spiked.set, "@@00", sep="")

                  spiked.in.exp.nr.mut <-
                      rowSums(sim.exp.nr.mut.i[transc.in.spiked.set, ])

                  altered.gene <-
                      sample(genes.in.spiked.set, 1,
                             prob = spiked.in.exp.nr.mut)

                  ##get altered transcript for altered gene
                  altered.gene <-
                      names(transcripts2genes[transcripts2genes == altered.gene])

                  genes.with.events <-
                      c(genes.with.events, altered.gene)

                  ##add stuff to sim.EventsBySample.i object
                  ##use weights to pick a mutation for this alteration
                  ##the weights are equal to the expected number of mutations
                  ##under the null
                  mutation.gene <- sample(names(passenger.rates),1,replace=FALSE,
                                          prob=sim.exp.nr.mut.i[altered.gene, ])

                  ##add stuff to sim.EventsBySample.i
                  sim.EventsBySample.i <- rbind(sim.EventsBySample.i,
                                                c("Mut", sample,
                                                  shorten.gene.names(altered.gene),
                                                  mutation.gene))
              }
          }
      }
      else
      {
          ##distribute events over genes for non-spiked-in gene-sets
          for(sample in names(nr.events.sample))
          {
              ##sample transcripts with events from all transcripts
              genes.with.events <-
                  sample(rownames(Coverage), nr.events.sample[sample], replace = FALSE)

              ##assign mutations to the genes
              ##mutation.genes <-
              ##    sample(EventsBySample$MutationClass[EventsBySample$Sample == sample])
              ##names(mutation.genes) <- genes.with.events

              ##see if sample is altered for each of the spiked-in gene-sets
              ##considered
              unifs <- matrix(runif(length(perc.samples)*length(spiked.set.sizes), 0, 1),
                              ncol = length(perc.samples), nrow = length(spiked.set.sizes))

              ##make matrix of prob. that a gene-set is altered in this sample
              prob.set.alt <-
                  matrix(rep(perc.samples, length(spiked.set.sizes)),
                         ncol = length(perc.samples),
                         byrow = TRUE)/100

              ##get the entries where unifs < prob.set.alt -> these will correspond to the altered spiked-in sets
              alt.spiked <- which(prob.set.alt > unifs, arr.ind = TRUE)
              alt.spiked[,"col"] <- perc.samples[alt.spiked[,"col"]]
              alt.spiked[,"row"] <- spiked.set.sizes[alt.spiked[,"row"]]
              ##get the actual names of the spiked-in sets
              alt.spiked.sets <- sapply(as.data.frame(t(alt.spiked)),
                                        function(vect) {paste("gene.set.",
                                                              paste(vect,collapse="."),
                                                              sep="")})
              ##randomly chose a gene to be altered from each spiked-in set
              altered.spiked.genes <- sapply(GeneSets[alt.spiked.sets], "sample", 1)
              ##now get transcripts for these genes
              genes.with.events.spiked <- paste(altered.spiked.genes, "@@00", sep="")

              ##mutation.genes.spiked <-
              ##    sample(names(passenger.rates),length(genes.with.events.spiked),
              ##           replace=FALSE,
              ##           prob=
              ##           table(EventsBySample$MutationClass)[names(passenger.rates)])
              ##names(mutation.genes.spiked) <- genes.with.events.spiked

              genes.with.events <- c(genes.with.events, genes.with.events.spiked)
              ##mutation.genes <- c(mutation.genes, mutation.genes.spiked)

              mutation.genes <- rep("", length(genes.with.events))
              names(mutation.genes) <- genes.with.events

              for(gene in genes.with.events)
              {
                  mutation.gene <- sample(names(passenger.rates),1,replace=FALSE,
                                          prob=sim.exp.nr.mut.i[gene, ])
                  mutation.genes[gene] <- mutation.gene
              }

              ##add stuff to sim.EventsBySample.i
              sim.EventsBySample.i <- rbind(sim.EventsBySample.i,
                                            cbind("Mut", sample,
                                                  shorten.gene.names(genes.with.events),
                                                  mutation.genes))

              mutation.genes <- rep("", length(genes.with.events))
              names(mutation.genes) <- genes.with.events
          }
      }

      colnames(sim.EventsBySample.i) <- colnames(EventsBySample)
      rownames(sim.EventsBySample.i) <- NULL

      sim.EventsBySample.i <- as.data.frame(sim.EventsBySample.i)

      ##make cma.alter, cma.cov, and cma.samp objects
      WT <- rep("", nrow(sim.EventsBySample.i))
      Mut <- rep("", nrow(sim.EventsBySample.i))
      Context <- rep("", nrow(sim.EventsBySample.i))

      for(j in 1:nrow(sim.EventsBySample.i))
        {
          if(length(grep("A.to.", sim.EventsBySample.i[j,4])) > 0)
            {
              WT[j] <- "A"
            }
          if(length(grep("G.in.", sim.EventsBySample.i[j,4])) > 0 |
             length(grep("G.not.in.", sim.EventsBySample.i[j,4])) > 0)
            {
              WT[j] <- "G"
            }
          if(length(grep("C.in.", sim.EventsBySample.i[j,4])) > 0 |
             length(grep("C.not.in.", sim.EventsBySample.i[j,4])) > 0)
            {
              WT[j] <- "C"
            }
          if(length(grep("T.to.", sim.EventsBySample.i[j,4])) > 0)
            {
              WT[j] <- "T"
            }
          if(length(grep("to.G", sim.EventsBySample.i[j,4])) > 0)
            {
              Mut[j] <- "G"
            }
          if(length(grep("to.C", sim.EventsBySample.i[j,4])) > 0)
            {
              Mut[j] <- "C"
            }
          if(length(grep("to.A", sim.EventsBySample.i[j,4])) > 0)
            {
              Mut[j] <- "A"
            }
          if(length(grep("to.T", sim.EventsBySample.i[j,4])) > 0)
            {
              Mut[j] <- "T"
            }
          if(length(grep("in.GpA.to", sim.EventsBySample.i[j,4])) > 0)
            {
              Context[j] <- "GpA"
            }
          if(length(grep("in.CpG.to", sim.EventsBySample.i[j,4])) > 0)
            {
              Context[j] <- "CpG"
            }
          if(length(grep("in.TpC.to", sim.EventsBySample.i[j,4])) > 0)
            {
              Context[j] <- "TpC"
            }
          if(length(grep("ins.del", sim.EventsBySample.i[j,4])) > 0)
            {
              Mut[j] <- "ins.del"
              Context[j] <- "All"
            }
        }

      sim.cma.alter.i <- data.frame(Gene = sim.EventsBySample.i$Symbol,
                                    Type = sim.EventsBySample.i$Event,
                                    Sample = sim.EventsBySample.i$Sample,
                                    Screen = "Disc",
                                    WTNuc = WT,
                                    Context = Context,
                                    MutNuc = Mut)

      GeneSymb <- shorten.gene.names(rownames(sim.Coverage[[i]]))

      cc <- sim.Coverage[[i]][!duplicated(GeneSymb),]
      rownames(cc) <- unique(GeneSymb)

      ##make coverage object
      sim.cma.cov.i <- data.frame(Gene = rep(rownames(cc), each=18),
                                  Screen = c("Disc","Prev"),
                                  WTNuc = rep(c("C","G","G","C",
                                    "A","C","G","T",""),
                                    each = 2),
                                  Context = rep(c("CpG","CpG","GpA","TpC",
                                    "","","","","All"),
                                    each = 2),
                                  Coverage =
                                  as.numeric(unlist(t(cc[,c(1,10,2,11,3,12,4,13,5,14,6,15,
                                                            7,16,8,17,9,18)]))))

      ##make object with number of samples in which each gene is sequenced
      sim.cma.samp.i <- data.frame(Gene = rep(rownames(cc), each=2),
                                   Screen = c("Disc","Prev"),
                                   NrSamp = as.numeric(unlist(t(cc[,19:20]))))

      sim.Scores[[i]] <- cma.scores(cma.alter = sim.cma.alter.i,
                                    cma.cov = sim.cma.cov.i,
                                    cma.samp = sim.cma.samp.i,
                                    passenger.rates = passenger.rates)

      sim.cma.alter[[i]] <- sim.cma.alter.i
      sim.cma.cov[[i]] <- sim.cma.cov.i
      sim.cma.samp[[i]] <- sim.cma.samp.i
    }

    return(list(GeneSets = GeneSets,
                cma.alter = sim.cma.alter,
                cma.cov = sim.cma.cov,
                cma.samp = sim.cma.samp,
                Scores = sim.Scores))
  }

