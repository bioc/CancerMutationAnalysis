pass.null.het.p.values <- function(GeneSets,
                                   GoodSamples,
                                   GeneSymb,
                                   nr.events.sample,
                                   passenger.rates,
                                   Coverage,
                                   EventsBySample,
                                   EventsBySampleBySet)
  {
    ##get sample-specific constants we're multiplying the passenger null by
    sample.constants <- matrix(rep(nr.events.sample/mean(nr.events.sample),
                                   each=length(passenger.rates)),ncol=length(GoodSamples))
    colnames(sample.constants) <- GoodSamples

    ##make matrix of sample-specific passenger rates
    passenger.rates.mat <- matrix(rep(as.numeric(passenger.rates),length(GoodSamples)),ncol=length(GoodSamples))
    colnames(passenger.rates.mat) <- GoodSamples
    rownames(passenger.rates.mat) <- names(passenger.rates)
    passenger.rates.mat <- passenger.rates.mat*sample.constants

    ##get probability of A not getting mutated, T not getting mutated etc., under the passenger null
    prob.nucl.not.alt.mat <- matrix(0,ncol=length(GoodSamples),nrow=9)
    colnames(prob.nucl.not.alt.mat) <- GoodSamples
    rownames(prob.nucl.not.alt.mat) <- c("C.in.CpG", "G.in.CpG", "G.in.GpA" , "C.in.TpC",
                                         "A", "C.not.in.CpG.or.TpC", "G.not.in.CpG.or.GpA",
                                         "T", "ins.del")
    
    prob.nucl.not.alt.mat["C.in.CpG",] <-
      1-passenger.rates.mat["C.in.CpG.to.G",]-
        passenger.rates.mat["C.in.CpG.to.A",]-
          passenger.rates.mat["C.in.CpG.to.T",]
    prob.nucl.not.alt.mat["G.in.CpG",] <-
      1-passenger.rates.mat["G.in.CpG.to.C",]-
        passenger.rates.mat["G.in.CpG.to.A",]-passenger.rates.mat["G.in.CpG.to.T",]
    prob.nucl.not.alt.mat["G.in.GpA",] <-
      1-passenger.rates.mat["G.in.GpA.to.C",]-
        passenger.rates.mat["G.in.GpA.to.A",]-passenger.rates.mat["G.in.GpA.to.T",]
    prob.nucl.not.alt.mat["C.in.TpC",] <-
      1-passenger.rates.mat["C.in.TpC.to.G",]-
        passenger.rates.mat["C.in.TpC.to.A",]-passenger.rates.mat["C.in.TpC.to.T",]
    prob.nucl.not.alt.mat["C.not.in.CpG.or.TpC",] <-
      1-passenger.rates.mat["C.not.in.CpG.or.TpC.to.G",]-
        passenger.rates.mat["C.not.in.CpG.or.TpC.to.A",]-
          passenger.rates.mat["C.not.in.CpG.or.TpC.to.T",]
    prob.nucl.not.alt.mat["G.not.in.CpG.or.GpA",] <-
      1-passenger.rates.mat["G.not.in.CpG.or.GpA.to.C",]-
        passenger.rates.mat["G.not.in.CpG.or.GpA.to.A",]-
          passenger.rates.mat["G.not.in.CpG.or.GpA.to.T",]
    prob.nucl.not.alt.mat["A",] <-
      1-passenger.rates.mat["A.to.C",]-passenger.rates.mat["A.to.G",]-passenger.rates.mat["A.to.T",]
    prob.nucl.not.alt.mat["T",] <-
      1-passenger.rates.mat["T.to.C",]-passenger.rates.mat["T.to.G",]-passenger.rates.mat["T.to.A",]
    prob.nucl.not.alt.mat["ins.del",] <-
      1-passenger.rates.mat["ins.del",]

    ##declare object with probability of gene being altered under passenger null
    prob.pass.null.mat <- matrix(1, nrow = nrow(Coverage),
                                 ncol = length(GoodSamples))
    rownames(prob.pass.null.mat) <- names(GeneSymb)
    colnames(prob.pass.null.mat) <- GoodSamples

    ##calculate probability of gene NOT having any mutations under the null first
    Coverage.sample <- Coverage[,1:9]/length(GoodSamples)
    colnames(Coverage.sample) <- rownames(prob.nucl.not.alt.mat) 
    for(context in rownames(prob.nucl.not.alt.mat))
      {
        for(sample in GoodSamples)
          {
            prob.pass.null.mat[,sample] <-
              prob.pass.null.mat[,sample]*
                prob.nucl.not.alt.mat[context,sample]^Coverage.sample[,context]
          }
      }
    
    ##now subtract from 1 to get probability of gene having AT LEAST 1 somatic mutation under the null
    prob.pass.null.mat <-
          1-prob.pass.null.mat

    ##get lengths (sizes) of the sets
    lengths.sets <- sapply(GeneSets, length)

    ##create matrix of probabilities q_si (see pdf document)
    q <- matrix(0, nrow = length(GeneSets), ncol = length(GoodSamples))
    rownames(q) <- names(GeneSets)
    colnames(q) <- GoodSamples

    ##make data-frame out of gene-sets
    sets.frame <- data.frame(sets = rep(names(GeneSets),
                             lengths.sets),
                             genes = unlist(GeneSets))
    rownames(sets.frame) <- NULL

    nr.GoodSamples <- length(GoodSamples)
    for(sample in 1:nr.GoodSamples)
      {
          ##print(sample)
          prob.pass.null.sample <- prob.pass.null.mat[,sample]
          ##make data-frame out of prob.pass.null.sample
          prob.pass.null.sample.frame <- data.frame(genes = names(prob.pass.null.sample),
                                                    probs = prob.pass.null.sample)

          ##make data-frame at set level
          probs.for.sets <- prob.pass.null.sample[as.character(sets.frame$genes)]

          probs.for.sets.list <- list()
          curr.length <- 0
          for(set in 1:length(GeneSets))
          {
              length.s <- lengths.sets[set]
              probs.for.sets.list[[set]] <-
                  probs.for.sets[(curr.length+1):
                                 (curr.length+length.s)]
              curr.length <- curr.length+length.s
          }
          names(probs.for.sets.list) <- names(GeneSets)

          q[,sample] <- 1-sapply(probs.for.sets.list, function(p) {prod(1-p)})

      }

    ##get statistics for each set
    EE <- EventsBySampleBySet
    EE[EE > 0] <- 1
    Ts.new <- rowSums(EE*1/q)

    ##calculate p-values directly
    q.temp <- as.data.frame(t(q))
    one.over.q <- as.data.frame(t(1/q))
    q <- q.temp

    p.vals.direct <- mapply(get.p.vals.from.q.probs.het,
                            q, one.over.q, Ts.new)

    p.vals.direct
  }
