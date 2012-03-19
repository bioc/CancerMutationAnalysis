##function to get p-values for "heterogeneous" subject-based method (i.e. the statistic is a linear combination of the Bernoulli random variables, with the weights being the inverse probabilities of success)
perm.null.het.p.values <- function(EventsBySampleBySet = EventsBySampleBySet,
                                   GoodSamples = GoodSamples,
                                   GeneSymb = GeneSymb,
                                   nr.events.sample = nr.events.sample,
                                   GeneSets = GeneSets,
                                   return.q = FALSE)
{
    ##get the null probability that a gene is not altered in each of the samples
    nr.genes <- length(GeneSymb)
    pr.gene.not.altered <- 1 - nr.events.sample/nr.genes
    names(pr.gene.not.altered) <- GoodSamples
    ##get lengths (sizes) of the sets
    lengths.sets <- sapply(GeneSets, length)

    ##create matrix of probabilities q_si (see pdf document)
    ##pr.gene.not.altered <- as.data.frame(t(pr.gene.not.altered))
    ##q <- 1-sapply(pr.gene.not.altered,
    ##              function(prob, len) {prob^len},
    ##              lengths.sets)
    
    ni <- nr.events.sample
    ms <- lengths.sets
    Ni <- t(t(rep(1,length(ms)))) %*% t(ni)
    Ms <- t(t(ms)) %*% t(rep(1,length(ni)))
    G <- nr.genes
    
    q <- 1-dhyper(0, Ni, G-Ni, Ms)

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

    names(p.vals.direct) <- names(GeneSets)
    return(p.vals.direct)
    if(return.q == TRUE)
      {
        return(list(p.vals = p.vals.direct,
                    q.matrix = q))
      }
  }

