perm.null.p.values <- function(EventsBySampleBySet = EventsBySampleBySet,
                               GoodSamples = GoodSamples,
                               nr.events.sample = nr.events.sample,
                               GeneSets = GeneSets,
                               GeneSymb,
                               return.q = FALSE)
{
    ##get nr. of altered samples for each set
    EE <- EventsBySampleBySet
    EE[EE > 0] <- 1
    AlteredSamplesPerSet <- rowSums(EE)

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

    ##calculate p-values directly
    q <- as.data.frame(t(q))
    AlteredSamplesPerSet <- as.data.frame(t(AlteredSamplesPerSet))

    p.vals.direct <- mapply(get.p.vals.from.q.probs,
                            q, AlteredSamplesPerSet)

    if(return.q == TRUE)
    {
        return(list(p.vals = p.vals.direct,
                    q.matrix = q))
    }
    else
    {
        return(p.vals.direct)
    }
  }

