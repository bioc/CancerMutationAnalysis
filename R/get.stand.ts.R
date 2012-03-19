get.stand.ts <-
function(EventsBySampleBySet = EventsBySampleBySet,
                         GoodSamples = GoodSamples,
                         nr.events.sample = nr.events.sample,
                         GeneSets = GeneSets,
                         GeneSymb = GeneSymb,
                         return.q = FALSE) 
  {
    ##function to get number of altered samples for each set
    get.nr.altered.samples <- function(EventsBySampleForSet)
      {
        sum(EventsBySampleForSet > 0)
      }
    ##get nr. of altered samples for each set
    AlteredSamplesPerSet = apply(EventsBySampleBySet, 1, get.nr.altered.samples)
    names(AlteredSamplesPerSet) <- names(GeneSets)
    
    ##get the null probability that a gene is not altered in each of the samples
    pr.gene.not.altered <- (length(GeneSymb) - nr.events.sample)/length(GeneSymb)
    names(pr.gene.not.altered) <- GoodSamples
    ##get lengths (sizes) of the sets
    lengths.sets <- sapply(GeneSets, length)
    ##create matrix of probabilities q_si (see pdf document)
    q <- matrix(0, nrow = length(GeneSets), ncol = length(GoodSamples))
    rownames(q) <- names(GeneSets)
    colnames(q) <- GoodSamples

    for(sample in GoodSamples)
      {
        q[, sample] <- 1 - pr.gene.not.altered[sample]^lengths.sets
      }
    
    ##calculate standardized T_s
    stand.ts <- rep(0, length = length(GeneSets))
    names(stand.ts) <- names(GeneSets)

    ##function to get standardized T_s from q.probs for a specific set
    get.stand.ts.from.q.probs <- function(set, q)
      {
          ##print(set)
          ##set all the probabilities corresponding to all the samples for gene-set
          ##"set"
          q.probs <- q[set, ]
          ##print(q.probs)
          
          ##get observed t_s
          ts <- AlteredSamplesPerSet[set]

          ##standardize t_s
          exp.ts <- length(GoodSamples)-sum(1-q.probs)
          var.ts <- sum(1-q.probs)-sum((1-q.probs)^2)
          stand.ts <- (ts - exp.ts)/sqrt(var.ts)
      }

    stand.ts <- lapply(names(GeneSets), get.stand.ts.from.q.probs, q)

    stand.ts <- unlist(stand.ts)
    names(stand.ts) <- names(GeneSets)
    return(stand.ts)
    if(return.q == TRUE)
      {
        return(list(stand.ts = stand.ts,
                    q.matrix = q))
      }
  }

