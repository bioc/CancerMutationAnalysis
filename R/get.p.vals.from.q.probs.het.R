##function to get p-values from q.probs for a specific set for "heterogeneous" subject-based method (i.e. the statistic is a linear combination of the Bernoulli random variables)
##lin.comb.mat = matrix of linear combinations
##obs.vect = vector of observed test statistics
get.p.vals.from.q.probs.het <- function(q, lin.comb, obs)
  {
    ##print("Oh yeah! Minus 3! Go me!")
    
    ##sum all the probabilities to get as many or more altered samples as observed (note that the probabilities in prob.altered.samples start at 0)
    p.vals.direct <-
      1 - get.cdf.het.berns(n = length(q),
                            theta.n = q,
                            alpha.n = lin.comb,
                            j = obs - 10^(-3))
    
    p.vals.direct
  }
