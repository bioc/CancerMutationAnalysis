##function to get p-values from q.probs for a specific set for subject-based method
get.p.vals.from.q.probs <- function(q, obs)
{
    prob.altered.samples <- poisbin(q)

    ##sum all the probabilities to get as many or more altered samples as observed (note that the probabilities in prob.altered.samples start at 0)
    p.vals.direct <- sum(prob.altered.samples[(obs+1):length(prob.altered.samples)])

    p.vals.direct
}
