poisbin <-
function(theta)
{
    ## theta is a length n vector of success probabilities
    ## for the Bernoulli trials.  
    n <- length(theta)
    cc <- rep(NA,n+1)
    cc[1] <- 1-theta[1]
    cc[2] <- theta[1]
    bb <- cc  # previous value
    for( i in 2:n )
    {
        cc[1] <- bb[1]*(1-theta[i]) 
        cc[i+1] <- bb[i] * theta[i] 
        for( j in 2:i )
        { 
            cc[j] <- bb[j-1]*theta[i]+bb[j]*(1-theta[i])
        }
        bb <- cc
        ##	 print(i)
    }
    return(cc)  # Poisson binomial prob mass function on 0:n
}

