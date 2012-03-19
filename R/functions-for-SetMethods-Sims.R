extract.sims.method <- function(object, method){
               results <- object@results

               tmp <- sapply(results, function(res,m) {res[,m]}, method)
               rownames(tmp) <- rownames(object@results[[1]])

               tmp
           }

combine.sims <- function(obj1, obj2)
{
    new.obj <- new("SetMethodsSims")
    if(!identical(obj1@null.dist, obj2@null.dist) |
       !identical(obj1@perc.samples, obj2@perc.samples) |
       !identical(obj1@spiked.set.sizes, obj2@spiked.set.sizes))
    {
        stop("Cannot combine these two simulations - Simulation parameters are not the same")
    }
    else
    {
       new.obj@null.dist <- obj1@null.dist
       new.obj@perc.samples <- obj1@perc.samples
       new.obj@spiked.set.sizes <- obj1@spiked.set.sizes
       new.obj@GeneSets <- obj1@GeneSets
       new.obj@cma.alter <- c(obj1@cma.alter, obj2@cma.alter)
       new.obj@cma.cov <- c(obj1@cma.cov, obj2@cma.cov)
       new.obj@cma.samp <- c(obj1@cma.samp, obj2@cma.samp)
       new.obj@results <- c(obj1@results, obj2@results)
   }
    new.obj
}
