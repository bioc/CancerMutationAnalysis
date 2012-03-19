setMethod("show", "SetMethodsSims",
          function(object){
              nr.sims <- length(object@results)
              nr.sets <- length(object@GeneSets)

              cat("Simulation results for gene-set analysis of mutations\n")

              cat("  Data-generating mechanism :", object@null.dist, "\n")
              cat("  Number of simulations     :", nr.sims, "\n")
              cat("  Number of gene-sets       :", nr.sets, "\n")
              cat("     Original  :",
                  nr.sets - length(object@perc.samples)*length(object@spiked.set.sizes), "\n")
              cat("     Spiked-in :",
                  length(object@perc.samples)*length(object@spiked.set.sizes), "\n")
              cat("  Spiked-in sets            :\n")
              cat("     Probability of being mutated in a set :",
                  object@perc.samples/100, "\n")
              cat("     Length (as number of genes)           :",
                  object@spiked.set.sizes, "\n")
              cat("\n")
          })


