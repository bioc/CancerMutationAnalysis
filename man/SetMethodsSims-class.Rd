\name{SetMethodsSims-class}
\docType{class}
\alias{SetMethodsSims-class}
\alias{SetMethodsSims-method}
\alias{show,SetMethodsSims-method}

\title{Class representation for depositing output from simulations.}
\description{Stores results from the \code{sim.data.p.values} function.}

\section{Objects from the class}{

 New objects can be created using calls of the form 
 \code{new("SetMethodsSims", null.dist, perc.samples, spiked.set.sizes,
       GeneSizes, GeneSets, Coverage, EventsBySample, Mutations, Scores,
       results)}

}

\section{Slots}{
  \describe{
    \item{\code{null.dist}:}{Object of class \code{"character"}.
      Can be either "Passenger null" or
      "Permutation null," depending on what 
      method is used to get the null data.} 
    \item{\code{perc.samples}:}{Object of class \code{"numeric"}.
      Vector representing the probabilities 
      of the spiked-in gene-sets being
      altered in any given sample, as percentages; 
      for example
      \code{perc.samples = c(75, 90)} 
      means that these probabilities are
      \code{0.75} and \code{0.90}.} 
    \item{\code{spiked.set.sizes}:}{Object of class \code{"numeric"}. 
      Vector representing the sizes, in genes, 
      of the spiked-in gene-sets;
      for example, if 
      \code{perc.samples = c(75, 90)} and
      \code{spiked.set.sizes = c(50, 100)}, 
      there would be 4 spiked-in sets,
      one with 50 genes and probability of 
      being altered of 0.75 in each sample,
      one with 50 genes and probability of 
      being altered of 0.90 in each sample,
      one with 100 genes and probability of 
      being altered of 0.75 in each sample, and
      one with 100 genes and probability of 
      being altered of 0.90 in each sample.} 
    \item{\code{GeneSets}:}{Object of class \code{"list"}.
      The entries of the list correspond to
      gene-sets and give the genes annotated to
      them.} 
    \item{\code{cma.alter}:}{Object of class \code{"list"}.
      The entries of the list are objects similar to
      the \code{GeneAlter} objects and correspond to the
      simulation iterations.}
    \item{\code{cma.cov}:}{Object of class \code{"list"}.
      The entries of the list are objects similar to
      the \code{GeneCov} objects and correspond to the
      simulation iterations.}
    \item{\code{cma.samp}:}{Object of class \code{"list"}.
      The entries of the list are objects similar to
      the \code{GeneSamp} objects and correspond to the
      simulation iterations.}
    \item{\code{Scores}:}{Object of class \code{"list"}.
      The entries of this list are the output of
      \code{cma.scores} and correspond to the
      simulation iterations.}
    \item{\code{results}:}{Object of class \code{"list"}.
      The entries of this list are the output of
      \code{cma.set.stat} and correspond to the
      simulation iterations.}
  }
}

\section{Methods}{
\describe{
\item{show}{\code{signature(object = "SetMethodsSims")}}
}
}

\references{
  Boca SM, Kinzler KW, Velculescu VE, Vogelstein B,
  Parmigiani G.
  Patient-oriented gene-set analysis for cancer mutation data.
  \emph{Genome Biology.} DOI: 10.1186/gb-2010-11-11-r112
}
\author{
Simina M. Boca, Giovanni Parmigiani.
}

\seealso{
\code{GeneCov}, \code{GeneSamp}, \code{GeneAlter},
\code{cma.set.sim}, \code{cma.set.stat},
\code{combine.sims}, \code{extract.sims.method}
}

\keyword{datasets}
