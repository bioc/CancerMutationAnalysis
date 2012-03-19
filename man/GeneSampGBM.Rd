\name{GeneSampGBM}
\docType{data}
\alias{GeneSampGBM}

\title{Data from the Parsons et al. 2008 study: Number of samples
  for each gene and screen type}
\description{Number of samples analyzed for each gene and 
  screen type from the
  Parsons et al. 2008 glioblastoma multiforme (GBM) study.}

\usage{data(ParsonsGBM08)}
\format{
  The number of samples in the GBM study from Parsons et al.,
  broken down by
  \emph{gene} and \emph{screen} (Discovery and Prevalence).
  The object is a data frame, with the variables:
  Gene, Screen, and NrSamp (number of samples).
  The two possible values for Screen are Disc ("Discovery")
  and Prev ("Prevalence").
}

\references{
  Parsons DW, Jones S, Zhang X, Lin JCH, Leary RJ, Angenendt P, Mankoo P,
  Carter H, Siu I, et al. 
  An integrated genomic analysis of human glioblastoma multiforme. 
  \emph{Science.} DOI: 10.1126/science.1164382

  Parmigiani G, Lin J, Boca S, Sjoeblom T, Kinzler KW,
  Velculescu VE, Vogelstein B. Statistical methods for the analysis of
  cancer genome sequencing data. 
  \url{http://www.bepress.com/jhubiostat/paper126/}
}
\seealso{
  \code{cma.scores}, \code{cma.fdr},
  \code{cma.set.stat}, \code{cma.set.sim},
  \code{SimMethodsSims-class},
  \code{GeneAlterGBM}, 
  \code{GeneCovGBM}
}

\keyword{datasets}
