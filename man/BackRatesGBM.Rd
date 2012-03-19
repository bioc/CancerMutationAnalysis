\name{BackRatesGBM}
\docType{data}
\alias{BackRatesGBM}

\title{Data from the Parsons et al. 2008 study: Background mutation rates}
\description{Background rates for somatic mutations used in the
  Parsons et al. 2008 glioblastoma multiforme (GBM) study.}

\usage{data(ParsonsGBM08)}
\format{
  The background rates for somatic mutations used in the
  Parsons et al. GBM study, broken down by mutation type.
  The object is a data frame,
  with the variables representing the 25 different mutation types,
  and the rows specifying whether the estimates of the
  backgrond rates are "Upper," "Median," or "Lower."
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
  \code{GeneCovGBM}, 
  \code{GeneSampGBM},
  \code{GeneAlterGBM}
}

\keyword{datasets}
