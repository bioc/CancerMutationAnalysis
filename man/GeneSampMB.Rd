\name{GeneSampMB}
\docType{data}
\alias{GeneSampMB}

\title{Data from the Parsons et al. 2011 study: Number of samples
  for each gene and screen type}
\description{Number of samples analyzed for each gene and 
  screen type from the
  Parsons et al. 2011 medulloblastoma (MB) study.}

\usage{data(ParsonsMB11)}
\format{
  The number of samples in the MB study from Parsons et al.,
  broken down by
  \emph{gene} and \emph{screen} (Discovery and Prevalence).
  The object is a data frame, with the variables:
  Gene, Screen, and NrSamp (number of samples).
  The two possible values for Screen are Disc ("Discovery")
  and Prev ("Prevalence").
}

\references{
  Parsons DW, Li M, Zhang X, Jones S, Leary RJ, Lin J, Boca SM, Carter
  H, Samayoa J, Bettegowda C, et al.
  The genetic landscape of the childhood cancer medulloblastoma.
  \emph{Science.} DOI: 10.1126/science.1198056
  
  Parmigiani G, Lin J, Boca S, Sjoeblom T, Kinzler KW,
  Velculescu VE, Vogelstein B. Statistical methods for the analysis of
  cancer genome sequencing data. 
  \url{http://www.bepress.com/jhubiostat/paper126/}
}
\seealso{
  \code{cma.scores}, \code{cma.fdr},
  \code{cma.set.stat}, \code{cma.set.sim},
  \code{SimMethodsSims-class},
  \code{GeneAlterMB}, 
  \code{GeneCovMB}
}

\keyword{datasets}
