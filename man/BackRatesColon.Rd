\name{BackRatesColon}
\docType{data}
\alias{BackRatesColon}

\title{Data from the Wood et al. 2007 study: Background mutation rates}
\description{Background rates for somatic mutations used in the
  colon cancer portion of the Wood et al. 2007 study.}

\usage{data(WoodColon07)}
\format{
  The background rates for somatic mutations used in the
  colon cancer portion of the Wood et al. study,
  broken down by mutation type.
  The object is a data frame,
  with the variables representing the 25 different mutation types,
  and the rows specifying whether the estimates of the
  backgrond rates are "Lower," "Median," or "Upper," as well
  as whether or not the rates are separately estimated for
  the prevalence screen (denoted by "SepPrev").
}

\references{
  Wood LD, Parsons DW, Jones S, Lin J, Sjoblom T, Leary RJ, Shen D,
  Boca SM, Barber T, Ptak J, et al.
  The genomic landscapes of human breast and colorectal cancers.
  \emph{Science.} DOI:10.1126/science.1145720

  Parmigiani G, Lin J, Boca S, Sjoeblom T, Kinzler KW,
  Velculescu VE, Vogelstein B. Statistical methods for the analysis of
  cancer genome sequencing data. 
  \url{http://www.bepress.com/jhubiostat/paper126/}
}
\seealso{
  \code{cma.scores}, \code{cma.fdr},
  \code{cma.set.stat}, \code{cma.set.sim},
  \code{SimMethodsSims-class},
  \code{GeneCovColon}, 
  \code{GeneSampColon},
  \code{GeneAlterBreast}
}

\keyword{datasets}
