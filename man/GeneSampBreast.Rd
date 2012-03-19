\name{GeneSampBreast}
\docType{data}
\alias{GeneSampBreast}

\title{Data from the Wood et al. 2007 and Sjoeblom et al. 2006 studies: Number of samples
  for each gene and screen type}
\description{Number of samples analyzed for each gene and 
  screen type from the breast cancer portion of the Wood et al. 2007
  and Sjoeblom et al. 2006 studies.}

\usage{data(WoodBreast07)}
\format{
  The number of samples in the breast cancer portions of the
  Wood et al. and Sjoeblom et al. studies, broken down by
  \emph{gene} and \emph{screen} (Discovery and Prevalence).
  The object is a data frame, with the variables:
  Gene, Screen, and NrSamp (number of samples).
  The two possible values for Screen are Disc ("Discovery")
  and Prev ("Prevalence").
}

\references{
  Wood LD, Parsons DW, Jones S, Lin J, Sjoblom T, Leary RJ, Shen D,
  Boca SM, Barber T, Ptak J, et al.
  The genomic landscapes of human breast and colorectal cancers.
  \emph{Science.} DOI:10.1126/science.1145720

  Sjoeblom T, Jones S, Wood LD, Parsons DW, Lin J, Barber T, Mandelker D,
  Leary R, Ptak J, Silliman N, et al.
  The consensus coding sequences of human breast and colorectal cancers.
  \emph{Science.} DOI: 10.1126/science.1133427
  
  Parmigiani G, Lin J, Boca S, Sjoeblom T, Kinzler KW,
  Velculescu VE, Vogelstein B. Statistical methods for the analysis of
  cancer genome sequencing data. 
  \url{http://www.bepress.com/jhubiostat/paper126/}
}
\seealso{
  \code{cma.scores}, \code{cma.fdr},
  \code{cma.set.stat}, \code{cma.set.sim},
  \code{SimMethodsSims-class},
  \code{GeneAlterBreast}, 
  \code{GeneCovBreast}
}

\keyword{datasets}
