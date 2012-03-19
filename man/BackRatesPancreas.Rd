\name{BackRatesPancreas}
\docType{data}
\alias{BackRatesPancreas}

\title{Data from the Jones et al. 2008 study: Background mutation rates}
\description{Background rates for somatic mutations used in the
  Jones et al. 2008 pancreatic cancer study.}

\usage{data(JonesPancreas08)}
\format{
  The background rates for somatic mutations used in the
  Jones et al. pancreatic cancer study,
  broken down by mutation type.
  The object is a data frame,
  with the variables representing the 25 different mutation types,
  and the rows specifying whether the estimates of the
  backgrond rates are "Upper," "Median," or "Lower."
}

\references{
  Jones S, Zhang X, Parsons DW, Lin JC, Leary RJ, Angenendt P, Mankoo P,
  Carter H, Kamiyama H, Jimeno A, et al.
  Core signaling pathways in human pancreatic cancers revealed by global
  genomic analyses.
  \emph{Science.} DOI: 10.1126/science.1164368

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
