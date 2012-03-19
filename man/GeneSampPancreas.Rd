\name{GeneSampPancreas}
\docType{data}
\alias{GeneSampPancreas}

\title{Data from the Jones et al. 2008 study: Number of samples
  for each gene and screen type}
\description{Number of samples analyzed for each gene and 
  screen type from the
  Jones et al. 2008 pancreatic cancer study.}

\usage{data(JonesPancreas08)}
\format{
  The number of samples in the pancreatic cancer study from Jones et
  al., broken down by
  \emph{gene} and \emph{screen} (Discovery and Prevalence).
  The object is a data frame, with the variables:
  Gene, Screen, and NrSamp (number of samples).
  The two possible values for Screen are Disc ("Discovery")
  and Prev ("Prevalence").
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
  \code{GeneAlterPancreas}, 
  \code{GeneCovPancreas}
}

\keyword{datasets}
