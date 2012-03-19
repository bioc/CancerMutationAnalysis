\name{GeneCovBreast}
\docType{data}
\alias{GeneCovBreast}

\title{Data from the Wood et al. 2007 and Sjoeblom et al. 2006 studies:
  Total number of nucleotides "at risk" ("coverage")}
\description{Total numbers of nucleotides "at risk" that were 
  successfully sequenced in RefSeq genes in the breast cancer portion of
  the Wood et al. 2007 and Sjoeblom et al. 2006
  studies.}

\usage{data(WoodBreast07)}
\format{
  Total number of nucleotides available for mutations ("coverage") 
  in the breast cancer portion of the Wood et al. and Sjoeblom et al. studies, 
  broken down by \emph{gene}, \emph{screen} (Discovery or Prevalence),
  and \emph{mutation type}, composed of the wild type nucleotide
  and its context. The object is a data frame, with the
  variables: Gene, Screen, WtNuc (wild type nucleotide),
  Context, and Coverage. The two possible values
  for Screen are Disc ("Discovery") and Prev ("Prevalence").
  The nucleotides availables for indels are all the successfully
  sequenced nucleotides in a gene; the corresponding rows
  have a "" entry for WtNuc and an "All" entry for "Context."
  The nucleotides availables for other
  mutations are excluding nucleotides which can only give rise to
  synonymous mutations. 
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
  \code{GeneSampBreast}
}

\keyword{datasets}
