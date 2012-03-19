\name{GeneAlterBreast}
\docType{data}
\alias{GeneAlterBreast}

\title{Data from the Wood et al. 2007 and Sjoeblom et al. 2006
  studies: Alterations for every gene and sample}
\description{Somatic alterations for each gene and tumor sample from the
  breast cancer portion of the Wood et al. 2007 and Sjoeblom et al. 2006
  studies.}

\usage{data(WoodBreast07)}
\format{
  The somatic mutations in the breast cancer portions of the
  Wood et al. and Sjoeblom et al. studies,
  broken down by \emph{gene}, \emph{type} (point mutation, amplification,
  or deletion), \emph{sample},
  \emph{screen} (Discovery or Prevalence),
  and, for point mutations, \emph{mutation type},
  composed of the wild type nucleotide, its context,
  and the mutated nucleotide. The object is a data frame,
  with the variables: Gene, Type, Sample, Screen, WTNuc
  (wild type nucleotide), Context, and MutNuc (mutated
  nucleotide). The two possible values
  for Screen are Disc ("Discovery") and Prev ("Prevalence").
  The three possible values for Type are Mut (point mutations),
  Amp (large amplifications), and Del (large deletions.)
  Indels have a "" entry for WTNuc, an "All"
  entry for Context, and a "ins.del" entry for MutNuc.
  Large amplifications and deletions have "" entries for
  WTNuc, Context, and MutNuc. For this study,
  only point mutation are available.
  
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
  \code{GeneCovBreast}, 
  \code{GeneSampBreast}
}

\keyword{datasets}
