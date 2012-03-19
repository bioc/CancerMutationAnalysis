\name{GeneAlterGBM}
\docType{data}
\alias{GeneAlterGBM}

\title{Data from the Parsons et al. 2008 study: Alterations for
  every gene and sample}
\description{Somatic alterations for each gene and tumor sample from the
  Parsons et al. 2008 glioblastoma multiforme (GBM) study.}

\usage{data(ParsonsGBM08)}
\format{
  The somatic mutations in the GBM study from Parsons et al., broken
  down by \emph{gene}, \emph{type} (point mutation, amplification,
  or deletion), \emph{sample},
  \emph{screen} (Discovery or Prevalence),
  and, for point mutations, \emph{mutation type},
  composed of the wild type nucleotide, its context,
  and the mutated nucleotide. The object is a data frame,
  with the variables: Gene, Type, Sample, Screen, WTNuc
  (wild type nucleotide), Context, and MutNuc (mutated
  nucleotide). The two possible values
  for Screen are Disc ("Discovery") and Prev ("Prevalence").
  For this study, only the Discovery screen is considered.
  The three possible values for Type are Mut (point mutations),
  Amp (large amplifications), and Del (large deletions.)
  Indels have a "" entry for WTNuc, an "All"
  entry for Context, and a "ins.del" entry for MutNuc.
  Large amplifications and deletions have "" entries for
  WTNuc, Context, and MutNuc.
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
  \code{GeneSampGBM}
}

\keyword{datasets}
