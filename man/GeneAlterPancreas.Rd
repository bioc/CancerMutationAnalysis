\name{GeneAlterPancreas}
\docType{data}
\alias{GeneAlterPancreas}

\title{Data from the Jones et al. 2008 study: Alterations for
  every gene and sample}
\description{Somatic alterations for each gene and tumor sample from the
  Jones et al. 2008 pancreatic cancer study.}

\usage{data(JonesPancreas08)}
\format{
  The somatic mutations in the pancreatic cancer study from Jones et
  al., broken
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
  \code{GeneCovPancreas}, 
  \code{GeneSampPancreas}
}

\keyword{datasets}
