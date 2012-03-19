\name{GeneCovPancreas}
\docType{data}
\alias{GeneCovPancreas}

\title{Data from the Jones et al. 2008 study: Total number of
  nucleotides "at risk" ("coverage")}
\description{Total numbers of nucleotides "at risk" that were 
successfully sequenced in RefSeq genes in the Jones et al. 2008
pancreatic cancer study.}

\usage{data(JonesPancreas08)}
\format{
  Total number of nucleotides available for mutations ("coverage")
  in the pancreatic cancer study from Jones et al., 
  broken down by \emph{gene}, \emph{screen} (Discovery or Prevalence),
  and \emph{mutation type}, composed of the wild type nucleotide
  and its context. The object is a data frame, with the
  variables: Gene, Screen, WtNuc (wild type nucleotide),
  Context, and Coverage. The two possible values
  for Screen are Disc ("Discovery") and Prev ("Prevalence").
  For this study, only the Discovery screen is considered.
  The nucleotides availables for indels are all the successfully
  sequenced nucleotides in a gene; the corresponding rows
  have a "" entry for WtNuc and an "All" entry for "Context."
  The nucleotides availables for other
  mutations are excluding nucleotides which can only give rise to
  synonymous mutations. 
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
  \code{GeneSampPancreas}
}

\keyword{datasets}
