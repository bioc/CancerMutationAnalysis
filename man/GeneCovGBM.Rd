\name{GeneCovGBM}
\docType{data}
\alias{GeneCovGBM}

\title{Data from the Parsons et al. 2008 study: Total number of
  nucleotides "at risk" ("coverage")}
\description{Total numbers of nucleotides "at risk" that were 
successfully sequenced in RefSeq genes in the Parsons et al. 2008
glioblastoma multiforme (GBM) study.}

\usage{data(ParsonsGBM08)}
\format{
  Total number of nucleotides available for mutations ("coverage")
  in the GBM study from Parsons et al., 
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
  \code{GeneAlterGBM}, 
  \code{GeneSampGBM}
}

\keyword{datasets}
