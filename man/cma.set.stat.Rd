\name{cma.set.stat}
\alias{cma.set.stat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Implements gene-set analysis methods.
}
\description{
This function implements the gene-set analysis methods. It returns a data-frame
with p-values and q-values for all the methods selected.
}
\usage{
cma.set.stat(cma.alter,
                     cma.cov,
                     cma.samp,
                     GeneSets,
                     ID2name=NULL, 
		     Scores,
		     passenger.rates = t(data.frame(0.55*rep(1.0e-6,25))),
		     BH = TRUE,
		     gene.method = FALSE, 
		     perm.null.method = TRUE, 
		     perm.null.het.method = FALSE,
		     pass.null.method = FALSE, 
		     pass.null.het.method = FALSE,
                     score = "logLRT",
                     verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cma.alter}{
    Data frame with somatic mutation information, broken down by
    gene, sample, screen, and mutation type.
    See \code{GeneAlterBreast} for an example.
  }
  \item{cma.cov}{
    Data frame with the total number of nucleotides "at
    risk" ("coverage"), broken down by gene, screen, and mutation type.
    See \code{GeneCovBreast} for an example.
  }
  \item{cma.samp}{
    Data frame with the number of samples analyzed,
    broken down by gene and screen.
    See \code{GeneSampBreast} for an example.
  }
  \item{GeneSets}{
    An object which annotates genes to gene-sets; it can either be a list
    with each component representing a set, or an object of the class
    AnnDbBimap.
  }
  \item{ID2name}{
    Vector mapping the gene identifiers used in the GeneSets object
    to the gene names used in the other objects; if they are the same,
    this parameter is not needed. See \code{EntrezID2Name} for an example.
  }
  \item{Scores}{
    Data frame of gene scores. The logLRT scores are used for the gene.method
    option. It can be the output of \code{cma.scores}. 
    If the gene.method option is set
    to FALSE, this parameter is not needed.
  }
  \item{passenger.rates}{
    Data frame with 1 row and 25 columns, of passenger mutation rates per
    nucleotide, by type, or "context". Columns denote types and must be
    in the same order as the first 25 columns in the 
    \code{MutationsBrain} objects. 
  }
  \item{BH}{
    If set to \code{TRUE}, uses the Benjamini-Hochberg method to get q-values;
    if set to \code{FALSE}, uses the Storey method from the 
    \code{qvalue} package.
  }
  \item{gene.method}{
    If set to \code{TRUE}, implements gene-oriented method.
  }
  \item{perm.null.method}{
    If set to \code{TRUE}, implements patient-oriented method 
    with permutation null and no heterogeneity.
  }
  \item{perm.null.het.method}{
    If set to \code{TRUE}, implements patient-oriented method 
    with permutation null and heterogeneity.
  }
  \item{pass.null.method}{
    If set to \code{TRUE}, implements patient-oriented method 
    with passenger null and no heterogeneity.
  }
  \item{pass.null.het.method}{
    If set to \code{TRUE}, implements patient-oriented method 
    with passenger null and heterogeneity.
  }
  \item{score}{
    Can be any of the scores which result from \code{cma.scores}.
    Specifies the gene-scoring mechanism used in the gene-oriented
    method.
  }
  \item{verbose}{
    If \code{TRUE}, prints intermediate messages.
  }
}
\value{
  A data frame, with the rows representing set names and the columns
  representing the p-values and q-values corresponding to the different
  methods.
}
\references{
  Boca SM, Kinzler KW, Velculescu VE, Vogelstein B,
  Parmigiani G.
  Patient-oriented gene-set analysis for cancer mutation data.
  \emph{Genome Biology.} DOI: 10.1186/gb-2010-11-11-r112
    
  Parmigiani G, Lin J, Boca S, Sjoeblom T, Kinzler KW,
  Velculescu VE, Vogelstein B. Statistical methods for the analysis of
  cancer genome sequencing data. 
  \url{http://www.bepress.com/jhubiostat/paper126/}
  
  Benjamini Y and Hochberg Y. 
  Controlling the false discovery rate: a practical and powerful
  approach to multiple testing.
  \emph{Journal of the Royal Statistical Society B.}
  DOI: 10.2307/2346101 
  
  Storey JD and Tibshirani R.
  Statistical significance for genome-wide experiments.
  \emph{Proceedings of the National Academy of Sciences.}
  DOI: 10.1073/pnas.1530509100

  Schaeffer EM, Marchionni L, Huang Z, Simons B, Blackman A, Yu W,
  Parmigiani G, Berman DM.
  Androgen-induced programs for prostate epithelial growth and invasion arise in embryogenesis and are reactivated in cancer.
  \emph{Oncogene.}
  DOI: 10.1038/onc.2008.327

  Thomas MA, Taub AE. 
  Calculating binomial probabilities when the trial probabilities are
  unequal.
  \emph{Journal of Statistical Computation and Simulation.}
  DOI: 10.1080/00949658208810534
  
  Parsons DW, Jones S, Zhang X, Lin JCH, Leary RJ, Angenendt P, Mankoo P,
  Carter H, Siu I, et al. 
  An Integrated Genomic Analysis of Human Glioblastoma Multiforme. 
  \emph{Science.} DOI: 10.1126/science.1164382

  Wood LD, Parsons DW, Jones S, Lin J, Sjoeblom, Leary RJ, Shen D,
  Boca SM, Barber T, Ptak J, et al. The Genomic Landscapes of Human
  Breast and Colorectal Cancer. \emph{Science.} DOI: 10.1126/science.1145720
}
\author{
  Simina M. Boca, Giovanni Parmigiani, Luigi Marchionni, Michael A. Newton.
}

\seealso{
\code{GeneCov}, \code{GeneSamp}, \code{GeneAlter},
\code{BackRates}, \code{cma.scores}, \code{cma.set.sim}
}

\examples{
library(KEGG.db)
data(ParsonsGBM08)
data(EntrezID2Name)

setIDs <- c("hsa00250", "hsa05213")
SetResults <- cma.set.stat(cma.alter = GeneAlterGBM,
                                   cma.cov = GeneCovGBM,
                                   cma.samp = GeneSampGBM,
                                   GeneSets =  KEGGPATHID2EXTID[setIDs],
                                   ID2name = EntrezID2Name,
                                   perm.null.method = TRUE,
                                   pass.null.method = TRUE)

SetResults
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

