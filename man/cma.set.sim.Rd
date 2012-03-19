\name{cma.set.sim}
\alias{cma.set.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates data and performs gene-set analysis methods on the simulated datasets.
}
\description{
This function simulates data under the passenger or permutation null, either
under the null or including spiked-in gene-sets. It then calculates the
p-values and q-values for all the selected gene-set analysis methods.
}
\usage{
cma.set.sim(cma.alter,
                  cma.cov,
                  cma.samp,
		  GeneSets, 
		  passenger.rates = t(data.frame(0.55*rep(1.0e-6,25))),
		  ID2name=NULL,
		  BH = TRUE, 
		  nr.iter, 
		  pass.null = FALSE,
		  perc.samples = NULL, 
		  spiked.set.sizes = NULL,
		  gene.method = FALSE, 
		  perm.null.method = TRUE,
		  perm.null.het.method = FALSE, 
		  pass.null.method = FALSE,
		  pass.null.het.method = FALSE, 
		  show.iter,
		  KnownMountains = c("EGFR","SMAD4","KRAS",
		  "TP53","CDKN2A","MYC","MYCN","PTEN","RB1"),
		  exclude.mountains=TRUE,
                  verbose=TRUE)
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
  \item{passenger.rates}{
    Data frame with 1 row and 25 columns, of passenger mutation rates per
    nucleotide, by type, or "context". Columns denote types and must be
    in the same order as the first 25 columns in the 
    \code{MutationsBrain} objects. 
  }
  \item{ID2name}{
    Vector mapping the gene identifiers used in the GeneSets object
    to the gene names used in the other objects; if they are the same,
    this parameter is not needed. See \code{EntrezID2Name} for an example.
  }
  \item{BH}{
    If set to \code{TRUE}, uses the Benjamini-Hochberg method to get q-values;
    if set to \code{FALSE}, uses the Storey method from the 
    \code{qvalue} package.
  }
  \item{nr.iter}{
    The number of iterations to be simulated.
  }
  \item{pass.null}{
    If set to true \code{TRUE}, implements the passenger null hypothesis,
    using the rates from \code{passenger.rates}; otherwise, implements the
    permutation null, permuting mutational events. 
  }
  \item{perc.samples}{
    Vector representing the probabilities of the spiked-in gene-sets being
    altered in any given sample, as percentages; for example
    \code{perc.samples = c(75, 90)} means that these probabilities are
    \code{0.75} and \code{0.90}.
  }
  \item{spiked.set.sizes}{
    Vector representing the sizes, in genes, of the spiked-in gene-sets;
    for example, if \code{perc.samples = c(75, 90)} and
    \code{spiked.set.sizes = c(50, 100)}, there would be 4 spiked-in sets,
    one with 50 genes and probability of being altered of 0.75 in each sample,
    one with 50 genes and probability of being altered of 0.90 in each sample,
    one with 100 genes and probability of being altered of 0.75 in each sample, and
    one with 100 genes and probability of being altered of 0.90 in each sample.
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
  \item{show.iter}{
    If set to \code{TRUE} and \code{verbose} is also set to \code{TRUE}, shows what simulation is currently running.
  }
  \item{KnownMountains}{
    Vector of genes to be excluded from the permutation null simulations 
    if \code{exclude.mountains = TRUE}. 
  }
  \item{exclude.mountains}{
    If set to \code{TRUE}, excludes the genes in \code{KnownMountains.}
  }
  \item{verbose}{
    If \code{TRUE}, prints intermediate messages.
  }
}
\value{
An object of the class \code{SetMethodsSims}. See
\code{SetMethodsSims} for more details.
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
  \emph{Journal of the Royal Statistical Society B},
  DOI: 10.2307/2346101 

  Storey JD and Tibshirani R.
  Statistical significance for genome-wide experimens.
  \emph{Proceedings of the National Academy of Sciences.}
  DOI: 10.1073/pnas.1530509100
  
  Parsons DW, Jones S, Zhang X, Lin JCH, Leary RJ, Angenendt P, Mankoo P,
  Carter H, Siu I, et al. 
  An Integrated Genomic Analysis of Human Glioblastoma Multiforme. 
  \emph{Science.} DOI: 10.1126/science.1164382
  
  Wood LD, Parsons DW, Jones S, Lin J, Sjoeblom, Leary RJ, Shen D,
  Boca SM, Barber T, Ptak J, et al. The Genomic Landscapes of Human
  Breast and Colorectal Cancer. \emph{Science.} DOI: 10.1126/science.1145720
}
\author{
Simina M. Boca, Giovanni Parmigiani.
}

\seealso{
\code{SetMethodsSims-class},
\code{CoverageBrain},
\code{EventsBySampleBrain}, \code{GeneSizes08},
\code{MutationsBrain}, \code{ID2name},
\code{cma.set.stat},
\code{extract.sims.method}, 
\code{combine.sims}
}
\examples{
##Note that this takes a few minutes to run:
library(KEGG.db)
data(ParsonsGBM08)
data(EntrezID2Name)

setIDs <- c("hsa00250", "hsa05213")
set.seed(831984)
ResultsSim <- 
    cma.set.sim(cma.alter = GeneAlterGBM,
                      cma.cov = GeneCovGBM,
                      cma.samp = GeneSampGBM,
                      GeneSets =  KEGGPATHID2EXTID[setIDs],
                      ID2name = EntrezID2Name,
                      nr.iter = 2,
                      pass.null = TRUE,
                      perc.samples = c(75, 95),
                      spiked.set.sizes = 50,
                      perm.null.method = TRUE,
                      pass.null.method = TRUE)

ResultsSim
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{datagen}

