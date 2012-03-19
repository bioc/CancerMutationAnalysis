\name{cma.scores}
\alias{cma.scores}
\title{Gene-level scores for the analysis of somatic point mutations in cancer}
\description{
  Computes various gene-level scores for the analysis of somatic point mutations in cancer.
}
\usage{
cma.scores(cma.alter = NULL,
           cma.cov,
           cma.samp,
           scores = c("CaMP", "logLRT"),
           cma.data = NULL,
           coverage = NULL,
	   passenger.rates = t(data.frame(0.55*rep(1.0e-6,25))),
           allow.separate.rates = TRUE,
           filter.above=0, 
           filter.below=0, 
           filter.threshold=0,
	   filter.mutations=0,
           aa=1e-10, 
           bb=1e-10,
           priorH0=1-300/13020, 
           prior.a0=100,
           prior.a1=5,
           prior.fold=10)}
\arguments{
  \item{cma.alter}{Data frame with somatic mutation information, broken down by
    gene, sample, screen, and mutation type.
    See \code{GeneAlterBreast} for an example.}
  \item{cma.cov}{Data frame with the total number of nucleotides "at
    risk" ("coverage"), broken down by gene, screen, and mutation type.
    See \code{GeneCovBreast} for an example.}
  \item{cma.samp}{Data frame with the number of samples analyzed,
    broken down by gene and screen.
    See \code{GeneSampBreast} for an example.}
  \item{scores}{Vector with the scores which are to be computed.
    It can include: \code{CaMP} (Cancer Mutation Prevalence score),
    \code{logLRT} (log Likelihood Ratio Test score),
    \code{neglogPg}, \code{logLRT}, \code{logitBinomialPosteriorDriver},
    \code{PoissonlogBF}, \code{PoissonPosterior},
    \code{Poissonlmlik0}, \code{Poissonlmlik1}}
  \item{cma.data}{Provided for back-compatibility and internal
    operations. \code{cma.data} and \code{coverage} objects were
    used in prior versons of this package, and may be specified
    instead of \code{cma.alter}, \code{cma.cov}, and \code{cma.samp}.}
  \item{coverage}{Provided for back-compatibility and internal
    operations. \code{cma.data} and \code{coverage} objects were
    used in prior versons of this package, and may be specified
    instead of \code{cma.alter}, \code{cma.cov}, and \code{cma.samp}.}
  \item{passenger.rates}{Data frame of "passenger" (or "background")
    mutation rates per nucleotide, by type, or "context". If two
    rows are present, the first refers to the Discovery screen and
    the second to the Prevalence screen.}
  \item{allow.separate.rates}{If \code{TRUE}, allows for use separate rates for
    Discovery and Prevalence screens.}
  \item{filter.threshold}{This and the following three input control
    filtering of genes, allowing to exclude genes from analysis, by size
    and number of mutations. Different criteria can be set above and
    below this threshold. The threshold is a gene size in base pairs.}
  \item{filter.above}{Minimum number of mutations per
    Mb, applied to genes of size greater than \code{threshold.size}.}
  \item{filter.below}{Minimum number of mutations per
    Mb, applied to genes of size lower than \code{threshold.size}.}
  \item{filter.mutations}{Only consider genes
    whose total number of mutations is greater than or equal to
    \code{filter.mutations}.}   
  \item{aa}{Hyperparameter of beta prior used in compute.binomial.posterior.}
  \item{bb}{Hyperparameter of beta prior used in compute.binomial.posterior}
  \item{priorH0}{Prior probability of the null hypothesis, used to
    convert the BF in compute.poisson.BF to a posterior probability}
  \item{prior.a0}{Shape hyperparameter of gamma prior on passenger rates used in compute.poisson.BF}
  \item{prior.a1}{Shape hyperparameter of gamma prior on non-passenger rates used in compute.poisson.BF}
  \item{prior.fold}{Hyperparameter of gamma prior  on non-passenger
    rates used compute.poisson.BF. The mean of the gamma is set so that
    the ratio of the mean to the passenger rate is the specified
    \code{prior.fold} in each type.}
}

\details{  
  The scores computed by this function are relevant for two stage
  experiments like the one in the Sjoeblom et al. article. In this design genes
  are sequenced in a first "Discovery" sample. A non-random set of genes
  is then also sequenced in a subsequent "Prevalence"
  (or "Validation") screen. For instance, in Sjoeblom et al. and
  Wood et al., genes "pass" the Discovery screen if they are mutated
  at least once in it.
  The goal of this tool is to facilitate reanalysis of the
  Sjoeblom et al. 2006, Wood et al. 2007, Jones et al. 2008,
  Parsons et al. 2008, and Parsons et al. 2011 datasets.
  Application to other projects requires a detailed
  understanding of these projects.
}
\value{
  A data frame giving gene-by-gene values for each score. The columns in
  this data frame are:
  \item{CaMP}{The CaMP score of Sjoeblom and colleagues.}
  \item{neglogPg}{The negative log10 of Pg, where Pg represents the probability
  that a gene has its exact observed mutation profile under the null,
  i.e. assuming the given passenger rates.}
  \item{logLRT}{The log10 of the likelihood ratio test (LRT).}
  \item{logitBinomialPosteriorDriver}{logit of the posterior
    probability that a gene's  mutation rates above the
    specified passenger rates using a binomial model}
  \item{PoissonlogBF}{The log10 of the Bayes Factor (BF) using a
    Poisson-Gamma model.}
  \item{PoissonPosterior}{The posterior probability that a given
  gene is a driver, using a Poisson-Gamma model.}
  \item{Poissonlmlik0}{Marginal likelihood under the null hypothesis in
    the Poisson-Gamma model}
  \item{Poissonlmlik1}{Marginal likelihood under the alternative  hypothesis in
    the Poisson-Gamma model}
}
\references{
  Parmigiani G, Lin J, Boca S, Sjoeblom T, Kinzler KW,
  Velculescu VE, Vogelstein B. Statistical methods for the analysis of
  cancer genome sequencing data. 
  \url{http://www.bepress.com/jhubiostat/paper126/}

  Sjoeblom T, Jones S, Wood LD, Parsons DW, Lin J, Barber T,
  Mandelker D, Leary R, Ptak J, Silliman N, et al.  The
  consensus coding sequences of breast and colorectal cancers.
  \emph{Science.} DOI: 10.1126/science.1133427

  Wood LD, Parsons DW, Jones S, Lin J, Sjoeblom, Leary RJ, Shen D,
  Boca SM, Barber T, Ptak J, et al. The Genomic Landscapes of Human
  Breast and Colorectal Cancer. \emph{Science.} DOI:
  10.1126/science.1145720

  Parsons DW, Jones S, Zhang X, Lin JCH, Leary RJ, Angenendt P, Mankoo P,
  Carter H, Siu I, et al. 
  An Integrated Genomic Analysis of Human Glioblastoma Multiforme. 
  \emph{Science.} DOI: 10.1126/science.1164382
  
  Jones S, Zhang X, Parsons DW, Lin JC, Leary RJ, Angenendt P, Mankoo P,
  Carter H, Kamiyama H, Jimeno A, et al.
  Core signaling pathways in human pancreatic cancers revealed by global
  genomic analyses.
  \emph{Science.} DOI: 10.1126/science.1164368
  
  Parsons DW, Li M, Zhang X, Jones S, Leary RJ, Lin J, Boca SM, Carter
  H, Samayoa J, Bettegowda C, et al.
  The genetic landscape of the childhood cancer medulloblastoma.
  \emph{Science.} DOI: 10.1126/science.1198056
}
\author{Giovanni Parmigiani, Simina M. Boca}
\seealso{\code{GeneCov}, \code{GeneSamp}, \code{GeneAlter},
  \code{BackRates},
  \code{cma.set.stat}}
\examples{
data(ParsonsGBM08)
ScoresGBM <- cma.scores(cma.alter = GeneAlterGBM,
                        cma.cov = GeneCovGBM,
                        cma.samp = GeneSampGBM)
}
\keyword{htest}

