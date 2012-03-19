\name{cma.fdr}
\alias{cma.fdr}
\title{Gene-level Empirical Bayes (EB) false discovery rate (FDR)
  analysis for somatic mutations in cancer}
\description{
  Empirical Bayes estimates of the False Discovery Rate (FDR) and
  passenger probabilities in the analysis of somatic mutations in cancer.
}
\usage{
cma.fdr(cma.alter,
        cma.cov,
        cma.samp,
        scores = c("CaMP", "logLRT"),
        passenger.rates = t(data.frame(.55*rep(1.0e-6,25))),
        allgenes=TRUE,
        estimate.p0=FALSE,
        p0.step=1,
        p0=1,
        eliminate.noval=FALSE,
        filter.threshold=0, 
        filter.above=0,
        filter.below=0,
        filter.mutations=0, 
        aa=1e-10, 
        bb=1e-10,
        priorH0=1-500/13020, 
        prior.a0=100,
        prior.a1=5,
        prior.fold=10,
        M=2,
        DiscOnly=FALSE,
        PrevSamp="Sjoeblom06",
        KnownCANGenes=NULL,
        showFigure=FALSE,
        cutoffFdr=0.1)}
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
  \item{passenger.rates}{Data frame of passenger mutation rates per
    nucleotide, by type, or "context". If two
    rows are present, the first refers to the Discovery screen and
    the second to the Prevalence screen.}
  \item{allgenes}{If TRUE, genes where no mutations were found are
    considered in the analysis.}
  \item{estimate.p0}{If TRUE, estimates the percent of genes with only
    passenger mutations. Requires \code{allgenes=TRUE} }
  \item{p0.step}{Size of bins of histograms in the distribution of
    scores, to use in estimating p0 if estimate.p0 = TRUE. All scores
    are in the log 10 scale.}
  \item{p0}{Proportion of genes with only
    passenger mutations. Only used if \code{estimate.p0=FALSE}}
  \item{eliminate.noval}{If \code{TRUE}, the genes which are not validated are
    eliminated from the analysis. Validated genes are those where at
    least one mutation was found in both the Discovery and Prevalence
    (or Validation) screens.}
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
  \item{bb}{Hyperparameter of beta prior used in compute.binomial.posterior.}
  \item{priorH0}{Prior probability of the null hypothesis, used to
    convert the BF in compute.poisson.BF to a posterior probability}
  \item{prior.a0}{Shape hyperparameter of gamma prior on passenger rates used in compute.poisson.BF}
  \item{prior.a1}{Shape hyperparameter of gamma prior on non-passenger rates used in compute.poisson.BF}
  \item{prior.fold}{Hyperparameter of gamma prior on non-passenger
    rates used compute.poisson.BF. The mean of the gamma is set so that
    the ratio of the mean to the passenger rate is the specified
    \code{prior.fold} in each type.}
  \item{M}{The number of null
    datasets generated  to get the false discovery rates.
    Numbers on the order of 100 are recommended, but this will
    cause the function to run very slowly.}
  \item{DiscOnly}{If TRUE, only considers data from Discovery screen.}
  \item{PrevSamp}{If "Sjoeblom06", then the experimental design from
    Sjoeblom et al. or Wood et al. is used, namely, genes "pass" from the
    Discovery into the Prevalence (or Validation) screens if they are
    mutated at least once in the Discovery samples. If "Parsons11",
    the experimental design from Parsons et al. 2011 is approximated,
    namely, in the null datasets, a gene passes into the Prevalence
    screen if it is mutated at
    least once, and is found on a specified list of known cancer candidate
    (CAN) genes, or if it is mutated at least twice.}
  \item{KnownCANGenes}{Vector of known CAN genes, to be used if
    PrevSamp is not set to "Sjoeblom07".}
  \item{showFigure}{If TRUE, displays a figure for each score in
    \code{scores}, showing the right tail of the density of scores
    under the null, the right tail of the density of real scores
    as a rug (1-d) plot and the number of real genes and average number
    of null genes to the right of the cutoff chosen based on
    \code{cutoffFdr.}}
  \item{cutoffFdr}{If \code{showFigure} is set to TRUE, it gives the
     value at which we are interested in controlling the false discovery
     rate (Fdr). The corresponding score threshold is plotted on the figure,
     with the number of real genes greater than it and the average
     number of null genes greater than it specified. The estimated Fdr
     at that threshold is the ratio of the average number of null
     genes and the number of real genes, multiplied by \code{p0}, which
     is often taken to be 1.}
}
\value{
  A list of data frames. Each gives a gene gene-by-gene significance for
  one of the score requested. The columns in
  each data frame are:
  \item{score}{The score requested (e.g. the LRT).}
  \item{F}{Number of genes experimentally observed to give a larger
    score than the gene in question.}
  \item{F0}{Number of genes giving a larger
    score than the gene in question in datasets simulated from passenger
    mutation rates.}
  \item{Fdr}{The Empirical Bayes False Discovery Rate, as defined in
    Efron and Tibshirani 2002.}
  \item{fdr}{The Empirical Bayes Local False Discovery Rate, as defined in
    Efron and Tibshirani 2002.}
  \item{p0}{Scalar, Proportion of genes with only
    passenger mutations. Estimated or passed on from input (depending on
    whether estimate.p0 is \code{TRUE}} 
}
\references{
  Efron B, Tibshirani R.  Empirical {B}ayes methods and false discovery
  rates for microarrays.
  \emph{Genetic Epidemiology}. DOI: 10.1002/gepi.1124
  
  Parmigiani G, Lin J, Boca S, Sjoeblom T, Kinzler KW,
  Velculescu VE, Vogelstein B. Statistical methods for the analysis of
  cancer genome sequencing data, 2007. 
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
  \code{BackRates},\code{cma.scores}}
\examples{
data(ParsonsMB11)
set.seed(188310)
cma.fdr.out <- cma.fdr(cma.alter = GeneAlterMB,
                       cma.cov = GeneCovMB,
                       cma.samp = GeneSampMB,
                       allgenes = TRUE,
                       estimate.p0=FALSE,
                       eliminate.noval=FALSE,
                       filter.mutations=0,
                       M = 2)
names(cma.fdr.out)
}
\keyword{htest}
