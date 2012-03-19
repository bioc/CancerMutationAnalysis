\name{extract.sims.method}
\alias{extract.sims.method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts the p-values or q-values from a SetMethodsSims object
for a specific method.
}
\description{
This function is used to obtain a single data frame with the p-values or
q-values from one of the specific gene-set analysis methods, from a
\code{SetMethodsSims} object which has the results from simulated datasets.
}
\usage{
extract.sims.method(object, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Object of the class \code{SetMethodsSims}.
}
  \item{method}{
  Character string giving the method used for extraction, and whether p-values
  or q-values are extracted. The string should be one of the column names
  of the data frame resulting from the \code{cma.set.stat} function.  
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
  \emph{Genome Biology.} DOI:10.1186/gb-2010-11-11-r112
}
\author{
Simina M. Boca, Giovanni Parmigiani.
}

\seealso{
\code{SetMethodsSims-class},
\code{cma.set.sim},
\code{cma.set.stat}
}
\examples{
\dontrun{
##Note that this takes a few minutes to run:
library(KEGG.db)
data(ParsonsGBM08)
data(EntrezID2Name)

setIDs <- c("hsa00250", "hsa05213")
set.seed(831984)
ResultsSim <- 
    sim.data.p.values(cma.alter = GeneAlterGBM,
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

extract.sims.method(ResultsSim, "p.values.perm.null")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
