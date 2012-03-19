\name{combine.sims}
\alias{combine.sims}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Combines two SetMethodSims objects.
}
\description{
  This function is used to combine two \code{SetMethodSims} objects,
  which have the results from simulated datasets, provided that
  the values for \code{pass.null, perc.samples,} and \code{spiked.set.sizes} 
  match up when the objects are generated with the
  \code{sim.data.p.values} function.
}
\usage{
combine.sims(obj1, obj2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj1}{
    Object of the class \code{SetMethodsSims}.
  }
  \item{obj2}{
    Object of the class \code{SetMethodsSims}.
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
\code{cma.set.sim}
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

combine.sims(ResultsSim, ResultsSim)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
