##function to take whatever follows "@@" off gene names
shorten.gene.names <- function(GeneNames)
  {
    GeneSymb <- GeneNames
    shorten.name <- function(name)
    {
        gsub("(\\w*)@@(\\w*)", "\\1", name, perl=TRUE)
    }
    GeneSymb <- sapply(GeneSymb, shorten.name)
    names(GeneSymb) <- NULL

    GeneSymb
  }



