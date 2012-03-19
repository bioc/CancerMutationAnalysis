get.EventsBySampleBySet.from.EventsBySample <-
  function(EventsBySample,
           GeneSets,
           GoodSamples)
{
    EventsBySampleBySet <- matrix(0, nrow = length(GeneSets),
                                  ncol = length(GoodSamples))

    colnames(EventsBySampleBySet) <- paste("Events@",GoodSamples,sep="")
    rownames(EventsBySampleBySet) <- names(GeneSets)

    ##split EventsBySample by sample
    split.EventsBySample <- split(EventsBySample, EventsBySample$Sample)
    ##same order as GoodSamples
    split.EventsBySample <- split.EventsBySample[GoodSamples]

    nr.samples <- length(split.EventsBySample)

    for (ss in 1:nr.samples)
    {
        GenesAlteredinSample <- split.EventsBySample[[ss]][,3] # get genes altered in sample ss
        tmp <- lapply(GeneSets, intersect, GenesAlteredinSample) # get genes from each gene-set which are altered in sample ss
        AlteredSamples <- sapply(tmp, length)

        EventsBySampleBySet[, ss] <- AlteredSamples
    }

    EventsBySampleBySet
}
