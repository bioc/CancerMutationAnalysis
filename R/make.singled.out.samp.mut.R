make.singled.out.samp.mut <-
function(EventsBySample,
                                      Mutations,
                                      singled.out.samp)
  {
    ##create mutation object for samples which are singled out
    singled.out.Mutations <- Mutations

    ##first set all mutations to 0
    singled.out.Mutations[,1:50] <- 0

    singled.out.Mutations[,60:77] <-
      Mutations[,60:77]*length(singled.out.samp)/length(unique(EventsBySample$Sample))

    ##set number of tumors analyzed in discovery phase to number of tumors which are singled out
    singled.out.Mutations[,78] <- length(singled.out.samp)

    ##get shortened gene names
    GeneSymb <-
      shorten.gene.names(rownames(Mutations))

    names(GeneSymb) <- rownames(Mutations)

    ##change mutation class column in EventsBySample object to have
    ##the same form as in Mutations object
    mut.context <- as.character(unique(EventsBySample$MutationClass))

    names(mut.context) <- mut.context
    mut.context <- paste("MutationsDiscovery.",mut.context,sep="")
    names(mut.context) <- as.character(unique(EventsBySample$MutationClass))

    singled.out.EventsBySample <-
      EventsBySample[EventsBySample$Sample %in% singled.out.samp,]

    ##get all the alterations in the singled-out samples
    ##for(gene in rownames(singled.out.Mutations))
    for(gene in as.character(singled.out.EventsBySample$Symbol))
      {
        ##get the gene transcripts corresponding to the symbol "gene"
        transcripts <- names(GeneSymb[GeneSymb == gene])
        if(length(transcripts) == 0)
          {
            print(gene)
          }
        for(context in names(mut.context))
          {
            ##get the number of alterations of gene "gene" in context
            ##"context" in the singled-out samples
            nr.alt <-
              sum(singled.out.EventsBySample$Symbol == gene &
                  singled.out.EventsBySample$MutationClass == context)

            ##put that numer of alterations in the correct spot in the
            ##Mutations object
            singled.out.Mutations[transcripts, mut.context[context]] <-
              nr.alt
          }
      }

    rownames(singled.out.Mutations) <- rownames(Mutations)

    ##eliminate all the rows which have no mutations
    singled.out.Mutations <-
      singled.out.Mutations[rowSums(singled.out.Mutations[,1:50])!=0,]

    ##return new Mutations object
    singled.out.Mutations
  }

