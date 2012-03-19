runGSEA <-
function(AnnList, data, datacol, IDcol, IDinRowNames=TRUE, absolute=TRUE, type= "f", alternative="mixed", ranks.only=TRUE, nsim=NULL) {
	if (absolute==TRUE) {
          data[,datacol] <- abs(as.numeric(data[,datacol]))
	}
        
        if (IDinRowNames) {
          setID <- NULL
          gseaPval <- NULL
          if (sum(rownames(data)%in%AnnList)==dim(data)[1]) {
            setID <- names(AnnList[1])
            gseaPval <- NA
          }
          else  {
            if (TRUE%in%(rownames(data)%in%AnnList)==TRUE) {
              setID <- names(AnnList[1])
              gseaPval <- geneSetTest((rownames(data)%in%AnnList), as.numeric(data[,datacol]), alternative, type, ranks.only, nsim)
            }
            else {
              setID <- names(AnnList[1])
              gseaPval <- NA
            }
          }
        }
        
        if (!IDinRowNames) {
          setID <- NULL
          gseaPval <- NULL
          if (sum(data[,IDcol]%in%AnnList)==dim(data)[1]) {
            setID <- names(AnnList[1])
            gseaPval <- NA
          }
          else  {
            if (TRUE%in%(data[,IDcol]%in%AnnList)==TRUE) {
              setID <- names(AnnList[1])
              gseaPval <- geneSetTest((data[,IDcol]%in%AnnList), as.numeric(data[,datacol]), alternative, type, ranks.only, nsim)
            }
            else {
              setID <- names(AnnList[1])
              gseaPval <- NA
            }
          }
        }
        
	ans <- as.vector(gseaPval)
	names(ans) <- names(setID)
	return(ans)
      }

