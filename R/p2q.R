p2q <-
function(pvalue,BH=TRUE){
  pvalue[is.na(pvalue)] <- 1
  if (BH) {
    qval.out <- ( length(pvalue) * pvalue ) / rank(pvalue)
    qval.out[qval.out>1] <- 1
    }
  if (!BH) qval.out <- qvalue(pvalue)$qvalue
  return(qval.out)
}

