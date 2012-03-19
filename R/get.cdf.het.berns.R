get.cdf.het.berns <-
function(n, j, theta.n, alpha.n)
{  
  ##sort the alphas
  alpha.n <- alpha.n[order(alpha.n)]
  theta.n <- theta.n[order(alpha.n)]

  ##call the C code
  returnvalue <- .C("cdf",as.integer(n),as.double(j),
                    as.double(theta.n), as.double(alpha.n),
                    result=double(1),
                    PACKAGE="PatientGeneSets")
  cdf <- returnvalue$result
  
  cdf    
}

