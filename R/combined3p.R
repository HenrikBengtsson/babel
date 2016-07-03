combined3.p <- function(vec)
  {
    if(sum(is.na(vec)) > 0) output <- NA
    else
      {
        n <- length(vec)
        n.fac <- gamma(n+1)
        sum.vec <- sum(vec)
        floor.vec <- floor(sum.vec)
        output <- 0
        for(k in 0:floor.vec)
          {
            output <- output+((-1)^k)*(n.fac)*((sum.vec-k)^n)/(gamma(n-k+1)*gamma(k+
                                                                                  1))
#        print(output)
          }
        output <- output/n.fac
      }
    output
  }

