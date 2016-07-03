doCount <- function(rnv.pos,rnadisp,fit,rpdisp,match.rnv.urnv,rpv)
  {
    length.pos <- length(rnv.pos)
    len <- 10000*length.pos
    nbsample <- rnbinom(len,mu=rnv.pos,size=1/rnadisp)
    mu <- (nbsample*fit$coefficients)
    ref <- (matrix(rnbinom(len,mu=mu,size=1/rpdisp),nrow=length.pos,ncol=10000))[match.rnv.urnv,]
    count <- rowSums(ref>rpv)
    count
  }
    
