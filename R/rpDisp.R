rpDisp <- function(x,y,nbins=20,min.bin=10,trim.x=0.1,fixdisp=0.2,fixalpha=0.01)
  {
#Sort for later    
    order.x <- order(x)
    x <- x[order.x]
    y <- y[order.x]
#Lower quantile of x
    x.lower <- quantile(x,trim.x)
#Upper quantile of x
    x.upper <- quantile(x,1-trim.x)
    keepers <- which(x>=x.lower & x<=x.upper)
    x <- x[keepers]
    y <- y[keepers]
#Now fit by ls
    fit <- lsfit(x,y,intercept=FALSE)
#    fit <- rlm(x,y)    
    fitted.values <- x*fit$coefficients
#Do two runs, one with a guess of alpha, the other with the estimated alpha
    j <- 1
    while(j<3)
      {
        if(j==1) pvalues <- pnbinom(y,mu=fitted.values,size=1/fixdisp)
        else if (j==2) pvalues <- pnbinom(y,mu=fitted.values,size=1/disp)
        pvalues <- sapply(pvalues,one.to.two)
        keepers <- rep(TRUE,length(pvalues))
        keepers[which(pvalues<fixalpha)] <- FALSE
        x.keepers <- x[keepers]
        y.keepers <- y[keepers]
#    lowessfit <- rlm(x,y,f=f)    
#Starting and ending of bins based on quantiles so same number in every bin
        x.bins <- quantile(x.keepers,seq(0,1,length.out=nbins+1))
#    x.bins <- quantile(x[x>=x.lower & x<=x.upper],seq(0,1,length.out=nbins+1))    
        v <- rep(NA_real_,nbins)
        lambda <- rep(NA_real_,nbins)
        count.bins <- rep(NA_integer_,nbins)
        for(i in seq_len(nbins))
          {
#Which x are in the current bin
            which.bin <- which(x.keepers>=x.bins[i] & x.keepers<x.bins[i+1])
#How many are in the current bin
            count.bins[i] <- length(which.bin)
#Unnecessary check for enough in bin
            if(length(which.bin)>=min.bin)
              {
#Which ys are in bin            
                current.y <- y.keepers[which.bin]
##Variance of trimmed y value
#            v[i] <- mad(current.y)^2
                v[i] <- var(current.y)
#            v[i] <- mad(current.y)^2            
#            v[i] <- var(trimmed.y)            
#Which value of sorted x is in middle of bin
                x.middle <- median(x.keepers[which.bin])
#            lowessind <- findIndices(lowessfit$x,(x.bins[i]+x.bins[i+1])/2)            
#What is the y value at the middle bin position, that is the intensity
                lambda[i] <- x.middle*fit$coefficients
#            lambda[i] <- lowessfit$y[lowessind]            
              }
          }
        lambda <- lambda[!is.na(lambda)]
        v <- v[!is.na(v)]
        term1 <- sum(v*lambda*(1+lambda))
        term2 <- sum(lambda^3)
        term3 <- sum(lambda^4)
        disp <- (term1-term2)/term3
        j <- j+1
      }
    list(disp=disp,lambda=lambda,v=v,x.bins=x.bins,count.bins=count.bins,pvalues=sort(pvalues),fitted.values=fitted.values,keepers=keepers,x=x,y=y)
  }

