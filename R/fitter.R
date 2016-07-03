fitter <- function(rna.vector, rp.vector, method, trim.x=NULL, trim.y=FALSE, disp=0.2, alpha=0.001)
  {
    if(method=="rlm") fit <- lsfit(rna.vector,rp.vector)
#    if(method=="rlm") fit <- rlm(rna.vector,rp.vector)
    else
      {
        if(!is.null(trim.x))
          {
            quantile.low <- quantile(rna.vector,trim.x)
            quantile.high <- quantile(rna.vector,1-trim.x)
            keepers.x <- which(rna.vector>=quantile.low & rna.vector<=quantile.high)
            rna.vector <- rna.vector[keepers.x]
            rp.vector <- rp.vector[keepers.x]
          }
        if(trim.y)
          {
            fit1 <- lsfit(rna.vector,rp.vector,intercept=FALSE)
            fitted.values <- rna.vector*fit1$coefficients
            pvalue.vector <- pnbinom(rp.vector,mu=fitted.values,size=1/disp)
            keepers.y <- which(pvalue.vector>alpha & pvalue.vector<(1-alpha))
            rna.vector <- rna.vector[keepers.y]
            rp.vector <- rp.vector[keepers.y]
          }
        if(method=="ls")
          {
            fit <- lsfit(rna.vector,rp.vector,intercept=FALSE)
          }
        else if(method=="blue")
          {
            fit1 <- lsfit(rna.vector,rp.vector,intercept=FALSE)
            fit1.values <- rna.vector*fit1$coefficients
            fit1.phi <- rpDisp(rna.vector,rp.vector)$disp
            fit1.variances <- fit1.values*(1+fit1.values*fit1.phi)
            fit <- lsfit(rna.vector,rp.vector,wt=1/fit1.variances,intercept=FALSE)
          }
      }
    fit
  }

