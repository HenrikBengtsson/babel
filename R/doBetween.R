doBetween <- function(pmat,group,nSD=3,type=c("one-sided","two-sided"))
  {
    zmat <- qnorm(as.matrix(pmat))
    group <- unlist(group)
    unique.group <- unique(group)
    n <- length(unique.group)
    m <- choose(n,2)
    pval <- matrix(NA_real_,nrow(pmat),m)
    direction <- matrix(1,nrow(pmat),m)
    colnames(pval) <- 1:m
    counter <- 1
    for(i in seq_len(n-1))
      {
        which.i <- which(group==unique.group[i])
        zmat.i <- as.matrix(zmat[,which.i])
        num1 <- rowSums(zmat.i)
        which.del.i <- which(is.na(num1))
        for(j in (i+1):n)
          {
            which.j <- which(group==unique.group[j])
            zmat.j <- as.matrix(zmat[,which.j])
            num2 <- rowSums(zmat.j)
            which.del.j <- which(is.na(num2))
            num <- num1-num2
            p <- length(which.i)+length(which.j)
            which.del <- which(abs(num)>(nSD*sqrt(p)))
            which.del <- sort(unique(c(which.del.i,which.del.j,which.del)))
            zmat.ij <- cbind(zmat.i,zmat.j)
            if(length(which.del)>0) zmat.ij <- zmat.ij[-which.del,]
            covmat <- cov(zmat.ij)
            var.ij <- sum(diag(covmat))
            for(ii in seq_len(ncol(zmat.ij)-1))
                {
                  for(jj in (ii+1):ncol(zmat.ij))
                    {
                      if(ii<=length(which.i) & jj>length(which.i))
                        {
                          var.ij <- var.ij-2*covmat[ii,jj]
                        }
                      else
                        {
                          var.ij <- var.ij+2*covmat[ii,jj]
                        }
                    }
                }
            se.ij <- sqrt(var.ij)
            tstat <- num/se.ij
            if(type=="one-sided") pval[,counter] <- pnorm(tstat)
            else if (type=="two-sided")
              {
                pval[,counter] <- 2*(1-pnorm(abs(tstat)))
                which.switch <- which(tstat>0)
                direction[which.switch,counter] <- (-1)
              }
            colnames(pval)[counter] <- paste(unique.group[i],".vs.",unique.group[j],sep="")
            counter <- counter+1
          }
      }
    list(pval=pval,direction=direction)
  }

