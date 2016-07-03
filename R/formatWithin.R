formatWithin <- function(within,method,p,rnames,cnames,keeper.genes,n)
  {
    twosided <- 2*within
    alt <- 2*(1-within)
    which.switch <- which(alt<twosided)
    twosided[which.switch] <- alt[which.switch]
    qvalues <- matrix(NA_real_,n,p)
    for(i in seq_len(p))
      {
        current.keepers <- which(keeper.genes[,i])
        qvalues[current.keepers,i] <- p.adjust(twosided[current.keepers,i],method=method)
      }
#    qvalues <- apply(twosided,2,p.adjust,method=method)
    direction <- matrix(1,nrow(within),ncol(within))
    direction[which(within>0.5)] <- (-1)
    output.within <- vector("list",p)
    names(output.within) <- cnames
    for(i in seq_len(p))
      {
        current.keepers <- which(keeper.genes[,i])
        new.direction <- new.within <- new.twosided <- new.qvalues <- rep(NA_real_,n)
        new.direction[current.keepers] <- direction[current.keepers,i]
        new.within[current.keepers] <- within[current.keepers,i]
        new.twosided[current.keepers] <- twosided[current.keepers,i]
        new.qvalues[current.keepers] <- qvalues[current.keepers,i]
        output.within[[i]] <- cbind.data.frame(rnames,new.direction,new.within,new.twosided,new.qvalues)
#        output.within[[i]] <- cbind.data.frame(rnames,direction[,i],within[,i],twosided[,i],qvalues[,i])        
        rownames(output.within[[i]]) <- NULL
        colnames(output.within[[i]]) <- c("Gene","Direction","P-value (one-sided)","P-value (two-sided)","FDR")        
      }
    output.within
  }

