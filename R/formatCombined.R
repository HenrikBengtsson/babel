formatCombined <- function(combined,group,method,rnames,n)
{
  twosided <- 2*combined
  alt <- 2*(1-combined)
  which.switch <- which(alt<twosided)
  twosided[which.switch] <- alt[which.switch]
  qvalues <- apply(twosided,2,p.adjust,method=method)
  direction <- matrix(1,nrow(combined),ncol(combined))
  direction[which(combined>0.5)] <- (-1)
  ugroup <- unique(group)
  lgroup <- length(ugroup)
  output.combined <- vector("list",lgroup)
  names(output.combined) <- ugroup
  for(i in seq_len(lgroup))
    {
      current.keepers <- which(!is.na(combined[,i]))
      new.direction <- new.twosided <- new.qvalues <- rep(NA_real_,n)
      new.direction[current.keepers] <- direction[current.keepers,i]
      new.twosided[current.keepers] <- twosided[current.keepers,i]
      new.qvalues[current.keepers] <- qvalues[current.keepers,i]
      output.combined[[i]] <- cbind.data.frame(rnames,new.direction,new.twosided,new.qvalues)
#        output.combined[[i]] <- cbind.data.frame(rnames,direction[,i],twosided[,i],qvalues[,i])
      rownames(output.combined[[i]]) <- NULL
      colnames(output.combined[[i]]) <- c("Gene","Direction","P-value","FDR")        
    }
  output.combined
}

