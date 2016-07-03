doCombined <- function(pmat,group)
{
  unique.group <- unique(group)
  n <- length(unique.group)
  cmat <- matrix(NA_real_,nrow(pmat),n)
  for(i in seq_len(n))
    {
      which.i <- which(group==unique.group[i])
      p.i <- as.matrix(pmat[,which.i])
      cmat[,i] <- apply(p.i,1,combined3.p)
    }
  colnames(cmat) <- paste("combined.",unique.group,sep="")
  cmat
}

