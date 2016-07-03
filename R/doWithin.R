doWithin <- function(rna,rp,group,nreps,keeper.genes,trim.x=0.1,trim.var=0.1,rnadisp=NULL)
{
  n <- ncol(rna)
#  rna <- rna[keeper.genes,]
#  rp <- rp[keeper.genes,]
  smat <- matrix(NA_real_,nrow(rna),n)
  rownames(smat) <- rownames(rna)
  colnames(smat) <- colnames(rna)
  if(!is.null(rnadisp)) rnadisp <- rnadisp
  else
    {
      unique.group <- unique(group)
      count.group <- rep(NA_integer_,length(unique.group))
      for(i in seq_len(length(unique.group))) count.group[i] <- length(which(group==unique.group[i]))
      keeper.groups <- which(count.group>1)
      if(length(keeper.groups)<2)
         {
           warning("Must be at least two sanmples in at least two groups to estimate rna dispersion.  Fixing at 0.1.  May want to select value for rnadisp in argument to babel")
           rnadisp <- 0.1
         }
      else
         {
           keeper.samples <- which(!is.na(match(group,unique.group[keeper.groups])))
           rnadisp <- estDisp(rna[,keeper.samples],group[keeper.samples],NULL)
         }
    }
#
# Permute
#
  loops <- nreps%/%10000
  for(i in seq_len(n))
    {
      print(paste("Running",colnames(rna)[i]))
      current.keepers <- which(keeper.genes[,i])
      rnv <- rna[current.keepers,i]; # This sample's RNA
      rpv <- rp[current.keepers,i]; # This sample's Ribosome
      rpdisp <- rpDisp(rnv,rpv)$disp; # Defaults for rpDisp currently match our desired input
      urnv <- unique(rnv); # Unique values of RNA read counts
      match.rnv.urnv <- match(rnv,urnv)
      pos <- match(urnv,rnv); # Position of unique RNA read count values
      fit <- fitter(rnv,rpv,method="ls",trim.x=trim.x,trim.y=FALSE)
      rnv.pos <- rnv[pos]
      counts <- mclapply(1:loops, FUN=function(j) {
        doCount(rnv.pos=rnv.pos, rnadisp=rnadisp, fit=fit,                   
                rpdisp=rpdisp, match.rnv.urnv=match.rnv.urnv, rpv=rpv)
      })
          count <- Reduce(`+`,counts)
          smat[current.keepers,i] <- count
    }
  pmat <- (smat+1)/(nreps+2)
  pmat
}

