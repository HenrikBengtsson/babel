babel <- function(rna,rp,group,nreps,method.adjust="BH",min.rna=10,nSD=3,...)
  {
    if(sum(dim(rna)==dim(rp))<2) stop("rna and rp are different sizes")
    n <- nrow(rna)
    p <- ncol(rp)
    if(length(group)!=p) stop("group must have same length as the number of columns of rna")
    if(is.null(rownames(rna))) rownames(rna) <- 1:n
    if(is.null(rownames(rp))) rownames(rp) <- 1:n
    if(is.null(colnames(rna))) colnames(rna) <- 1:p
    if(is.null(colnames(rp))) colnames(rp) <- 1:p
    if(sum(rownames(rna)==rownames(rp))<n) stop("rownames of rna and rp must match")
    if(sum(colnames(rna)==colnames(rp))<p) stop("colnames of rna and rp must match")
    if((nreps%%10000)!=0) stop("nreps must be divisible by 10000")
    minreps <- getOption("babel.minreps", 100000L)
    if(nreps<minreps) stop(paste("nreps must at least",minreps))
    if(min.rna<1) stop("min.rna needs to be at least 1")
#    if(length(unique(group))!=2) stop("There must be exactly two groups")
#    mins.rna <- apply(rna,1,min)
    keeper.genes <- matrix(FALSE,nrow(rna),ncol(rna))
    keeper.genes[which(rna>=min.rna)] <- TRUE
    within <- doWithin(rna=rna,rp=rp,group=group,nreps=nreps,keeper.genes=keeper.genes,...)
    output.within <- formatWithin(within,method=method.adjust,p=p,rnames=rownames(rna),cnames=colnames(rna),keeper.genes=keeper.genes,n=n)
    combined <- doCombined(within,group)
    output.combined <- formatCombined(combined,group=group,method=method.adjust,rnames=rownames(rna),n=n)
    between <- doBetween(within,group,nSD=nSD,type="two-sided")
    output.between <- formatBetween(between,rna=rna,group=group,method=method.adjust,n=n)
    list(within=output.within,combined=output.combined,between=output.between)
  }

