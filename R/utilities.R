#Turn one-sided p-values into two-sided p-values
one.to.two <- function(p)
  {
    2*min(p,1-p)
  }

#group is a vector of group labels
#lib.size is the effective library size

estDisp <- function(counts,group,lib.size)
  {
    estimateCommonDisp(DGEList(counts=counts,group=group,lib.size=lib.size))$common.dispersion
  }    

