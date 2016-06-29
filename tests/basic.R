library(babel)
data(ribo.prof)
test.rna <- ribo.prof$test.rna
test.rp <- ribo.prof$test.rp
test.group=c("A","B","A","B")
set.seed(12345)
options(mc.cores=1)
options(babel.minreps=10000)
test.babel <- babel(test.rna,test.rp,group=test.group,nreps=10000,min.rna=10)
within.babel <- test.babel$within
combined.babel <- test.babel$combined
between.babel <- test.babel$between
length.within <- length(within.babel)
length.combined <- length(combined.babel)
length.between <- length(between.babel)
if(length.within!=4) stop(paste("within is length",length.within,",it should be length 4",sep=""))
if(length.combined!=2) stop(paste("combined is length",length.combined,",it should be length 2",sep=""))
if(length.between!=1) stop(paste("between is length",length.between,",it should be length 1",sep=""))
                           

