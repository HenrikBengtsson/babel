formatBetween <- function(between,rna,group,method,n=n)
  {
    min.group <- min(table(group))
    ugroup <- unique(group)
    lgroup <- length(ugroup)
    m <- choose(lgroup,2)
    output.between <- vector("list",m)
    counter <- 1
    for(i in seq_len(lgroup-1))
      {
        for(j in (i+1):lgroup)
          {
            names(output.between)[counter] <- paste(ugroup[i],".vs.",ugroup[j],sep="")
            qvalues <- p.adjust(between$pval[,counter],method=method)
            if(min.group<2)
              {
                new.between <- new.qvalues <- new.direction <- rep(NA_real_,n)
                new.between <- between$pval[,counter]
                new.qvalues <- qvalues
                new.direction <- between$direction[,counter]
                output.between[[counter]] <-cbind.data.frame(rownames(rna),new.between,new.qvalues,new.direction)
#                output.between[[counter]] <-cbind.data.frame(rownames(rna),between$pval[,counter],qvalues,between$direction[,counter])
                colnames(output.between[[counter]]) <- c("Gene","P-value","FDR","Direction")
              }
            else
              {
                which.ij <- which(group==ugroup[i]|group==ugroup[j])
                dge <- suppressMessages(DGEList(counts=rna[,which.ij],group=group[which.ij]))
                dsp <- estimateCommonDisp(dge)
                dsp <- estimateTagwiseDisp(dsp)
                res <- topTags(exactTest(dsp,pair=c(ugroup[j],ugroup[i])),n=nrow(rna))$table
                match.genes <- match(rownames(rna),rownames(res))
                mrna.qval <- res$FDR[match.genes]
                mrna.logfc <- res$logFC[match.genes]
                change.type <- rep("both",nrow(rna))
                change.type[mrna.qval>0.05 & abs(mrna.logfc)<1.5] <- "translational_only"
                new.between <- new.qvalues <- new.direction <- rep(NA_real_,n)
                new.between <- between$pval[,counter]
                new.qvalues <- qvalues
                new.direction <- between$direction[,counter]
                output.between[[counter]] <-cbind.data.frame(rownames(rna),mrna.logfc,mrna.qval,change.type,new.between,new.qvalues,new.direction)
#                output.between[[counter]] <-cbind.data.frame(rownames(rna),mrna.logfc,mrna.qval,change.type,between$pval[,counter],qvalues,between$direction[,counter])
                colnames(output.between[[counter]]) <- c("Gene","mRNA_logFC","mRNA_FDR","Change_type","P-value","FDR","Direction")
              }
            counter <- counter+1
          }
      }
    output.between
  }
                          
