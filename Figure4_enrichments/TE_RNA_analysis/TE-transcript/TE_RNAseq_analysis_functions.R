options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
#library('bitops')
library('limma')
library('affy')

#######################################################################################################################################
process_TE_rnaseq <- function(my.tissue, my.matrix, reps.3=3, reps.12=3, reps.29=3, my.TE.names) {
  
  ncols <- dim(my.matrix)[2]
  
  # get output file prefix
  my.outprefix <- paste(Sys.Date(),my.tissue,"DESeq2_TE_LINEAR_model_with_age",sep="_")
  
  # get the genes with no reads out
  my.null <- which(apply(my.matrix[,2:ncols], 1, sum) <= 1) # see deseq2 vignetter
  # Now pull out the spike in genes
  spikes.idx <- grep("ERCC-", rownames(my.matrix))
  my.exclude <- union(my.null,spikes.idx)
  
  my.filtered.matrix <- my.matrix[-my.exclude,2:ncols]
  rownames(my.filtered.matrix) <- my.matrix[-my.exclude,1]
  
  age <- as.numeric(c(rep(3,reps.3),rep(12,reps.12),rep(29,reps.29) )) # age in months
  
  # design matrix
  dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), age = age )
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                                colData = dataDesign,
                                design = ~ age)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  res <- results(dds.deseq, name= "age")
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  
  # normalized expression value
  tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
  colnames(tissue.cts) <- c(paste("3m",1:reps.3,sep=""),paste("12m",1:reps.12,sep=""),paste("29m",1:reps.29,sep=""))
  
  # get TE indeces
  TE.idx <- rownames(tissue.cts) %in% my.TE.names$V1
  
  ##### do MDS analysis
  mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  mds.result.te <- cmdscale(1-cor(tissue.cts[TE.idx,],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x.te <- mds.result.te[, 1]
  y.te <- mds.result.te[, 2]
  
  my.colors <- c(rep("coral",reps.3), rep("blueviolet", reps.12),rep("dodgerblue",reps.29))
  
  my.mds.out <- paste(my.outprefix,"_MDS_plots.pdf")
  
  pdf(my.mds.out, width = 10, height = 5)
  plot(x.te, y.te, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling (TE Only)",cex=2)
  points(x.te, y.te, pch=16,col=my.colors,cex=2)
  legend("topleft",c("3m","12m","29m"),col=c("coral","blueviolet","dodgerblue"),pch=16,bty='n',pt.cex=2)
  par(mfrow=c(1,1))
  dev.off()
  
    
  # expression range
  my.exp.out <- paste(my.outprefix,"_Normalized_counts_boxplot.pdf")
  
  pdf(my.exp.out)
  boxplot(tissue.cts,col=c(rep("coral",reps.3),rep("blueviolet",reps.12),rep("dodgerblue",reps.29)),
          cex=0.5,ylab="Log2 DESeq2 Normalized counts", main = my.tissue) 
  boxplot(tissue.cts[TE.idx,],col=c(rep("coral",reps.3),rep("blueviolet",reps.12),rep("dodgerblue",reps.29)),
          cex=0.5,ylab="Log2 DESeq2 Normalized counts", main = paste(my.tissue, "TE_only") )  
  dev.off()
  
  ### get the heatmap of aging changes at FDR5
  ## exclude NA
  res <- res[!is.na(res$padj),]
    
  # exlude non TE
  TE.idx.2 <- rownames(res) %in% my.TE.names$V1
  res <- res[TE.idx.2,]
  
  genes.aging <- rownames(res)[res$padj < 0.05]
  my.num.aging <- length(genes.aging)
  
  if (my.num.aging > 0) {
    # heatmap drawing - only if there is at least one gene
    my.heatmap.out <- paste(my.outprefix,"_Heatmap_significant_TEs.pdf")
    
    pdf(my.heatmap.out, onefile=F)
    my.heatmap.title <- paste(my.tissue," aging singificant (FDR<5%), ",my.num.aging, " genes",sep="")
    pheatmap(tissue.cts[genes.aging,],
             cluster_cols = F,
             cluster_rows = T,
             colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
             show_rownames = T, scale="row",
             main = my.heatmap.title, cellwidth = 30)
    dev.off()
  }

    
  
  # output result tables to files
  my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_ALL.txt")
  my.out.stats <- paste(my.outprefix,"_all_TE_statistics.txt")
  my.out.fdr5 <- paste(my.outprefix,"_FDR5_TE_statistics.txt")
  my.out.rdata <- paste(my.outprefix,"_statistics.RData")
  
  write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
  write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
  write.table(res[genes.aging,], file = my.out.fdr5, sep = "\t" , row.names = T, quote=F)
  
  #save(res,my.out.rdata)
 
  return(list(res,tissue.cts))
}

#######################################################################################################################################
