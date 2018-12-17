options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')


#######################################################################################################################################
## This function takes in the RNAseq matrix and processes it through DEseq modeling with respect to age in months
# INPUT: my.matrix: count matrix from subreads, columns selected, rownames implemented
#        reps.y, reps.12, reps.o: replicates from each age, default is 3
#        my.tissue: name of the tissue

# OUTPUT: list with [[1]] DEseq result object and [[2]] normalized log2 count matrix

process_aging_rnaseq <- function(my.tissue, my.matrix, my.ages) {
  
  ncols <- dim(my.matrix)[2]
  
  # get output file prefix
  my.outprefix <- paste(Sys.Date(),my.tissue,"DESeq2_LINEAR_model_with_age",sep="_")
  
  # get the genes with no reads out
  my.null <- which(apply(my.matrix[,3:ncols], 1, sum) <= 1) # see deseq2 vignetter
  # Now pull out the spike in genes
  spikes.idx <- grep("ERCC-", rownames(my.matrix))
  my.exclude <- union(my.null,spikes.idx)
  
  my.filtered.matrix <- my.matrix[-my.exclude,3:ncols]
  rownames(my.filtered.matrix) <- my.matrix[-my.exclude,1]
  
  age <- as.numeric(my.ages) # age in months
  
  # design matrix
  dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), age = age )
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                                colData = dataDesign,
                                design = ~ age)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds)
  res <- results(dds.deseq, name= "age") # added the name of the tested variable: doesn't seem to be taken correctly by default for FC
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  
  # normalized expression value
  tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
  colnames(tissue.cts) <- paste(my.ages,"m",1:length(my.ages),sep="")
  
  # do MDS analysis
  mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
  x <- mds.result[, 1]
  y <- mds.result[, 2]
  
  my.pos <- which(my.ages == min(my.ages))
  
  my.colors <- c(rep("coral",max(my.pos)),rep("dodgerblue",length(my.ages) - max(my.pos) + 1))
  
  my.mds.out <- paste(my.outprefix,"_MDS_plot.pdf")
  
  pdf(my.mds.out)
  plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2", main="Multi-dimensional Scaling",cex=2)
  points(x, y, pch=16,col=my.colors,cex=2)
  legend("topleft",c(paste(min(my.ages),"m",sep=""),paste(max(my.ages),"m")),col=c("coral","dodgerblue"),pch=16,bty='n',pt.cex=2)
  dev.off()
  
  # expression range
  my.exp.out <- paste(my.outprefix,"_Normalized_counts_boxplot.pdf")
  
  pdf(my.exp.out)
  boxplot(tissue.cts,col=my.colors,
          cex=0.5,ylab="Log2 DESeq2 Normalized counts", main = my.tissue)  
  dev.off()
  
  ### get the heatmap of aging changes at FDR5
  ## exclude NA
  res <- res[!is.na(res$padj),]
    
  genes.aging <- rownames(res)[res$padj < 0.05]
  my.num.aging <- length(genes.aging)
  
  if (my.num.aging > 2) {
    # heatmap drawing - only if there is at least one gene
    my.heatmap.out <- paste(my.outprefix,"_Heatmap_significant_genes.pdf")
    
    pdf(my.heatmap.out)
    my.heatmap.title <- paste(my.tissue," aging singificant (FDR<5%), ",my.num.aging, " genes",sep="")
    pheatmap(tissue.cts[genes.aging,],
             cluster_cols = F,
             cluster_rows = T,
             colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
             show_rownames = F, scale="row",
             main = my.heatmap.title, cellwidth = 30)
    dev.off()
  }

  
  # do clustering
  my.pv <- pvclust(tissue.cts,nboot=50)
  my.heatmap.out <- paste(my.outprefix,"_PVCLUST_result.pdf")
  
  pdf(my.heatmap.out)
  plot(my.pv)
  dev.off()
  
  
  # output result tables to files
  my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix.txt")
  my.out.stats <- paste(my.outprefix,"_all_genes_statistics.txt")
  my.out.fdr5 <- paste(my.outprefix,"_FDR5_genes_statistics.txt")
  my.out.rdata <- paste(my.outprefix,"_statistics.RData")
  
  write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
  write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
  write.table(res[genes.aging,], file = my.out.fdr5, sep = "\t" , row.names = T, quote=F)
  
  return(list(res,tissue.cts))
}

#######################################################################################################################################