options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
library('limma')
library('affy')

#######################################################################################################################################
## This function takes in the RNAseq matrix and processes it through DEseq modeling with respect to age in months
# INPUT: my.matrix: count matrix from subreads, columns selected, rownames implemented
#        reps.3, reps.12, reps.29: replicates from each age, default is 3
#        my.tissue: name of the tissue

# OUTPUT: list with [[1]] DEseq result object and [[2]] normalized log2 count matrix

process_aging_rnaseq <- function(my.tissue, my.matrix, reps.3=3, reps.12=3, reps.29=3) {
  
  ncols <- dim(my.matrix)[2]
  
  # get output file prefix
  my.outprefix <- paste(Sys.Date(),my.tissue,"DESeq2_LINEAR_model_with_age_with_TEs",sep="_")
  
  # get the genes with no reads out
  my.exclude <- which(apply(my.matrix[,3:ncols], 1, sum) <= 1) # see deseq2 vignetter

  
  my.filtered.matrix <- my.matrix[-my.exclude,3:ncols]
  rownames(my.filtered.matrix) <- my.matrix[-my.exclude,1]
  my.filtered.matrix[is.na(my.filtered.matrix)] <- 0 # need to deal with NA
  
  my.genes <- rownames(my.filtered.matrix)[c(grep('NM',rownames(my.filtered.matrix)), grep('NR',rownames(my.filtered.matrix)))]
  my.repeats <- rownames(my.filtered.matrix)[-c(grep('NM',rownames(my.filtered.matrix)), grep('NR',rownames(my.filtered.matrix)))]
  
  # age vector  
  age <- as.numeric(c(rep(3,reps.3),rep(12,reps.12),rep(29,reps.29) )) # age in months
  
  # design matrix
  dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), age = age )
  
  # get matrix using age as a modeling covariate
  dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                                colData = dataDesign,
                                design = ~ age)
  
  # run DESeq normalizations and export results
  dds.deseq <- DESeq(dds,fitType = "local") # parametric was being skewed!!!
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  res <- results(dds.deseq, name= "age") # added the name of the tested variable: doesn't seem to be taken correctly by default for FC
  
  # normalized expression value
  tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
  colnames(tissue.cts) <- c(paste("3m",1:reps.3,sep=""),paste("12m",1:reps.12,sep=""),paste("29m",1:reps.29,sep=""))

  ### get the heatmap of aging changes at FDR5
  ## exclude NA
  res <- res[!is.na(res$padj),]
  
  genes.aging <- rownames(res)[res$padj < 0.05]
  my.num.aging <- length(intersect(genes.aging,my.repeats) )
  
  if (my.num.aging > 0) {
    # heatmap drawing - only if there is at least one TEs
    my.heatmap.out <- paste(my.outprefix,"_Heatmap_significant_Repeats.pdf")

    pdf(my.heatmap.out, onefile = F)
    my.heatmap.title <- paste(my.tissue," aging significant (FDR<5%), ",my.num.aging, " TEs",sep="")
    pheatmap(tissue.cts[intersect(genes.aging,my.repeats),],
             cluster_cols = F,
             cluster_rows = T,
             colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
             show_rownames = T, scale="row",
             main = my.heatmap.title, cellwidth = 30)
    dev.off()
  }
  
  
  # output result tables to files
  my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix.txt")
  my.out.fdr5 <- paste(my.outprefix,"_FDR5_repeats_statistics.txt")
  my.out.rdata <- paste(my.outprefix,"_statistics.RData")
  
  write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
  write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
  write.table(res[intersect(genes.aging,my.repeats),], file = my.out.fdr5, sep = "\t" , row.names = T, quote=F)
  
  return(list(res,tissue.cts))
}

#######################################################################################################################################