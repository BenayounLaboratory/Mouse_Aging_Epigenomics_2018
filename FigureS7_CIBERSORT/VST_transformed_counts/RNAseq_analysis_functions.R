options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
#library('bitops')


#######################################################################################################################################
## This function takes in the RNAseq matrix and processes it through DEseq modeling with respect to age in months
# INPUT: my.matrix: count matrix from subreads, columns selected, rownames implemented
#        reps.3, reps.12, reps.29: replicates from each age, default is 3
#        my.tissue: name of the tissue

# OUTPUT: list with [[1]] DEseq result object and [[2]] normalized log2 count matrix

process_rnaseq_to_VST <- function(my.tissue, my.matrix, reps.3=3, reps.12=3, reps.29=3) {
  
  ncols <- dim(my.matrix)[2]
  
  # get output file prefix
  my.outprefix <- paste(Sys.Date(),my.tissue,"DESeq2_age_RNAseq_normalizations",sep="_")
  
  # get the genes with no reads out
  my.null <- which(apply(my.matrix[,3:ncols], 1, sum) <= 1) # see deseq2 vignetter
  # Now pull out the spike in genes
  spikes.idx <- grep("ERCC-", rownames(my.matrix))
  my.exclude <- union(my.null,spikes.idx)
  
  my.filtered.matrix <- my.matrix[-my.exclude,3:ncols]
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
  
  
  # plot dispersion
  my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")
  
  pdf(my.disp.out)
  plotDispEsts(dds.deseq)
  dev.off()
  
  # normalized expression value
  vsd <- getVarianceStabilizedData(dds.deseq)
  colnames(vsd) <- c(paste("3m",1:reps.3,sep=""),paste("12m",1:reps.12,sep=""),paste("29m",1:reps.29,sep=""))
  
  # expression range
  my.exp.out <- paste(my.outprefix,"_Normalized_counts_boxplot.pdf")
  
  pdf(my.exp.out)
  boxplot(vsd,col=c(rep("coral",reps.3),rep("blueviolet",reps.12),rep("dodgerblue",reps.29)),
          cex=0.5,ylab="Log2 DESeq2 VST-Normalized counts", main = my.tissue)  
  dev.off()
  
  # output result tables to files
  my.out.ct.mat <- paste(my.outprefix,"log2_VST_counts_matrix.txt")

  
  write.table(vsd, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
  
  return(vsd)
}

#######################################################################################################################################


preprocess_matrix <- function(my.matrix) {
  
  ncols <- dim(my.matrix)[2]
  
  # get the genes with no reads out
  my.null <- which(apply(my.matrix[,3:ncols], 1, sum) <= 1) # see deseq2 vignetter
  # Now pull out the spike in genes
  spikes.idx <- grep("ERCC-", rownames(my.matrix))
  my.exclude <- union(my.null,spikes.idx)
  
  my.filtered.matrix <- my.matrix[-my.exclude,3:ncols]
  rownames(my.filtered.matrix) <- my.matrix[-my.exclude,1]
    
  return(my.filtered.matrix)
}

#######################################################################################################################################