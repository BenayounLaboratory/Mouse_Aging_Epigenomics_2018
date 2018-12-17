library('DiffBind')
library("rtracklayer")
# diffbind has changed the way chromosomes are dealt with
# it renames them in teh new version (2.8.0), which creates issues

# my.tissue <- "Heart"
# my.mark <- "H3K4me3"
# my.csv.samples <- "Heart_FIXSEQ_exp_K4me3_Repeats.csv"


run_diffbind <- function (my.tissue, my.mark, my.csv.samples) {
  
  my.prefix <- paste(Sys.Date(),my.tissue,my.mark,"aging_Repeats_DiffBind",sep = "_")
  
  my.dba.obj <- dba(sampleSheet=my.csv.samples,skipLines=1,attributes=c(DBA_ID,DBA_CONDITION))
  my.dba.obj <- dba.count(my.dba.obj)
  
  my.heatmap.out <- paste(my.prefix,"diffbind_heatmap.pdf", sep = "_")
  pdf(my.heatmap.out)
  plot(my.dba.obj, colScheme="Reds")
  dev.off()
  
  my.dba.obj <- dba.contrast(my.dba.obj, categories=DBA_CONDITION, minMembers=2)
  my.dba.obj <- dba.analyze(my.dba.obj,method=DBA_DESEQ2)
  my.dba.obj
  
  my.pca.out <- paste(my.prefix,"PCA_DESeq2.pdf", sep = "_")
  pdf(my.pca.out)
  dba.plotPCA(my.dba.obj,DBA_CONDITION,method=DBA_DESEQ2,cex.pt=0.5)
  dev.off()
  
  # attributes have chnged. Allvectors is not used anymore, access matrix with "binding" instead
  
  # ALSO: because of new implementation of sorting (https://support.bioconductor.org/p/63295/), need to map chromosome names
  # I double checked, this was not the case with the version used for most of the paper (!!!!!)
  # There, all chromosome names are present
  
  my.chromosomes = my.dba.obj$chrmap[my.dba.obj$binding[,1]]
  # my.chromosomes CHR
  # 1               chr1   1
  # 1458           chr10   2
  # 2840           chr11   3
  # 4869           chr12   4
  # 5751           chr13   5
  # 6773           chr14   6
  # 7620           chr15   7
  # 8615           chr16   8
  # 9390           chr17   9
  # 10567          chr18  10
  # 11313          chr19  11
  # 12157           chr2  12
  # 14107           chr3  13
  # 15353           chr4  14
  # 17013           chr5  15
  # 18762           chr6  16
  # 20027           chr7  17
  # 21792           chr8  18
  # 23156           chr9  19
  # 24507           chrX  20
  # 25188           chrY  21
  
  # output result table
  my.res.table <- data.frame(cbind(my.chromosomes,my.dba.obj$binding[,-1])) # remove internal mapping
  write.table(my.res.table,file=paste(my.prefix,"normalized_counts.txt", sep = "_"),quote=F,row.names=F,sep="\t")
  
  # output bed for annotations
  my.bed.table <- my.res.table[,1:3]
  my.bed.table$name <- paste(my.res.table$my.chromosomes,my.res.table$START,my.res.table$END,sep = "-")
  write.table(my.bed.table,file=paste(my.prefix,"peaks.bed", sep = "_"),quote=F,row.names=F,col.names=F,sep="\t")
  
}

