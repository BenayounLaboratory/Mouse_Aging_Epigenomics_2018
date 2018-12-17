setwd('/Volumes/BB_Backup_3/BD_aging_project/2018-09_revision_analyses/Other_aging_transcriptomes/From_Raw/DEseq2')
options(stringsAsFactors = FALSE)

source('RNAseq_analysis_functions.R')
# 2018-09-24
# analyze other published RNA-seq for revision

####################################    White    #################################### 
# read in subread count matrix
my.white1 <- read.table('../STAR/White_aging_Liver_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.white <- my.white1[,c(1,6:12)]
rownames(my.white) <- my.white[,1]

my.white.RNAseq.process <- process_aging_rnaseq("Liver_White_2015", my.white, c(rep(4,3),rep(28,3)))
save(my.white.RNAseq.process, file = "2018-09-24_White_Liver_RNAseq.RData")
##################################################################################### 

####################################    Bochkis    #################################### 
# read in subread count matrix
my.bochkis1 <- read.table('../STAR/Bochkis_aging_Liver_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.bochkis <- my.bochkis1[,c(1,6:13)]
rownames(my.bochkis) <- my.bochkis[,1]

my.bochkis.RNAseq.process <- process_aging_rnaseq("Liver_Bochkis_2015", my.bochkis, c(rep(3,4),rep(21,3)))
save(my.bochkis.RNAseq.process, file = "2018-09-24_Bochkis_Liver_RNAseq.RData")
#####################################################################################  

####################################    Boisvert    #################################### 
# read in subread count matrix
my.boisvert1 <- read.table('../STAR/Boisvert_aging_cerebellum_astrocytes_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.boisvert <- my.boisvert1[,c(1,6:14)]
rownames(my.boisvert) <- my.boisvert[,1]

my.boisvert.RNAseq.process <- process_aging_rnaseq("Cereb_astrocytes_Boisvert_2015", my.boisvert, c(rep(4,4),rep(24,4)))
save(my.boisvert.RNAseq.process, file = "2018-09-24_Boisvert_Cereb_astrocytes_RNAseq.RData")
#####################################################################################  