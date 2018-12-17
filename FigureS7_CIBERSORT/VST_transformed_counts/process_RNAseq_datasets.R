setwd('/Volumes/MyBook_3/BD_aging_project/RNAseq/All_tissues_analysis/CIBERSORT/VST_transformed_coutns/')
source('RNAseq_analysis_functions.R')

library('DESeq2')

# 2016-10-19
# run to extract VST transformed counts

# 2017-04-05
# add mouse alzheimer's data and 40Hz treated hippocampus
# run to extract VST transformed counts

####################################    Liver    #################################### 
# read in subread count matrix
my.liver1 <- read.table('/Volumes/MyBook_3/BD_aging_project/RNAseq/Liver/STAR/Aging_Liver_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.liver <- my.liver1[,c(1,6:15)]
rownames(my.liver) <- my.liver[,1]

# process RNAseq data and save RData object
my.liver.RNAseq.process <- process_rnaseq_to_VST("Liver", my.liver)
save(my.liver.RNAseq.process, file="RNA_seq_result_Liver_2016-10-18.RData")
##################################################################################### 

####################################   Heart   #################################### 
# read in subread count matrix
my.heart1 <- read.table('/Volumes/MyBook_3/BD_aging_project/RNAseq/Heart/STAR/Aging_Heart_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.heart <- my.heart1[,c(1,6:15)]
rownames(my.heart) <- my.heart[,1]

# process RNAseq data and save RData object
my.heart.RNAseq.process <- process_rnaseq_to_VST("Heart", my.heart)
save(my.heart.RNAseq.process, file="RNA_seq_result_Heart_2016-10-18.RData")
###################################################################################

#################################### Cerebellum #################################### 
# read in subread count matrix
# there were 2 nextseq runs based on poor clustering on flow cell
# will sum up count matrices
my.cereb1 <- read.table('/Volumes/MyBook_3/BD_aging_project/RNAseq/Cereb/1st_run/STAR/Aging_cerebellum_counts_genes.txt',skip=1,header=T,sep="\t")
my.cereb2 <- read.table('/Volumes/MyBook_3/BD_aging_project/RNAseq/Cereb/2nd_run/STAR/Aging_cerebellum_v2_counts_genes.txt',skip=1,header=T,sep="\t")

my.cereb <- my.cereb1[,c(1,6:15)]
my.cereb[,3:11] <- my.cereb[,3:11] + my.cereb2[,7:15]
rownames(my.cereb) <- my.cereb[,1]

# process RNAseq data and save RData object
my.cereb.RNAseq.process <- process_rnaseq_to_VST("Cerebellum", my.cereb)
save(my.cereb.RNAseq.process, file="RNA_seq_result_cereb_2016-10-18.RData")
####################################################################################

#################################### Olfactory Bulb #################################
# one of the 12mths samples was not analyzed
# read in subread count matrix
my.ob1 <- read.table('/Volumes/MyBook_3/BD_aging_project/RNAseq/OB/STAR/Aging_OlfactoryBulb_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.ob <- my.ob1[,c(1,6:14)]
rownames(my.ob) <- my.ob[,1]

# process RNAseq data and save RData object
my.ob.RNAseq.process <- process_rnaseq_to_VST("OlfactoryBulb", my.ob, reps.3=3, reps.12=2, reps.29=3)
save(my.ob.RNAseq.process, file="RNA_seq_result_OB_2016-10-18.RData")
##################################################################################### 

####################################  NPCs pools  ###################################
# read in subread count matrix
my.npc1 <- read.table('/Volumes/MyBook_3/BD_aging_project/RNAseq/NPC_Pool/STAR/Aging_NPCs_pool_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.npc <- my.npc1[,c(1,6:12)]
rownames(my.npc) <- my.npc[,1]

# process RNAseq data and save RData npcject
my.npc.RNAseq.process <- process_rnaseq_to_VST("NPCs", my.npc, reps.3=2, reps.12=2, reps.29=2)
save(my.npc.RNAseq.process, file="RNA_seq_result_NPCs_2016-10-18.RData")
##################################################################################### 


####################################################################################################################################
# Also normalize publically available data

####################################  Alzheimer's  ###################################
# read in subread count matrix
my.azh1 <- read.table('/Volumes/LaCie/Disease_model/Alzheimer_model/RNAseq/STAR/Alzheimer_model_expression_Hippocampus_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.azh <- my.azh1[,c(1,6:18)]
rownames(my.azh) <- my.azh[,1]

# clean up input
my.azh.proc <- preprocess_matrix(my.azh)

# design matrix
my.age = rep(c(rep(2,3),rep(6,3) ) ,2)
my.genotype = c(rep("WT",6),rep("CKp25",6))
dataDesign = data.frame( row.names = colnames( my.azh.proc ), age = my.age, genotype = my.genotype )

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.azh.proc, colData = dataDesign,design = ~ age + genotype)

# run DESeq and export normalized expression values
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- getVarianceStabilizedData(dds)
# output result tables to files
my.out.ct.mat <- paste(Sys.Date(),"Alzheimers_model_log2_VST_counts_matrix.txt")
write.table(vsd, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
##################################################################################### 

####################################  Hippo 40Hz  ###################################
# read in subread count matrix
my.hippo1 <- read.table('/Volumes/LaCie/Disease_model/Hippocampus_40Hz/STAR/Hippocampus_gamma40Hz_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.hippo <- my.hippo1[,c(1,6:12)]
rownames(my.hippo) <- my.hippo[,1]

# clean up input
my.hippo.proc <- preprocess_matrix(my.hippo)

# design matrix
my.treatment = c(rep("CTL",3),rep("GAMMA",3))
dataDesign = data.frame( row.names = colnames( my.hippo.proc ), treatment = my.treatment )

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.hippo.proc, colData = dataDesign,design = ~ treatment)

# run DESeq and export normalized expression values
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- getVarianceStabilizedData(dds)
# output result tables to files
my.out.ct.mat <- paste(Sys.Date(),"40Hz_gamma_hippocampus_model_log2_VST_counts_matrix.txt")
write.table(vsd, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
##################################################################################### 