setwd('/Volumes/BB_USC_1/TE_aging/HOMER/DESeq2')
source('TE_DE_functions_forAnne.R.R')


# 2017-12-15
# Use homer to get counts and perfrom DE on TE elements

# 2018-03-05
# continue

# 2018-03-13
# run for Anne a

####################################    Liver    #################################### 
# read in subread count matrix
my.liver <- read.csv('Count_table_HOMER/2017-12-12_Liver_RNA_genes_and_repeats_matrix.txt',header=T,sep="\t",stringsAsFactors=F)
rownames(my.liver) <- my.liver[,1]

# process RNAseq data and save RData object
my.liver.RNAseq.process <- process_aging_rnaseq("Liver", my.liver)
##################################################################################### 

####################################   Heart   #################################### 
# read in subread count matrix
my.heart <- read.csv('Count_table_HOMER/2017-12-12_heart_RNA_genes_and_repeats_matrix.txt',header=T,sep="\t",stringsAsFactors=F)
rownames(my.heart) <- my.heart[,1]

# process RNAseq data and save RData object
my.heart.RNAseq.process <- process_aging_rnaseq("Heart", my.heart)
###################################################################################

#################################### Cerebellum ####################################
# read in subread count matrix
my.cereb <- read.csv('Count_table_HOMER/2017-12-12_Cereb_RNA_genes_and_repeats_matrix.txt',header=T,sep="\t",stringsAsFactors=F)
rownames(my.cereb) <- my.cereb[,1]

# process RNAseq data and save RData object
my.cereb.RNAseq.process <- process_aging_rnaseq("Cerebellum", my.cereb)
####################################################################################

#################################### Olfactory Bulb #################################
# one of the 12mths samples was not analyzed
# read in subread count matrix
my.ob <- read.csv('Count_table_HOMER/2017-12-12_OB_RNA_genes_and_repeats_matrix.txt',header=T,sep="\t",stringsAsFactors=F)
rownames(my.ob) <- my.ob[,1]

# process RNAseq data and save RData object
my.ob.RNAseq.process <- process_aging_rnaseq("OlfactoryBulb", my.ob, reps.3=3, reps.12=2, reps.29=3)
#####################################################################################

####################################  NPCs pools  ###################################
# read in subread count matrix
my.npc <- read.csv('Count_table_HOMER/2017-12-12_NPC_RNA_genes_and_repeats_matrix.txt',header=T,sep="\t",stringsAsFactors=F)
rownames(my.npc) <- my.npc[,1]

# process RNAseq data and save RData object
my.npc.RNAseq.process <- process_aging_rnaseq("NPCs", my.npc, reps.3=2, reps.12=2, reps.29=2)
#####################################################################################

