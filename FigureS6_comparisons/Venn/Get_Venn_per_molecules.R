setwd('/Volumes/BB_Backup_3/BD_aging_project/2018-09_revision_analyses/Other_aging_transcriptomes/From_Raw/Venn/')
options(stringsAsFactors=F)

library('Vennerable')

# 2018-09-24
# run per tissue Venn diagramms
# with the reanalysis through our pipelime

my.liver.rna.sig.bb       <- read.csv('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Separate/2015-11-19_Liver_DESeq2_LINEAR_model_with_age _FDR5_genes_statistics.txt',sep="\t",header=T)
my.liver.rna.sig.bochkis  <- read.csv('../DEseq2/2018-09-24_Liver_Bochkis_2015_DESeq2_LINEAR_model_with_age _FDR5_genes_statistics.txt',sep="\t",header=T)
my.liver.rna.sig.white    <- read.csv('../DEseq2/2018-09-24_Liver_White_2015_DESeq2_LINEAR_model_with_age _FDR5_genes_statistics.txt',sep="\t",header=T)

my.cereb.rna.sig.bb       <- read.csv('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Separate/2015-11-19_Cerebellum_DESeq2_LINEAR_model_with_age _FDR5_genes_statistics.txt',sep="\t",header=T)
my.cereb.rna.sig.boisvert <- read.csv('../DEseq2/2018-09-24_Cereb_astrocytes_Boisvert_2015_DESeq2_LINEAR_model_with_age _FDR5_genes_statistics.txt',sep="\t",header=T)


my.liver.rna.all.bb       <- read.csv('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Separate/2015-11-19_Liver_DESeq2_LINEAR_model_with_age _all_genes_statistics.txt',sep="\t",header=T)
my.liver.rna.all.bochkis  <- read.csv('../DEseq2/2018-09-24_Liver_Bochkis_2015_DESeq2_LINEAR_model_with_age _all_genes_statistics.txt',sep="\t",header=T)
my.liver.rna.all.white    <- read.csv('../DEseq2/2018-09-24_Liver_White_2015_DESeq2_LINEAR_model_with_age _all_genes_statistics.txt',sep="\t",header=T)
my.liver.bckd <- unique(c(rownames(my.liver.rna.all.bb),rownames(my.liver.rna.all.bochkis),rownames(my.liver.rna.all.white)))

my.cereb.rna.all.bb       <- read.csv('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Separate/2015-11-19_Cerebellum_DESeq2_LINEAR_model_with_age _all_genes_statistics.txt',sep="\t",header=T)
my.cereb.rna.all.boisvert <- read.csv('../DEseq2/2018-09-24_Cereb_astrocytes_Boisvert_2015_DESeq2_LINEAR_model_with_age _all_genes_statistics.txt',sep="\t",header=T)
my.cereb.bckd <- unique(c(rownames(my.cereb.rna.all.boisvert),rownames(my.cereb.rna.all.bb)))

get_contigency_mat <- function (list1, list2, background) {
  my.both <- intersect(list1,list2)
  my.1.not.2 <- setdiff(list1,list2)
  my.2.not.1 <- setdiff(list2,list1)
  my.neither <- setdiff(background,unique(c(list1,list2)))
  
  my.contingency.mat <- matrix(c(length(my.both),
                                 length(my.1.not.2),
                                 length(my.2.not.1),
                                 length(my.neither)),
                               2,2)
  return(my.contingency.mat)
}


#####################################################################################################
######################################      Liver results        ######################################  

my.liver.rna_up.sig.bochkis <- rownames(my.liver.rna.sig.bb)[my.liver.rna.sig.bb$log2FoldChange > 0]
my.liver.rna_up.sig.white   <- rownames(my.liver.rna.sig.bochkis)[my.liver.rna.sig.bochkis$log2FoldChange > 0]
my.liver.rna_up.sig.bb      <- rownames(my.liver.rna.sig.white)[my.liver.rna.sig.white$log2FoldChange > 0]

my.liver.rna_dwn.sig.bochkis <- rownames(my.liver.rna.sig.bb)[my.liver.rna.sig.bb$log2FoldChange < 0]
my.liver.rna_dwn.sig.white   <- rownames(my.liver.rna.sig.bochkis)[my.liver.rna.sig.bochkis$log2FoldChange < 0]
my.liver.rna_dwn.sig.bb      <- rownames(my.liver.rna.sig.white)[my.liver.rna.sig.white$log2FoldChange < 0]


######### 1. Bochkis et al
# UP 
my.RNA.up <- list("Bochkis" = unique(my.liver.rna_up.sig.bochkis),
                  "Benayoun" = unique(my.liver.rna_up.sig.bb))
Venn.RNA.up <- Venn(my.RNA.up)

pdf("2018-09-24_Venn_Aging_Liver_RNA_UP_Bochkis.pdf")
plot(Venn.RNA.up, doWeights=T)
dev.off()

#### 
my.up.boch.mat <- get_contigency_mat(my.liver.rna_up.sig.bochkis, my.liver.rna_up.sig.bb, my.liver.bckd)
my.test <- fisher.test(my.up.boch.mat, alternative = "greater")
my.test$p.value
# p-value 3.268652e-67

# DWN
my.RNA.dwn <- list("Bochkis" = unique(my.liver.rna_dwn.sig.bochkis),
                  "Benayoun" = unique(my.liver.rna_dwn.sig.bb))
Venn.RNA.dwn <- Venn(my.RNA.dwn)

pdf("2018-09-24_Venn_Aging_Liver_RNA_DWN_Bochkis.pdf")
plot(Venn.RNA.dwn, doWeights=T)
dev.off()
#### 
my.dwn.boch.mat <- get_contigency_mat(my.liver.rna_dwn.sig.bochkis, my.liver.rna_dwn.sig.bb, my.liver.bckd)
my.test <- fisher.test(my.dwn.boch.mat, alternative = "greater")
my.test$p.value
# p-value 4.261916e-23


######### 2. White et al
# UP 
my.RNA.up <- list("White" = unique(my.liver.rna_up.sig.white),
                  "Benayoun" = unique(my.liver.rna_up.sig.bb))
Venn.RNA.up <- Venn(my.RNA.up)

pdf("2018-09-24_Venn_Aging_Liver_RNA_UP_White.pdf")
plot(Venn.RNA.up, doWeights=T)
dev.off()

#### 
my.up.white.mat <- get_contigency_mat(my.liver.rna_up.sig.white, my.liver.rna_up.sig.bb, my.liver.bckd)
my.test <- fisher.test(my.up.white.mat, alternative = "greater")
my.test$p.value
# p-value 3.523543e-50

# DWN
my.RNA.dwn <- list("White" = unique(my.liver.rna_dwn.sig.white),
                   "Benayoun" = unique(my.liver.rna_dwn.sig.bb))
Venn.RNA.dwn <- Venn(my.RNA.dwn)

pdf("2018-09-24_Venn_Aging_Liver_RNA_DWN_White.pdf")
plot(Venn.RNA.dwn, doWeights=T)
dev.off()
#### 
my.dwn.white.mat <- get_contigency_mat(my.liver.rna_dwn.sig.white, my.liver.rna_dwn.sig.bb, my.liver.bckd)
my.test <- fisher.test(my.dwn.white.mat, alternative = "greater")
my.test$p.value
# p-value 6.961207e-21


#####################################################################################################
######################################      Cereb results        ######################################  
my.cereb.rna_up.sig.boisvert <- rownames(my.cereb.rna.sig.boisvert)[my.cereb.rna.sig.boisvert$log2FoldChange > 0]
my.cereb.rna_up.sig.bb       <- rownames(my.cereb.rna.sig.bb)[my.cereb.rna.sig.bb$log2FoldChange > 0]

my.cereb.rna_dwn.sig.boisvert <- rownames(my.cereb.rna.sig.boisvert)[my.cereb.rna.sig.boisvert$log2FoldChange < 0]
my.cereb.rna_dwn.sig.bb       <- rownames(my.cereb.rna.sig.bb)[my.cereb.rna.sig.bb$log2FoldChange < 0]


# UP
my.RNA.up <- list("Boisvert" = unique(my.cereb.rna_up.sig.boisvert),
                  "Benayoun" = unique(my.cereb.rna_up.sig.bb))
Venn.RNA.up <- Venn(my.RNA.up)

pdf("2018-09-24_Venn_Aging_Cereb_RNA_UP.pdf")
plot(Venn.RNA.up, doWeights=T)
dev.off()

#### 
my.up.bois.mat <- get_contigency_mat(my.cereb.rna_up.sig.boisvert, my.cereb.rna_up.sig.bb, my.cereb.bckd)
my.test <- fisher.test(my.up.bois.mat, alternative = "greater")
my.test$p.value
# p-value 4.642509e-49


# DWN
my.RNA.dwn <- list("Boisvert" = unique(my.cereb.rna_dwn.sig.boisvert),
                   "Benayoun" = unique(my.cereb.rna_dwn.sig.bb))
Venn.RNA.dwn <- Venn(my.RNA.dwn)

pdf("2018-09-24_Venn_Aging_Cereb_RNA_DWN.pdf")
plot(Venn.RNA.dwn, doWeights=T)
dev.off()

my.dwn.bois.mat <- get_contigency_mat(my.cereb.rna_dwn.sig.boisvert, my.cereb.rna_dwn.sig.bb, my.cereb.bckd)
my.test <- fisher.test(my.dwn.bois.mat, alternative = "greater")
my.test$p.value
# p-value = 0.001284

