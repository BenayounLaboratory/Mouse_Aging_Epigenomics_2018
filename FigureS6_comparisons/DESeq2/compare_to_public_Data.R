setwd('/Volumes/BB_Backup_3/BD_aging_project/2018-09_revision_analyses/Other_aging_transcriptomes/From_Raw/DEseq2')
options(stringsAsFactors = FALSE)

# 2018-09-24
# compare our results to public datasets suggested to reviewer

load("2018-09-24_Bochkis_Liver_RNAseq.RData")
load("2018-09-24_Boisvert_Cereb_astrocytes_RNAseq.RData")
load("2018-09-24_White_Liver_RNAseq.RData")

my.bochkis.RNAseq.process[[1]] <- data.frame(my.bochkis.RNAseq.process[[1]])
my.boisvert.RNAseq.process[[1]] <- data.frame(my.boisvert.RNAseq.process[[1]])
my.white.RNAseq.process[[1]] <- data.frame(my.white.RNAseq.process[[1]])


my.bochkis.RNAseq.process[[1]]$GeneName <- rownames(my.bochkis.RNAseq.process[[1]])
my.boisvert.RNAseq.process[[1]]$GeneName <- rownames(my.boisvert.RNAseq.process[[1]])
my.white.RNAseq.process[[1]]$GeneName <- rownames(my.white.RNAseq.process[[1]])


# Load corresponding RData RNA with age results
load("/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Separate/RNA_seq_result_Liver_2015-11-19.RData")
load("/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Separate/RNA_seq_result_cereb_2015-11-19.RData")

my.liver.RNAseq.process[[1]] <- data.frame(my.liver.RNAseq.process[[1]])
my.cereb.RNAseq.process[[1]] <- data.frame(my.cereb.RNAseq.process[[1]])

my.liver.RNAseq.process[[1]]$GeneName <- rownames(my.liver.RNAseq.process[[1]])
my.cereb.RNAseq.process[[1]]$GeneName <- rownames(my.cereb.RNAseq.process[[1]])


### merge datasets
my.liver.merged.bochkis <- merge(my.liver.RNAseq.process[[1]],
                                 my.bochkis.RNAseq.process[[1]],
                                 by = "GeneName",
                                 suffixes = c(".benayoun",".bochkis") )

my.liver.merged.white <- merge(my.liver.RNAseq.process[[1]],
                               my.white.RNAseq.process[[1]],
                               by = "GeneName",
                               suffixes = c(".benayoun",".white") )

my.cereb.merged.boisvert <- merge(my.liver.RNAseq.process[[1]],
                                  my.boisvert.RNAseq.process[[1]],
                                  by = "GeneName",
                                  suffixes = c(".benayoun",".boisvert") )


### Bochkis
my.cor <- cor.test(my.liver.merged.bochkis$log2FoldChange.benayoun,my.liver.merged.bochkis$log2FoldChange.bochkis, method = "spearman")
my.cor$p.value # 8.643218e-177

my.cor$estimate # 0.2660927

# p-value < 2.2e-16
# rho 0.2660927 

pdf("Bochkis_liver_correlation.pdf")
smoothScatter(my.liver.merged.bochkis$log2FoldChange.benayoun,
              my.liver.merged.bochkis$log2FoldChange.bochkis,
              ylim = c(-0.2,0.2),
              xlim = c(-0.06,0.06))
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-0.04,0.17,paste0("Rho = ",signif(my.cor$estimate, 4) ))
text(-0.04,0.15,paste0("p = ",signif(my.cor$p.value, 2) ))
dev.off()


### White
my.cor <- cor.test(my.liver.merged.white$log2FoldChange.benayoun,my.liver.merged.white$log2FoldChange.white, method = "spearman")
my.cor$p.value # 0
my.cor$estimate # 0.4082354


pdf("White_liver_correlation.pdf")
smoothScatter(my.liver.merged.white$log2FoldChange.benayoun,
              my.liver.merged.white$log2FoldChange.white,
              ylim = c(-0.2,0.2),
              xlim = c(-0.1,0.1))
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-0.04,0.17,paste0("Rho = ",signif(my.cor$estimate, 4) ))
text(-0.04,0.15,paste0("p = ",signif(my.cor$p.value, 2) ))
dev.off()


### Boisvert
my.cor <- cor.test(my.cereb.merged.boisvert$log2FoldChange.benayoun,my.cereb.merged.boisvert$log2FoldChange.boisvert, method = "spearman")
my.cor$p.value # 1.973564e-34
my.cor$estimate # 0.1154348

pdf("Boisvert_liver_correlation.pdf")
smoothScatter(my.cereb.merged.boisvert$log2FoldChange.benayoun,
              my.cereb.merged.boisvert$log2FoldChange.boisvert,
              ylim = c(-0.2,0.2),
              xlim = c(-0.1,0.1))
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-0.04,0.17,paste0("Rho = ",signif(my.cor$estimate, 4) ))
text(-0.04,0.15,paste0("p = ",signif(my.cor$p.value, 2) ))
dev.off()

