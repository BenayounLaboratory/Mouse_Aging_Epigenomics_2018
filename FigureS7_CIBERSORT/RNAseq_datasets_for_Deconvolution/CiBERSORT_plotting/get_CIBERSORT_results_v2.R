setwd('/Volumes/BB_Backup_3//RNAseq_datasets_for_Deconvolution/Process_for_deconvolution/CiBERSORT_plotting')

library('pheatmap')



ciber.training <- read.csv('../2017-01-18/Training_CIBERSORT.Output_Job43.csv',row.names=1)

ciber.liver <- read.csv('../2017-01-18/Liver_CIBERSORT.Output_Job47.csv',row.names=1)
ciber.heart <- read.csv('../2017-01-18/Heart_CIBERSORT.Output_Job48.csv',row.names=1)
ciber.NPCs <- read.csv('../2017-01-18/NPCs_CIBERSORT.Output_Job46.csv',row.names=1)
ciber.OB <- read.csv('../2017-01-18/OB_CIBERSORT.Output_Job45.csv',row.names=1)
ciber.Cereb <- read.csv('../2017-01-18/Cerebellum_Cerebellum_CIBERSORT.Output_Job44.csv',row.names=1)
ciber.mix <- read.csv('../InSilicoMixtures/Mixtures_CIBERSORT.Output_Job49.csv',row.names=1)

colnames(ciber.liver)
# [1] "Input.Sample"        "Adipocytes"          "aNSCs"               "Astrocytes"          "Endothelial_cells"   "Erythrocytes"        "Granulocytes"       
# [8] "Macrophages"         "Microglia"           "Monocytes"           "Neurons"             "NK_cells"            "NPCs"                "Oligodendrocytes"   
# [15] "qNSCs"               "T_cells"             "Cardiomyocytes"      "Hepatocytes"         "B_cells"             "Cardiac_Fibroblasts" "Dermal_Fibroblasts" 
# [22] "Ependymal"           "P.value"             "Pearson.Correlation" "RMSE"

my.cell.cols.order <- c("Cardiomyocytes",
                        "Hepatocytes",
                        "aNSCs","NPCs","qNSCs",
                        "Astrocytes","Oligodendrocytes","Neurons","Ependymal",
                        "Cardiac_Fibroblasts", "Dermal_Fibroblasts" ,
                        "Monocytes","Macrophages","Microglia","Granulocytes",
                        "B_cells","T_cells","NK_cells",
                        "Endothelial_cells",
                        "Adipocytes",
                        "Erythrocytes"
                        )


################################################################################################################
##### heatmaps WITHOUT numbers
pheatmap(ciber.training[my.cell.cols.order,my.cell.cols.order], cluster_rows = F,  cluster_cols = F,
         display_numbers=F, cellwidth = 20 , cellheight = 20,
         filename = paste(Sys.Date(), "CIBERSORT_results_Training_cell_types.pdf", sep="_") )

pheatmap(ciber.liver[,my.cell.cols.order],cluster_rows = F, cluster_cols = F,
         display_numbers=F, cellwidth = 20 , cellheight = 20,
         filename = paste(Sys.Date(), "CIBERSORT_results_Liver_aging.pdf", sep="_") )

pheatmap(ciber.heart[,my.cell.cols.order],cluster_rows = F, cluster_cols = F, 
         display_numbers=F, cellwidth = 20 , cellheight = 20,
         filename = paste(Sys.Date(), "CIBERSORT_results_Heart_aging.pdf", sep="_") )

pheatmap(ciber.NPCs[,my.cell.cols.order],cluster_rows = F, cluster_cols = F,
         display_numbers=F, cellwidth = 20 , cellheight = 20,
         filename = paste(Sys.Date(), "CIBERSORT_results_NPCs_aging.pdf", sep="_") )

pheatmap(ciber.OB[,my.cell.cols.order],cluster_rows = F, cluster_cols = F, 
         display_numbers=F, cellwidth = 20 , cellheight = 20,
         filename = paste(Sys.Date(), "CIBERSORT_results_OB_aging.pdf", sep="_") )

pheatmap(ciber.Cereb[,my.cell.cols.order],cluster_rows = F, cluster_cols = F, 
         display_numbers=F, cellwidth = 20 , cellheight = 20,
         filename = paste(Sys.Date(), "CIBERSORT_results_Cerebellum_aging.pdf", sep="_") )

pheatmap(ciber.mix[,my.cell.cols.order],cluster_rows = F, cluster_cols = F, 
         display_numbers=F, cellwidth = 20 , cellheight = 20,
         filename = paste(Sys.Date(), "CIBERSORT_results_InSilicoMix.pdf", sep="_") )


########################################################################################################################
# beeplotS
library('beeswarm')

ciber.liver$age <- factor(c(rep("3m",3),rep("12m",3),rep("29m",3)))
ciber.heart$age <- factor(c(rep("3m",3),rep("12m",3),rep("29m",3)))
ciber.NPCs$age <- factor(c(rep("3m",2),rep("12m",2),rep("29m",2)))
ciber.OB$age <- factor(c(rep("3m",3),rep("12m",2),rep("29m",3)))
ciber.Cereb$age <- factor(c(rep("3m",3),rep("12m",3),rep("29m",3)))

my.cols <- c("coral","blueviolet","dodgerblue")

pdf("CIBERSORT_Inflammatory_cells_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(Microglia + Macrophages + Monocytes ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.05))
beeswarm(Microglia + Macrophages + Monocytes ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.05))
beeswarm(Microglia + Macrophages + Monocytes ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.05))
beeswarm(Microglia + Macrophages + Monocytes ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.05))
beeswarm(Microglia + Macrophages + Monocytes ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.05))
par(mfrow=c(1,1))
dev.off()

ciber.liver$age <- c(rep(3,3),rep(12,3),rep(29,3))
ciber.heart$age <- c(rep(3,3),rep(12,3),rep(29,3))
ciber.NPCs$age <- c(rep(3,2),rep(12,2),rep(29,2))
ciber.OB$age <- c(rep(3,3),rep(12,2),rep(29,3))
ciber.Cereb$age <- c(rep(3,3),rep(12,3),rep(29,3))

cereb.lm <- lm(Microglia + Macrophages+ Monocytes~ age, data = ciber.Cereb)
summary(cereb.lm) # p-value: 0.08821
heart.lm <- lm(Microglia + Macrophages+ Monocytes~ age, data = ciber.heart)
summary(heart.lm) # p-value: 0.1522
liver.lm <- lm(Microglia + Macrophages + Monocytes~ age, data = ciber.liver)
summary(liver.lm) # p-value: 0.9774
OB.lm <- lm(Microglia + Macrophages+ Monocytes~ age, data = ciber.OB)
summary(OB.lm) # p-value: 0.1791
npcs.lm <- lm(Microglia + Macrophages+ Monocytes~ age, data = ciber.NPCs)
summary(npcs.lm) # p-value: 0.8195


wilcox.test(c(ciber.Cereb$Microglia + ciber.Cereb$Macrophages + ciber.Cereb$Monocytes)[ciber.Cereb$age %in% '3'],
            c(ciber.Cereb$Microglia + ciber.Cereb$Macrophages + ciber.Cereb$Monocytes)[ciber.Cereb$age %in% '29'],
            alternative = "less")
#p-value = 0.2

wilcox.test(c(ciber.heart$Microglia + ciber.heart$Macrophages + ciber.heart$Monocytes)[ciber.heart$age %in% '3'],
            c(ciber.heart$Microglia + ciber.heart$Macrophages + ciber.heart$Monocytes)[ciber.heart$age %in% '29'],
            alternative = "less")
#p-value = 0.1

wilcox.test(c(ciber.liver$Microglia + ciber.liver$Macrophages + ciber.liver$Monocytes)[ciber.liver$age %in% '3'],
            c(ciber.liver$Microglia + ciber.liver$Macrophages + ciber.liver$Monocytes)[ciber.liver$age %in% '29'],
            alternative = "less")
#p-value = 0.65

wilcox.test(c(ciber.OB$Microglia + ciber.OB$Macrophages + ciber.OB$Monocytes)[ciber.OB$age %in% '3'],
            c(ciber.OB$Microglia + ciber.OB$Macrophages + ciber.OB$Monocytes)[ciber.OB$age %in% '29'],
            alternative = "less")
#p-value = 0.1

############  Lymphocytes
my.cols <- c("coral","blueviolet","dodgerblue")

pdf("CIBERSORT_B_cells_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(B_cells ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.01))
par(mfrow=c(1,1))
dev.off()


wilcox.test(c(ciber.Cereb$B_cells)[ciber.Cereb$age %in% '3'],
            c(ciber.Cereb$B_cells)[ciber.Cereb$age %in% '29'],
            alternative = "less")
#p-value = 0.908

wilcox.test(c(ciber.heart$B_cells)[ciber.heart$age %in% '3'],
            c(ciber.heart$B_cells)[ciber.heart$age %in% '29'],
            alternative = "less")
#p-value = 0.05

wilcox.test(c(ciber.liver$B_cells)[ciber.liver$age %in% '3'],
            c(ciber.liver$B_cells)[ciber.liver$age %in% '29'],
            alternative = "less")
#p-value = 0.1768

wilcox.test(c(ciber.OB$B_cells)[ciber.OB$age %in% '3'],
            c(ciber.OB$B_cells)[ciber.OB$age %in% '29'],
            alternative = "less")
#p-value = 0.2398


pdf("CIBERSORT_T_cells_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(T_cells ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.01))
beeswarm(T_cells ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.01))
beeswarm(T_cells ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.005))
beeswarm(T_cells ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.01))
beeswarm(T_cells ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.01))
par(mfrow=c(1,1))
dev.off()


pdf("CIBERSORT_NK_cells_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(NK_cells ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.01))
par(mfrow=c(1,1))
dev.off()


########## Fibrosis?
pdf("CIBERSORT_Fibroblasts_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(Cardiac_Fibroblasts + Dermal_Fibroblasts ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.05))
beeswarm(Cardiac_Fibroblasts + Dermal_Fibroblasts ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.05))
beeswarm(Cardiac_Fibroblasts + Dermal_Fibroblasts ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.05))
beeswarm(Cardiac_Fibroblasts + Dermal_Fibroblasts ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.05))
beeswarm(Cardiac_Fibroblasts + Dermal_Fibroblasts ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.05))
par(mfrow=c(1,1))
dev.off()


ciber.liver$age <- c(rep(3,3),rep(12,3),rep(29,3))
ciber.heart$age <- c(rep(3,3),rep(12,3),rep(29,3))
ciber.NPCs$age <- c(rep(3,2),rep(12,2),rep(29,2))
ciber.OB$age <- c(rep(3,3),rep(12,2),rep(29,3))
ciber.Cereb$age <- c(rep(3,3),rep(12,3),rep(29,3))

cereb.lm <- lm(Cardiac_Fibroblasts + Dermal_Fibroblasts~ age, data = ciber.Cereb)
summary(cereb.lm) # p-value: 0.2676
heart.lm <- lm(Cardiac_Fibroblasts + Dermal_Fibroblasts~ age, data = ciber.heart)
summary(heart.lm) # p-value: 0.1884
liver.lm <- lm(Cardiac_Fibroblasts + Dermal_Fibroblasts~ age, data = ciber.liver)
summary(liver.lm) # p-value: 0.07582
OB.lm <- lm(Cardiac_Fibroblasts + Dermal_Fibroblasts~ age, data = ciber.OB)
summary(OB.lm) # p-value: 0.8968
npcs.lm <- lm(Cardiac_Fibroblasts + Dermal_Fibroblasts~ age, data = ciber.NPCs)
summary(npcs.lm) # p-value: 0.2289


wilcox.test(ciber.heart$Cardiac_Fibroblasts[7:9],ciber.heart$Cardiac_Fibroblasts[1:3],alternative='greater')
#W = 7, p-value = 0.2
wilcox.test(ciber.liver$Dermal_Fibroblasts[7:9],ciber.liver$Dermal_Fibroblasts[1:3],alternative='greater')
#W = 7, p-value = 0.2

##############################################################################################################################
# 2017-04-05
# Disease models

hippo.mix <- read.csv('/Volumes/MyBook_3/BD_aging_project/RNAseq/All_tissues_analysis/CIBERSORT/Disease_models/Hippo_40Hz_CIBERSORT.Output_Job57.csv',row.names=1)
alz.mix <- read.csv('/Volumes/MyBook_3/BD_aging_project/RNAseq/All_tissues_analysis/CIBERSORT/Disease_models/Alzh_CIBERSORT.Output_Job58.csv',row.names=1)

hippo.mix$condition <- factor(c(rep("CTL",3),rep("40hz",3)))
alz.mix$condition <- factor(c(rep("CTL",6),rep("CKp25",6)))

my.cols <- rev(c("coral","dodgerblue"))

pdf("CIBERSORT_Inflammatory_cells_beeswarm_Hippocampus.pdf", width = 9, height = 6)
beeswarm(Microglia + Macrophages + Monocytes ~ condition, data = hippo.mix,pch = 16, col = my.cols, cex = 2, main = "Hippocampus 40Hz", ylim = c(0,0.2))
dev.off()

wilcox.test(c(hippo.mix$Microglia + hippo.mix$Macrophages + hippo.mix$Monocytes)[hippo.mix$condition %in% 'CTL'],
            c(hippo.mix$Microglia + hippo.mix$Macrophages + hippo.mix$Monocytes)[!(hippo.mix$condition %in% 'CTL')],
            alternative = "less")
#W = 0, p-value = 0.05 (hypothesis: there is more microglia in )



pdf("CIBERSORT_Inflammatory_cells_beeswarm_alzheimer.pdf", width = 9, height = 6)
beeswarm(Microglia + Macrophages + Monocytes ~ condition, data = alz.mix,pch = 16, col = my.cols, cex = 2, main = "alzheimer all", ylim = c(0,0.3))
dev.off()

wilcox.test(c(alz.mix$Microglia + alz.mix$Macrophages + alz.mix$Monocytes)[alz.mix$condition %in% 'CTL'],
            c(alz.mix$Microglia + alz.mix$Macrophages + alz.mix$Monocytes)[!(alz.mix$condition %in% 'CTL')],
            alternative = "less")
#p-value = 0.02056; all ages



############  Lymphocytes
my.cols <- c("coral","blueviolet","dodgerblue")

pdf("CIBERSORT_B_cells_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(B_cells ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.01))
beeswarm(B_cells ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.01))
par(mfrow=c(1,1))
dev.off()


pdf("CIBERSORT_T_cells_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(T_cells ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.01))
beeswarm(T_cells ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.01))
beeswarm(T_cells ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.005))
beeswarm(T_cells ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.01))
beeswarm(T_cells ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.01))
par(mfrow=c(1,1))
dev.off()


pdf("CIBERSORT_NK_cells_beeswarm.pdf", width = 9, height = 6)
par(mfrow=c(2,3))
beeswarm(NK_cells ~ age, data = ciber.Cereb,pch = 16, col = my.cols, cex = 2, main = "Cerebellum", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.heart,pch = 16, col = my.cols, cex = 2, main = "Heart", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.liver,pch = 16, col = my.cols, cex = 2, main = "Liver", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.OB,pch = 16, col = my.cols, cex = 2, main = "Olfactory Bulb", ylim = c(0,0.01))
beeswarm(NK_cells ~ age, data = ciber.NPCs, pch = 16, col = my.cols, cex = 2, main = "NPCs", ylim = c(0,0.01))
par(mfrow=c(1,1))
dev.off()

