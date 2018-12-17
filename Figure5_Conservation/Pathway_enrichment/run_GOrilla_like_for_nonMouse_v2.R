#setwd('/Volumes/BB_Backup_3//BD_aging_project/Public_datasets/RNA_profile_other_species/Pathway_enrichment')
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # PPS - setwd to current file location
library('mHG')
options(stringsAsFactors=F)

source('GOrilla_statistics_functions_for_nonMouse.R')

my.gmt.sets <- c("../../Figure4_enrichments/PATHWAYS/h.all.v5.1.symbols.gmt",
                 "../../Figure4_enrichments/PATHWAYS/M.musculus_KEGG-pathways_04-2017_NO_DISEASES_UC.gmt",
                 "../../Figure4_enrichments/PATHWAYS/TF-LOF_Expression_from_GEO_WITH_FOXO_2016-08-02_TF_targets.gmt")
                 
my.gmt.set.names <- c("MSIgDB_Hallmark_Datasets",
                    "KEGG_2017_no_diseases_UC",
                    "TF-LOF_Expression_from_GEO_WITH_FOXO_parsed-and-aggregated")
 
cbind(my.gmt.sets, my.gmt.set.names)

my.run <- 1:3


##################################################################################### 
#Load RData RNA with age results
load("../Comparisons/Rat/DEseq2_output/2017-05-19_Rat_aging_brain_RNAseq.RData")
load("../Comparisons/Rat/DEseq2_output/2017-05-19_Rat_aging_Heart_RNAseq.RData")
load("../Comparisons/Rat/DEseq2_output/2017-05-19_Rat_aging_Liver_RNAseq.RData")

load("../Comparisons/GTex/Output/2016-12-28_Cerebellum_GTEx_data_DEseq2_aging_genename_ONLYmale_no_batch.RData")
load("../Comparisons/GTex/Output/2016-12-28_Liver_GTEx_data_DEseq2_aging_genename_ONLYmale_no_batch.RData")
load("../Comparisons/GTex/Output/2016-12-29_Heart_GTEx_data_DEseq2_aging_genename_ONLYmale_no_batch.RData")

load("../Comparisons/Killifish_aging_RNAseq/Brain/Output/2017-05-15_Killifish_RNAseq_aging__statistics.RData")
my.brain.killi.process <- res; rm(res)
load("../Comparisons/Killifish_aging_RNAseq/Liver/Output/2017-05-15_Killifish_RNAseq_aging__statistics.RData")
my.liver.killi.process <- res; rm(res)

###########
# get ortholog names
# A. Killifish - use Param's BLAST results
my.orthology.killi <- read.csv('Orthology/BestHits_nfur-mmus_1e-3.txt', sep = "\t", header = F)
my.nfur.names <- unlist(lapply(strsplit(my.orthology.killi$V1,"|",fixed = T),get_first))
my.mouse.names <- unlist(lapply(strsplit(my.orthology.killi$V2,"|",fixed = T),get_first))
my.orth.table.killi <- data.frame(cbind(my.nfur.names,my.mouse.names))
colnames(my.orth.table.killi) <- c("Nfur_Symbol","Mouse_Symbol")

# B. Rat. Biomart
my.orthology.rat <- read.csv("Orthology/2017-05-15_Rat_mouse_biomart_orthology.txt", sep = "\t", header = T)
my.orth.table.rat <- unique(data.frame(cbind(my.orthology.rat$Gene.name,my.orthology.rat$Mouse.gene.name)))
colnames(my.orth.table.rat) <- c("Rat_Symbol","Mouse_Symbol")
my.orth.table.rat <- unique(my.orth.table.rat)

# C. Human Gencode/ Biomart
# read orthology on Gencode v19
my.orthology.human <- read.table("Orthology/2016-12-16_Correspondence_GeneName_Human_Mouse_Orthologs.txt",header=T,sep="\t")
my.orth.table.human <- unique(data.frame(cbind(my.orthology.human$Human_Symbol,my.orthology.human$Mouse_Symbol)))
colnames(my.orth.table.human) <- c("Human_Symbol","Mouse_Symbol")

##########
# My data
run_pathway_enrich("Heart_RAT", my.rat.heart.process, my.orth.table = my.orth.table.rat, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])
run_pathway_enrich("Liver_RAT", my.rat.liver.process, my.orth.table = my.orth.table.rat, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])
run_pathway_enrich("Brain_RAT", my.rat.brain.process, my.orth.table = my.orth.table.rat, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])

run_pathway_enrich("Heart_Human", my.heart.gtex.process, my.orth.table = my.orth.table.human, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])
run_pathway_enrich("Liver_Human", my.liver.gtex.process, my.orth.table = my.orth.table.human, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])
run_pathway_enrich("Cereb_Human", my.cereb.gtex.process, my.orth.table = my.orth.table.human, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])

run_pathway_enrich("Liver_Killifish", my.liver.killi.process, my.orth.table = my.orth.table.killi, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])
run_pathway_enrich("Brain_Killifish", my.brain.killi.process, my.orth.table = my.orth.table.killi, my.gmt.files = my.gmt.sets[my.run], my.set.names = my.gmt.set.names[my.run])

