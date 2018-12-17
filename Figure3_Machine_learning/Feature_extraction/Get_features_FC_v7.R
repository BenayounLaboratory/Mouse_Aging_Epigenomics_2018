setwd('/Volumes/MyBook_3/BD_aging_project/Machine_learning_aging/Predict_Fold_change')
options(stringsAsFactors = FALSE)

source('Features_extraction_functions_FC_v6.R')

# This code uses a large number of input features
# This code will take a couple of hours to run throughout

##############################################################################################
##################                1. get basic features                     ##################
##############################################################################################

##################     Liver     ##################
my.genes.liver <- get_RNA_features('Feature_folders/RNA_FPKM/2015-12-08_LIVER_FPKM_matrix_forML.txt',"Liver")

load("Feature_folders/RNAseq_DEseq2_results/RNA_seq_result_Liver_2015-11-19.RData")
get_DE_features(my.genes.liver,my.liver.RNAseq.process,"Liver",my.fdr = 0.1)

get_breadth_features_slopes("Liver", my.genes.liver, 
                            'Feature_folders/H3K4me3_breadth_data/Merged_ALL_AGES_MERGED_Liver_H3K4me3.PARSED_INTERSECTIONS.xls',
                            my.ages = c(3,3,12,12,29,29))

get_age_slopes("H3K4me3", "Liver", 
               'Feature_folders/PROMOTER_K4me3/Liver_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.liver, c(3,3,12,12,29,29))
               
get_age_slopes("H3K27ac", "Liver", 
               'Feature_folders/PROMOTER_K27ac/Liver_K27ac_cor_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.liver, c(3,3,12,12,29,29))

get_Nuc_features("Liver", my.genes.liver, 
                 'Feature_folders/Nucleosome_density_change/HOMER_Liver_DiNUP_DANPOS_GAINED.xls', 
                 'Feature_folders/Nucleosome_density_change/HOMER_Liver_DiNUP_DANPOS_LOST.xls',
                 'Feature_folders/Nucleosome_density_change/HOMER_2016-11-17_Liver_DiNUP_CHANGED_nucleosomes_FDR0.05.xls'
                 )

get_SELite_stat("Liver", my.genes.liver, 
                    'Feature_folders/Super_Enhancer_files/Annotated_3m1_Liver_H3K27ac.SuperEnhancers-Lite.xls',
                    'Feature_folders/Super_Enhancer_files/Annotated_3m2_Liver_H3K27ac.SuperEnhancers-Lite.xls')
get_SE_max_Score_feature("Liver", my.genes.liver,
                         'Feature_folders/Stitched_Enhancers_coverage/STITCHED_MACS2_ENHANCERS.Liver_SampleCoverage_Combined.xls',
                         c(3,3,12,12,29,29))


##################     Heart     ##################
my.genes.heart <- get_RNA_features('Feature_folders/RNA_FPKM/2015-12-08_HEART_FPKM_matrix_forML.txt',"Heart")

load("Feature_folders/RNAseq_DEseq2_results/RNA_seq_result_Heart_2015-11-19.RData")
get_DE_features(my.genes.heart,my.heart.RNAseq.process,"Heart",my.fdr = 0.1)

get_breadth_features_slopes("Heart", my.genes.heart, 
                            'Feature_folders/H3K4me3_breadth_data/Merged_ALL_AGES_MERGED_Heart_H3K4me3.PARSED_INTERSECTIONS.xls',
                            my.ages = c(3,3,12,12,29,29))

get_age_slopes("H3K4me3", "Heart", 
               'Feature_folders/PROMOTER_K4me3/Heart_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.heart, c(3,3,12,12,29,29))
get_age_slopes("H3K27ac", "Heart", 
               'Feature_folders/PROMOTER_K27ac/Heart_K27ac_cor_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.heart, c(3,3,12,12,29,29))

get_Nuc_features("Heart", my.genes.heart, 
                 'Feature_folders/Nucleosome_density_change/HOMER_Heart_DiNUP_DANPOS_GAINED.xls', 
                 'Feature_folders/Nucleosome_density_change/HOMER_Heart_DiNUP_DANPOS_LOST.xls',
                 'Feature_folders/Nucleosome_density_change/HOMER_2016-11-17_Heart_DiNUP_CHANGED_nucleosomes_FDR0.05.xls'
)

get_SELite_stat("Heart", my.genes.heart, 
                    'Feature_folders/Super_Enhancer_files/Annotated_3m1_Heart_H3K27ac.SuperEnhancers-Lite.xls',
                    'Feature_folders/Super_Enhancer_files/Annotated_3m2_Heart_H3K27ac.SuperEnhancers-Lite.xls')
get_SE_max_Score_feature("Heart", my.genes.heart, 
                         'Feature_folders/Stitched_Enhancers_coverage/STITCHED_MACS2_ENHANCERS.Heart_SampleCoverage_Combined.xls',
                         c(3,3,12,12,29,29))


##################     cereb     ##################
my.genes.cereb <- get_RNA_features('Feature_folders/RNA_FPKM/2015-12-08_CEREBELLUM_FPKM_matrix_forML.txt',"Cerebellum")

load("Feature_folders/RNAseq_DEseq2_results/RNA_seq_result_cereb_2015-11-19.RData")
get_DE_features(my.genes.cereb,my.cereb.RNAseq.process,"Cerebellum",my.fdr = 0.1)

get_breadth_features_slopes("Cerebellum", my.genes.cereb, 
                            'Feature_folders/H3K4me3_breadth_data/Merged_ALL_AGES_MERGED_Cerebellum_H3K4me3.PARSED_INTERSECTIONS.xls',
                            my.ages = c(3,3,12,12,29,29))

get_age_slopes("H3K4me3", "Cerebellum", 
               'Feature_folders/PROMOTER_K4me3/Cerebellum_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.cereb, c(3,3,12,12,29,29))
get_age_slopes("H3K27ac", "Cerebellum", 
               'Feature_folders/PROMOTER_K27ac/Cerebellum_K27ac_cor_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.cereb, c(3,3,12,12,29,29))


get_Nuc_features("Cerebellum", my.genes.cereb, 
                 'Feature_folders/Nucleosome_density_change/HOMER_Cerebellum_DiNUP_DANPOS_GAINED.xls', 
                 'Feature_folders/Nucleosome_density_change/HOMER_Cerebellum_DiNUP_DANPOS_LOST.xls',
                 'Feature_folders/Nucleosome_density_change/HOMER_2017-03-20_Cerebellum_DiNUP_CHANGED_nucleosomes_FDR0.05.xls'
)


get_SELite_stat("Cerebellum", my.genes.cereb, 
                    'Feature_folders/Super_Enhancer_files/Annotated_3m2_Cerebellum_H3K27ac.SuperEnhancers-Lite.xls',
                    'Feature_folders/Super_Enhancer_files/Annotated_3m1_Cerebellum_H3K27ac.SuperEnhancers-Lite.xls')
get_SE_max_Score_feature("Cerebellum", my.genes.cereb, 
                         'Feature_folders/Stitched_Enhancers_coverage/STITCHED_MACS2_ENHANCERS.Cerebellum_SampleCoverage_Combined.xls',
                         c(3,3,12,12,29,29))


##################     OB     ##################
my.genes.ob <- get_RNA_features('Feature_folders/RNA_FPKM/2015-12-08_OLFACTORY_BULB_FPKM_matrix_forML.txt',
                                   "OB", my.3=1:3, my.12=4:5, my.29=6:8)

load("Feature_folders/RNAseq_DEseq2_results/RNA_seq_result_OB_2015-11-19.RData")
get_DE_features(my.genes.ob,my.ob.RNAseq.process, "OB",my.fdr = 0.1)

get_breadth_features_slopes("OB", my.genes.ob, 
                            'Feature_folders/H3K4me3_breadth_data/Merged_ALL_AGES_MERGED_OB_H3K4me3.PARSED_INTERSECTIONS.xls',
                            my.ages = c(3,3,12,12,29,29))

get_age_slopes("H3K4me3", "OB", 
               'Feature_folders/PROMOTER_K4me3/OB_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.ob, c(3,3,29,29,29))
get_age_slopes("H3K27ac", "OB", 
               'Feature_folders/PROMOTER_K27ac/OB_K27ac_cor_normalized_promoter_densities_NO_HEADER.txt',
               my.genes.ob, c(3,3,29,29,29))

get_Nuc_features("OB", my.genes.ob, 
                 'Feature_folders/Nucleosome_density_change/HOMER_OlfactoryBulb_DiNUP_DANPOS_GAINED.xls', 
                 'Feature_folders/Nucleosome_density_change/HOMER_OlfactoryBulb_DiNUP_DANPOS_LOST.xls',
                 'Feature_folders/Nucleosome_density_change/HOMER_2016-11-17_OlfactoryBulb_DiNUP_CHANGED_nucleosomes_FDR0.05.xls'
)

get_SELite_stat("OB", my.genes.ob, 
                    'Feature_folders/Super_Enhancer_files/Annotated_3m1_Olfactory_Bulb_H3K27ac.SuperEnhancers-Lite.xls',
                    'Feature_folders/Super_Enhancer_files/Annotated_3m2_Olfactory_Bulb_H3K27ac.SuperEnhancers-Lite.xls')
get_SE_max_Score_feature("OB", my.genes.ob,
                         'Feature_folders/Stitched_Enhancers_coverage/STITCHED_MACS2_ENHANCERS.Olfactory.Bulb_SampleCoverage_Combined.xls',
                         c(3,3,12,12,29,29))

####
save(my.genes.liver, my.genes.heart, my.genes.cereb, my.genes.ob, file="2016-11-16_gene_lists_per_tissue.RData")


##############################################################################################
##################            2. read in "young only" features              ##################
##############################################################################################

##### Pol2 peak number per gene in young
get_Pol2_numbers("Liver", my.genes.liver, "Feature_folders/Pol2_peaks/HOMER_liver_8weeks_Pol2_peaks.xls")
get_Pol2_numbers("Heart", my.genes.heart, "Feature_folders/Pol2_peaks/HOMER_heart_8weeks_Pol2_peaks.xls")
get_Pol2_numbers("Cerebellum", my.genes.cereb, "Feature_folders/Pol2_peaks/HOMER_cerebellum_8weeks_Pol2_peaks.xls")
get_Pol2_numbers("OB", my.genes.ob, "Feature_folders/Pol2_peaks/HOMER_olfactoryBulb_8weeks_Pol2_peaks.xls")
    #### genes with no pol2 peaks are marked as 0 peaks, 500000 abs dist, and 0 max score


##### promoter intensity (all genes, get highest TR for each gene)
get_prom_features("Feature_folders/2016-01-20_Tissues_promoter_Pol2_normalized_promoter_densities_NO_HEADER.txt",
                  "Pol2", c("Cerebellum","Heart","Liver","OB","NPCs"))
get_prom_features("Feature_folders/2016-01-20_Tissues_promoter_H3K4me1_normalized_promoter_densities_NO_HEADER.txt",
                  "H3K4me1", c("Cerebellum","Heart","Liver","OB","NPCs"))
get_prom_features("Feature_folders/2016-01-20_Tissues_promoter_CTCF_normalized_promoter_densities_NO_HEADER.txt",
                  "CTCF", c("Cerebellum","Heart","Liver","OB","NPCs"))
get_prom_features("Feature_folders/2016-11-15_Tissues_promoter_Accessibility_normalized_promoter_densities_NO_HEADER.txt",
                  "DNAse-ATAC", c("Cerebellum","Heart","Liver","OB","NPCs"))
get_prom_features("Feature_folders/2016-11-16_Tissues_promoter_H3K27me3_normalized_promoter_densities_NO_HEADER.txt",
                  "H3K27me3", c("Cerebellum","Heart","Liver","OB","NPCs"))


##### traveling ratio
get_TR_features("Feature_folders/Pol2_TR/2016-1-22_Liver_mm9_TravelingRatios.txt", my.genes.liver, "Liver")
get_TR_features("Feature_folders/Pol2_TR/2016-1-22_Heart_mm9_TravelingRatios.txt", my.genes.heart, "Heart")
get_TR_features("Feature_folders/Pol2_TR/2016-1-22_Cerebellum_mm9_TravelingRatios.txt", my.genes.cereb, "Cerebellum")
get_TR_features("Feature_folders/Pol2_TR/2016-1-22_OB_mm9_TravelingRatios.txt", my.genes.ob, "OB")


##### SE distance to TSS
my.cereb.SE.dists <- get_SE_dist('Feature_folders/Meta_Super_Enhancers/HOMER_ALL_AGES_MERGED_Cerebellum_H3K27ac_broad_peaks.stitched.xls',
                                 'Feature_folders/Meta_Super_Enhancers/Annotated_ALL_AGES_MERGED_Cerebellum_H3K27ac_broad_peaksWITH-SELite-SCORE.xls',
                                 "Cerebellum")
my.liver.SE.dists <- get_SE_dist('Feature_folders/Meta_Super_Enhancers/HOMER_ALL_AGES_MERGED_Liver_H3K27ac_broad_peaks.stitched.xls',
                                 'Feature_folders/Meta_Super_Enhancers/Annotated_ALL_AGES_MERGED_Liver_H3K27ac_broad_peaksWITH-SELite-SCORE.xls',
                                 "Liver")
my.OB.SE.dists <- get_SE_dist('Feature_folders/Meta_Super_Enhancers/HOMER_ALL_AGES_MERGED_OB_H3K27ac_PooledLanes_broad_peaks.stitched.xls',
                              'Feature_folders/Meta_Super_Enhancers/Annotated_ALL_AGES_MERGED_OB_H3K27ac_PooledLanes_broad_peaksWITH-SELite-SCORE.xls',
                              "OlfactoryBulb")
my.heart.SE.dists <- get_SE_dist('Feature_folders/Meta_Super_Enhancers/HOMER_ALL_AGES_MERGED_Heart_H3K27ac_broad_peaks.stitched.xls',
                                 'Feature_folders/Meta_Super_Enhancers/Annotated_ALL_AGES_MERGED_Heart_H3K27ac_broad_peaksWITH-SELite-SCORE.xls',
                                 "Heart")

save(my.cereb.SE.dists, my.liver.SE.dists, my.OB.SE.dists, my.heart.SE.dists, file=paste(Sys.Date(),"All_SE_TSS_distance_feature_objects.RData",sep="_"))

##### Bivalent domains
get_bivalent("Liver", my.genes.liver, 'Feature_folders/Bivalent/HOMER_Bivalent_domains_Liver.xls')
get_bivalent("Heart", my.genes.heart, 'Feature_folders/Bivalent/HOMER_Bivalent_domains_Heart.xls')
get_bivalent("Cerebellum", my.genes.cereb, 'Feature_folders/Bivalent/HOMER_Bivalent_domains_Cerebellum.xls')
get_bivalent("OB", my.genes.ob, 'Feature_folders/Bivalent/HOMER_Bivalent_domains_OB_vs_WholeBrain.xls')


##### CTCF distance
get_CTCF_feats("Liver", my.genes.liver, 'Feature_folders/BCTCF_dists/HOMER_CTCF-liver_8weeks.xls')
get_CTCF_feats("Heart", my.genes.heart, 'Feature_folders/BCTCF_dists/HOMER_CTCF-heart_8weeks.xls')
get_CTCF_feats("Cerebellum", my.genes.cereb, 'Feature_folders/BCTCF_dists/HOMER_CTCF-cerebellum_8weeks.xls')
get_CTCF_feats("OB", my.genes.ob, 'Feature_folders/BCTCF_dists/HOMER_CTCF-olfactory-bulb.xls')

##############################################################################################
##################              3. process constant features                ##################
##############################################################################################

my.homer.prom.motifs <- read.csv("Feature_folders/Constant_features/2016-01-19_parsed_for_ML_mm9_masked_Homer_genome_wide_motifs.txt", header=T,sep="\t")

# alternative TSSs will create some issues -> pick the first one of each (duplicates)
my.genes <- unique(my.homer.prom.motifs$GeneName)

my.homer.prom.motifs.v2 <- data.frame(matrix(NA,length(my.genes),dim(my.homer.prom.motifs)[2]))
colnames(my.homer.prom.motifs.v2) <- colnames(my.homer.prom.motifs)

for ( i in 1:length(my.genes)) {
  my.first.ix <- (which(my.homer.prom.motifs$GeneName %in% my.genes[i]))[1]
  my.homer.prom.motifs.v2[i,] <- my.homer.prom.motifs[my.first.ix,]
  
}
#### Supplied to run functions, but TFBS features are not used in paper!!!!



my.homer.prom.CG <- read.csv("Feature_folders/Constant_features/mm9_masked_Homer_genome_wide_promoter_CG_CpG.homer.txt", header=T,sep="\t")
colnames(my.homer.prom.CG)[1] <- "GeneName"

my.homer.prom.CG.v2 <- data.frame(matrix(NA,length(my.genes),dim(my.homer.prom.CG)[2]))
colnames(my.homer.prom.CG.v2) <- colnames(my.homer.prom.CG)

for ( i in 1:length(my.genes)) {
  my.first.ix <- (which(my.homer.prom.motifs$GeneName %in% my.genes[i]))[1]
  my.homer.prom.CG.v2[i,] <- my.homer.prom.CG[my.first.ix,]
  
}

save(my.homer.prom.motifs.v2,my.homer.prom.CG.v2, file = "2016-10-05_HOMER_motifs_CG_features.RData")


######### CpG Islands
my.genes.all <- unique(c(my.genes.liver,my.genes.heart,my.genes.cereb,my.genes.ob))

my.ucsc.cpg <- read.csv("Feature_folders/Constant_features/HOMER_2016-11-07_UCSC_mm9_CpG_Islands.xls", header=T,sep="\t")
my.ucsc.cpg.v2 <- aggregate(Distance.to.TSS ~ Gene.Name, data = my.ucsc.cpg, FUN = function(x){NROW(x)})
colnames(my.ucsc.cpg.v2)[2] <- "CpG_islands"

# add non CpGs
my.non.cpg <- which(!(my.genes.all %in% my.ucsc.cpg.v2$Gene.Name))
my.non.cpg.data <- cbind(my.genes.all[my.non.cpg],rep(0,length(my.non.cpg)))
colnames(my.non.cpg.data) <- colnames(my.ucsc.cpg.v2)
my.ucsc.cpg.all <- rbind(my.ucsc.cpg.v2,my.non.cpg.data)

save(my.ucsc.cpg.all, file = "2016-11-16_UCSC_CpG_islands.RData")


####### Exons number per gene
my.genes.all <- unique(c(my.genes.liver,my.genes.heart,my.genes.cereb,my.genes.ob))

my.exon.data <- read.csv("Feature_folders/Constant_features/2016-09-29_Gene_exons_for_ML_Biomart_Ens85.txt", header=T,sep="\t")
my.exon.data.v2 <- aggregate(Constitutive.Exon ~ Ensembl.Gene.ID, data = my.exon.data, FUN = function(x){NROW(x)})
my.exon.data.ann <- unique(my.exon.data[,c(1,8)])
my.exon.data.final <- (merge(my.exon.data.ann,my.exon.data.v2,by='Ensembl.Gene.ID'))[,-1]

my.exon.data.final.v2 <- my.exon.data.final[my.exon.data.final$Associated.Gene.Name %in% my.genes.all,]
  
save(my.exon.data.final.v2, file = "2016-11-16_ENSEMBL_exon_numbers.RData")


##############################################################################################
##################            4. build feature matrux by tissue             ##################
##############################################################################################

#!!!!!!! Names of RData files NEED to be updated when code is independently rerun

### features for all tissues
load('./Extracted_features_RData/2016-10-05_CTCF_prom_feature_object.RData')
my.ctcf.data <- my.prom.av
my.ctcf.data$GeneName <- rownames(my.ctcf.data)
rm(my.prom.av)

load('./Extracted_features_RData/2016-10-05_H3K4me1_prom_feature_object.RData')
my.k4me1.data <- my.prom.av
my.k4me1.data$GeneName <- rownames(my.k4me1.data)
rm(my.prom.av)

load('./Extracted_features_RData/2016-10-05_Pol2_prom_feature_object.RData')
my.Pol2.data <- my.prom.av
my.Pol2.data$GeneName <- rownames(my.Pol2.data)
rm(my.prom.av)

load('./Extracted_features_RData/2016-11-16_DNAse-ATAC_prom_feature_object.RData')
my.access.data <- my.prom.av
my.access.data$GeneName <- rownames(my.access.data)
rm(my.prom.av)

load('./Extracted_features_RData/2016-11-16_H3K27me3_prom_feature_object.RData')
my.K27me3.data <- my.prom.av
my.K27me3.data$GeneName <- rownames(my.K27me3.data)
rm(my.prom.av)


#### load age/cell independent data
# GG/CpG
load("./Extracted_features_RData/2016-10-05_HOMER_motifs_CG_features.RData")
load("./Extracted_features_RData/2016-11-16_UCSC_CpG_islands.RData") #my.ucsc.cpg.all

# TF targets (CheA/ENCODE/GEO)
load("../EnrichR_data/2016-08-03_TF_targets_feature_matrices.RData") #### Supplied to run function, but features are not used in paper

# SE distance
load('Extracted_features_RData/2016-11-16_All_SE_TSS_distance_feature_objects.RData') ### SE dist

# exon number
load("./Extracted_features_RData/2016-11-16_ENSEMBL_exon_numbers.RData") # my.exon.data.final.v2



## Heart ##
load('./Extracted_features_RData/2016-10-04_Heart__SE_categorical_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_Heart__SE_scores_feature_object.RData')
load('./Extracted_features_RData/2016-11-17_Heart_Differential_Nucleosomes_3vs29m_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_Heart_FOR_FC_feature_object.RData')
my.exp.feat.mat <- my.feat.mat[,-4]
load('./Extracted_features_RData/2016-10-05_Heart_K4breadth_aging_with_ALL_features_object.RData')
load('./Extracted_features_RData/2016-10-04_Heart_RNA_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_Heart_Pol2_young_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_Heart__TravelingRatio_feature_object.RData')
load('./Extracted_features_RData/2016-11-17_Heart_H3K4me3_chromatin_aging_feature_object.RData')
my.K4me3.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-17_Heart_H3K27ac_chromatin_aging_feature_object.RData')
my.K27ac.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-16_Heart_CTCF_young_feature_object.RData')
load('./Extracted_features_RData/2016-11-16_Heart_Bivalency_young_feature_object.RData')

my.heart.features <- get_feat_matsv4("Heart",
                                     my.rna.med,
                                     my.exp.feat.mat,
                                     my.K4me3.feat.mat,
                                     my.K27ac.feat.mat,
                                     my.SE,
                                     my.SE.score,
                                     my.heart.SE.dists,
                                     my.K4breadth.max,
                                     my.tr.av,
                                     my.nuc.feats,
                                     my.k4me1.data,
                                     my.ctcf.data,
                                     my.Pol2.data,
                                     my.Pol2.feats,
                                     my.ctcf.feats,
                                     my.exon.data.final.v2,
                                     my.ucsc.cpg.all,
                                     my.homer.prom.CG.v2,
                                     my.ChIP.feat.matrix,
                                     my.LOF.feat.matrix)

my.heart.features.v2 <- get_additional_features("Heart",
                                                my.heart.features,
                                                my.K27me3.data,
                                                my.access.data,
                                                my.biv.feats)

rm(my.rna.med,my.exp.feat.mat,my.K4me3.feat.mat,my.K27ac.feat.mat,
   my.SE,my.SE.score,my.K4breadth.max,my.tr.av,my.Pol2.feats,my.nuc.feats,
   my.biv.feats,my.ctcf.feats)

## Liver ##
load('./Extracted_features_RData/2016-10-04_Liver__SE_categorical_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_Liver__SE_scores_feature_object.RData')
load('./Extracted_features_RData/2016-11-17_Liver_Differential_Nucleosomes_3vs29m_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_Liver_FOR_FC_feature_object.RData')
my.exp.feat.mat <- my.feat.mat[,-4]
load('./Extracted_features_RData/2016-10-05_Liver_K4breadth_aging_with_ALL_features_object.RData')
load('./Extracted_features_RData/2016-10-04_Liver_RNA_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_Liver_Pol2_young_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_Liver__TravelingRatio_feature_object.RData')
load('./Extracted_features_RData/2016-11-17_Liver_H3K4me3_chromatin_aging_feature_object.RData')
my.K4me3.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-17_Liver_H3K27ac_chromatin_aging_feature_object.RData')
my.K27ac.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-16_Liver_CTCF_young_feature_object.RData')
load('./Extracted_features_RData/2016-11-16_Liver_Bivalency_young_feature_object.RData')


my.liver.features <- get_feat_matsv4("Liver",
                                     my.rna.med,
                                     my.exp.feat.mat,
                                     my.K4me3.feat.mat,
                                     my.K27ac.feat.mat,
                                     my.SE,
                                     my.SE.score,
                                     my.liver.SE.dists,
                                     my.K4breadth.max,
                                     my.tr.av,
                                     my.nuc.feats,
                                     my.k4me1.data,
                                     my.ctcf.data,
                                     my.Pol2.data,
                                     my.Pol2.feats,
                                     my.ctcf.feats,
                                     my.exon.data.final.v2,
                                     my.ucsc.cpg.all,
                                     my.homer.prom.CG.v2,
                                     my.ChIP.feat.matrix,
                                     my.LOF.feat.matrix)

my.liver.features.v2 <- get_additional_features("Liver",
                                                my.liver.features,
                                                my.K27me3.data,
                                                my.access.data,
                                                my.biv.feats)

rm(my.rna.med,my.exp.feat.mat,my.K4me3.feat.mat,my.K27ac.feat.mat,
   my.SE,my.SE.score,my.K4breadth.max,my.tr.av,my.Pol2.feats,my.nuc.feats,
   my.biv.feats,my.ctcf.feats)

## Cerebellum ##
load('./Extracted_features_RData/2016-10-04_Cerebellum__SE_categorical_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_Cerebellum__SE_scores_feature_object.RData')
load('./Extracted_features_RData/2017-03-20_Cerebellum_Differential_Nucleosomes_3vs29m_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_Cerebellum_FOR_FC_feature_object.RData')
my.exp.feat.mat <- my.feat.mat[,-4]
load('./Extracted_features_RData/2016-10-05_Cerebellum_K4breadth_aging_with_ALL_features_object.RData')
load('./Extracted_features_RData/2016-10-04_Cerebellum_RNA_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_Cerebellum_Pol2_young_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_Cerebellum__TravelingRatio_feature_object.RData')
load('./Extracted_features_RData/2016-11-17_Cerebellum_H3K4me3_chromatin_aging_feature_object.RData')
my.K4me3.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-17_Cerebellum_H3K27ac_chromatin_aging_feature_object.RData')
my.K27ac.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-16_Cerebellum_CTCF_young_feature_object.RData')
load('./Extracted_features_RData/2016-11-16_Cerebellum_Bivalency_young_feature_object.RData')

my.cereb.features <- get_feat_matsv4("Cerebellum",
                                     my.rna.med,
                                     my.exp.feat.mat,
                                     my.K4me3.feat.mat,
                                     my.K27ac.feat.mat,
                                     my.SE,
                                     my.SE.score,
                                     my.cereb.SE.dists,
                                     my.K4breadth.max,
                                     my.tr.av,
                                     my.nuc.feats,
                                     my.k4me1.data,
                                     my.ctcf.data,
                                     my.Pol2.data,
                                     my.Pol2.feats,
                                     my.ctcf.feats,
                                     my.exon.data.final.v2,
                                     my.ucsc.cpg.all,
                                     my.homer.prom.CG.v2,
                                     my.ChIP.feat.matrix,
                                     my.LOF.feat.matrix)

my.cereb.features.v2 <- get_additional_features("Cerebellum",
                                                my.cereb.features,
                                                my.K27me3.data,
                                                my.access.data,
                                                my.biv.feats)

rm(my.rna.med,my.exp.feat.mat,my.K4me3.feat.mat,my.K27ac.feat.mat,
   my.SE,my.SE.score,my.K4breadth.max,my.tr.av,my.Pol2.feats,my.nuc.feats,
   my.biv.feats,my.ctcf.feats)

## OB ##
load('./Extracted_features_RData/2016-10-04_OB__SE_categorical_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_OB__SE_scores_feature_object.RData')
load('./Extracted_features_RData/2016-11-17_OB_Differential_Nucleosomes_3vs29m_feature_object.RData')
load('./Extracted_features_RData/2016-10-04_OB_FOR_FC_feature_object.RData')
my.exp.feat.mat <- my.feat.mat[,-4]
load('./Extracted_features_RData/2016-10-05_OB_K4breadth_aging_with_ALL_features_object.RData')
load('./Extracted_features_RData/2016-10-04_OB_RNA_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_OB_Pol2_young_feature_object.RData')
load('./Extracted_features_RData/2016-10-05_OB__TravelingRatio_feature_object.RData')
load('./Extracted_features_RData/2016-11-17_OB_H3K4me3_chromatin_aging_feature_object.RData')
my.K4me3.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-17_OB_H3K27ac_chromatin_aging_feature_object.RData')
my.K27ac.feat.mat <- my.feat.mat
load('./Extracted_features_RData/2016-11-16_OB_CTCF_young_feature_object.RData')
load('./Extracted_features_RData/2016-11-16_OB_Bivalency_young_feature_object.RData')

my.ob.features <- get_feat_matsv4("OB",
                                  my.rna.med,
                                  my.exp.feat.mat,
                                  my.K4me3.feat.mat,
                                  my.K27ac.feat.mat,
                                  my.SE,
                                  my.SE.score,
                                  my.cereb.SE.dists,
                                  my.K4breadth.max,
                                  my.tr.av,
                                  my.nuc.feats,
                                  my.k4me1.data,
                                  my.ctcf.data,
                                  my.Pol2.data,
                                  my.Pol2.feats,
                                  my.ctcf.feats,
                                  my.exon.data.final.v2,
                                  my.ucsc.cpg.all,
                                  my.homer.prom.CG.v2,
                                  my.ChIP.feat.matrix,
                                  my.LOF.feat.matrix)

my.ob.features.v2 <- get_additional_features("OB",
                                             my.ob.features,
                                             my.K27me3.data,
                                             my.access.data,
                                             my.biv.feats)

rm(my.rna.med,my.exp.feat.mat,my.K4me3.feat.mat,my.K27ac.feat.mat,
   my.SE,my.SE.score,my.K4breadth.max,my.tr.av,my.Pol2.feats,my.nuc.feats,
   my.biv.feats,my.ctcf.feats)


#################
# save ALL feature matrices
save(my.heart.features, my.heart.features.v2,
     my.liver.features, my.liver.features.v2,my.liver.features.v3,
     my.cereb.features, my.cereb.features.v2,
     my.ob.features, my.ob.features.v2,
     file = paste(Sys.Date(),"Complete_feature_matrices_FOLD_CHANGE.RData",sep="_"))

save(my.cereb.features, my.cereb.features.v2,
     file = paste(Sys.Date(),"Complete_feature_matrices_CEREB_FOLD_CHANGE.RData",sep="_"))


my.heart.features <- rm_nas(my.heart.features)
my.heart.features.v2 <- rm_nas(my.heart.features.v2)
my.liver.features <- rm_nas(my.liver.features)
my.liver.features.v2 <- rm_nas(my.liver.features.v2)
my.cereb.features <- rm_nas(my.cereb.features)
my.cereb.features.v2 <- rm_nas(my.cereb.features.v2)
my.ob.features <- rm_nas(my.ob.features)
my.ob.features.v2 <- rm_nas(my.ob.features.v2)
my.liver.features.v3 <- rm_nas(my.liver.features.v3)

# save ALL feature matrices without NAs
save(my.heart.features,my.heart.features.v2,
     my.liver.features,my.liver.features.v2,my.liver.features.v3,
     my.cereb.features,my.cereb.features.v2,
     my.ob.features, my.ob.features.v2,
     file = paste(Sys.Date(),"Complete_feature_matrices_FOLD_CHANGE_NA_RM.RData",sep="_"))


my.cereb.features <- rm_nas(my.cereb.features)
my.cereb.features.v2 <- rm_nas(my.cereb.features.v2)
save(my.cereb.features, my.cereb.features.v2,
     file = paste(Sys.Date(),"Complete_feature_matrices_CEREB_FOLD_CHANGE_NA_RM.RData",sep="_"))

