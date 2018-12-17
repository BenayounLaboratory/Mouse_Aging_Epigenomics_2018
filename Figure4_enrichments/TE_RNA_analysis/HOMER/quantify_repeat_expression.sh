# 2017-12-11
# use homer to quantify repeats

RNATAGS='/Volumes/BB_Backup_3/BD_aging_project/RNAseq/All_tissues_analysis/Bedgraphs'

# over annotated repeats
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_3m1_TAGs          > Cereb_3m1_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_3m2_TAGs          > Cereb_3m2_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_3m3_TAGs          > Cereb_3m3_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_12m1_TAGs         > Cereb_12m1_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_12m2_TAGs         > Cereb_12m2_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_12m3_TAGs         > Cereb_12m3_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_29m2_TAGs         > Cereb_29m2_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_29m3_TAGs         > Cereb_29m3_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/Cereb_29m4_TAGs         > Cereb_29m4_RNA_repeats.txt 
#
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_3m1_ATCACG_TAGs      > OB_3m1_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_3m2_CGATGT_TAGs      > OB_3m2_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_3m3_TTAGGC_TAGs      > OB_3m3_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_12m1_TGACCA_TAGs     > OB_12m1_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_12m3_GCCAAT_TAGs     > OB_12m3_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_29m2_ACTTGA_TAGs     > OB_29m2_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_29m3_GATCAG_TAGs     > OB_29m3_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/OB_29m1_CAGATC_TAGs     > OB_29m1_RNA_repeats.txt 
#
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/NPC_3m5-CGATGT_TAGs     > NPC_3m5_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/NPC_3m6-TGACCA_TAGs     > NPC_3m6_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/NPC_12m5-ACAGTG_TAGs    > NPC_12m5_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/NPC_12m6-GCCAAT_TAGs    > NPC_12m6_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/NPC_29m6-CTTGTA_TAGs    > NPC_29m6_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/NPC_29m5-CAGATC_TAGs    > NPC_29m5_RNA_repeats.txt 
#
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/3m4_heart_ATCACG_TAGs   > 3m4_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/3m5_heart_CGATGT_TAGs   > 3m5_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/3m6_heart_TTAGGC_TAGs   > 3m6_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/12m4_heart_TGACCA_TAGs  > 12m4_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/12m5_heart_ACAGTG_TAGs  > 12m5_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/12m6_heart_GCCAAT_TAGs  > 12m6_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/29m4_heart_CAGATC_TAGs  > 29m4_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/29m5_heart_ACTTGA_TAGs  > 29m5_heart_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/29m6_heart_GATCAG_TAGs  > 29m6_heart_RNA_repeats.txt 
#
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/3m1_Liver_ATCACG_TAGs   > 3m1_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/3m2_Liver_CGATGT_TAGs   > 3m2_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/3m3_Liver_TTAGGC_TAGs   > 3m3_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/12m2_Liver_ACAGTG_TAGs  > 12m2_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/12m3_Liver_GCCAAT_TAGs  > 12m3_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/12m1_Liver_TGACCA_TAGs  > 12m1_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/29m3_Liver_GATCAG_TAGs  > 29m3_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/29m2_Liver_ACTTGA_TAGs  > 29m2_Liver_RNA_repeats.txt 
#analyzeRepeats.pl repeats mm9 -raw -d $RNATAGS/29m1_Liver_CAGATC_TAGs  > 29m1_Liver_RNA_repeats.txt 
#

# over genes, to normalize libraries
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_3m1_TAGs          > Cereb_3m1_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_3m2_TAGs          > Cereb_3m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_3m3_TAGs          > Cereb_3m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_12m1_TAGs         > Cereb_12m1_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_12m2_TAGs         > Cereb_12m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_12m3_TAGs         > Cereb_12m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_29m2_TAGs         > Cereb_29m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_29m3_TAGs         > Cereb_29m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/Cereb_29m4_TAGs         > Cereb_29m4_RNA_genes.txt
#
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_3m1_ATCACG_TAGs      > OB_3m1_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_3m2_CGATGT_TAGs      > OB_3m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_3m3_TTAGGC_TAGs      > OB_3m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_12m1_TGACCA_TAGs     > OB_12m1_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_12m3_GCCAAT_TAGs     > OB_12m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_29m2_ACTTGA_TAGs     > OB_29m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_29m3_GATCAG_TAGs     > OB_29m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/OB_29m1_CAGATC_TAGs     > OB_29m1_RNA_genes.txt
#
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/NPC_3m5-CGATGT_TAGs     > NPC_3m5_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/NPC_3m6-TGACCA_TAGs     > NPC_3m6_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/NPC_12m5-ACAGTG_TAGs    > NPC_12m5_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/NPC_12m6-GCCAAT_TAGs    > NPC_12m6_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/NPC_29m6-CTTGTA_TAGs    > NPC_29m6_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/NPC_29m5-CAGATC_TAGs    > NPC_29m5_RNA_genes.txt
#
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/3m4_heart_ATCACG_TAGs   > 3m4_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/3m5_heart_CGATGT_TAGs   > 3m5_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/3m6_heart_TTAGGC_TAGs   > 3m6_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/12m4_heart_TGACCA_TAGs  > 12m4_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/12m5_heart_ACAGTG_TAGs  > 12m5_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/12m6_heart_GCCAAT_TAGs  > 12m6_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/29m4_heart_CAGATC_TAGs  > 29m4_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/29m5_heart_ACTTGA_TAGs  > 29m5_heart_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/29m6_heart_GATCAG_TAGs  > 29m6_heart_RNA_genes.txt
#
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/3m1_Liver_ATCACG_TAGs   > 3m1_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/3m2_Liver_CGATGT_TAGs   > 3m2_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/3m3_Liver_TTAGGC_TAGs   > 3m3_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/12m2_Liver_ACAGTG_TAGs  > 12m2_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/12m3_Liver_GCCAAT_TAGs  > 12m3_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/12m1_Liver_TGACCA_TAGs  > 12m1_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/29m3_Liver_GATCAG_TAGs  > 29m3_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/29m2_Liver_ACTTGA_TAGs  > 29m2_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d $RNATAGS/29m1_Liver_CAGATC_TAGs  > 29m1_Liver_RNA_genes.txt
#


# 2017-12-12
# now combine genes and repeats
#Parse_homer_counts.pl OB_29m3_RNA_repeats.txt OB_29m3_RNA_genes.txt
#Parse_homer_counts.pl OB_29m2_RNA_repeats.txt OB_29m2_RNA_genes.txt
#Parse_homer_counts.pl OB_29m1_RNA_repeats.txt OB_29m1_RNA_genes.txt
#Parse_homer_counts.pl OB_12m3_RNA_repeats.txt OB_12m3_RNA_genes.txt
#Parse_homer_counts.pl OB_12m1_RNA_repeats.txt OB_12m1_RNA_genes.txt
#Parse_homer_counts.pl OB_3m3_RNA_repeats.txt OB_3m3_RNA_genes.txt
#Parse_homer_counts.pl OB_3m2_RNA_repeats.txt OB_3m2_RNA_genes.txt
#Parse_homer_counts.pl OB_3m1_RNA_repeats.txt OB_3m1_RNA_genes.txt
#Parse_homer_counts.pl NPC_29m6_RNA_repeats.txt NPC_29m6_RNA_genes.txt
#Parse_homer_counts.pl NPC_29m5_RNA_repeats.txt NPC_29m5_RNA_genes.txt
#Parse_homer_counts.pl NPC_12m6_RNA_repeats.txt NPC_12m6_RNA_genes.txt
#Parse_homer_counts.pl NPC_12m5_RNA_repeats.txt NPC_12m5_RNA_genes.txt
#Parse_homer_counts.pl NPC_3m6_RNA_repeats.txt NPC_3m6_RNA_genes.txt
#Parse_homer_counts.pl NPC_3m5_RNA_repeats.txt NPC_3m5_RNA_genes.txt
#Parse_homer_counts.pl Cereb_29m4_RNA_repeats.txt Cereb_29m4_RNA_genes.txt
#Parse_homer_counts.pl Cereb_29m3_RNA_repeats.txt Cereb_29m3_RNA_genes.txt
#Parse_homer_counts.pl Cereb_29m2_RNA_repeats.txt Cereb_29m2_RNA_genes.txt
#Parse_homer_counts.pl Cereb_12m3_RNA_repeats.txt Cereb_12m3_RNA_genes.txt
#Parse_homer_counts.pl Cereb_12m2_RNA_repeats.txt Cereb_12m2_RNA_genes.txt
#Parse_homer_counts.pl Cereb_12m1_RNA_repeats.txt Cereb_12m1_RNA_genes.txt
#Parse_homer_counts.pl Cereb_3m3_RNA_repeats.txt Cereb_3m3_RNA_genes.txt
#Parse_homer_counts.pl Cereb_3m2_RNA_repeats.txt Cereb_3m2_RNA_genes.txt
#Parse_homer_counts.pl Cereb_3m1_RNA_repeats.txt Cereb_3m1_RNA_genes.txt
#Parse_homer_counts.pl 29m6_heart_RNA_repeats.txt 29m6_heart_RNA_genes.txt
#Parse_homer_counts.pl 29m5_heart_RNA_repeats.txt 29m5_heart_RNA_genes.txt
#Parse_homer_counts.pl 29m4_heart_RNA_repeats.txt 29m4_heart_RNA_genes.txt
#Parse_homer_counts.pl 29m3_Liver_RNA_repeats.txt 29m3_Liver_RNA_genes.txt
#Parse_homer_counts.pl 29m2_Liver_RNA_repeats.txt 29m2_Liver_RNA_genes.txt
#Parse_homer_counts.pl 29m1_Liver_RNA_repeats.txt 29m1_Liver_RNA_genes.txt
#Parse_homer_counts.pl 12m6_heart_RNA_repeats.txt 12m6_heart_RNA_genes.txt
#Parse_homer_counts.pl 12m5_heart_RNA_repeats.txt 12m5_heart_RNA_genes.txt
#Parse_homer_counts.pl 12m4_heart_RNA_repeats.txt 12m4_heart_RNA_genes.txt
#Parse_homer_counts.pl 12m3_Liver_RNA_repeats.txt 12m3_Liver_RNA_genes.txt
#Parse_homer_counts.pl 12m2_Liver_RNA_repeats.txt 12m2_Liver_RNA_genes.txt
#Parse_homer_counts.pl 12m1_Liver_RNA_repeats.txt 12m1_Liver_RNA_genes.txt
#Parse_homer_counts.pl 3m6_heart_RNA_repeats.txt 3m6_heart_RNA_genes.txt
#Parse_homer_counts.pl 3m5_heart_RNA_repeats.txt 3m5_heart_RNA_genes.txt
#Parse_homer_counts.pl 3m4_heart_RNA_repeats.txt 3m4_heart_RNA_genes.txt
#Parse_homer_counts.pl 3m3_Liver_RNA_repeats.txt 3m3_Liver_RNA_genes.txt
#Parse_homer_counts.pl 3m2_Liver_RNA_repeats.txt 3m2_Liver_RNA_genes.txt
#Parse_homer_counts.pl 3m1_Liver_RNA_repeats.txt 3m1_Liver_RNA_genes.txt

#parse_HOMER_mappings.pl 2017-12-12_Liver_RNA_genes_and_repeats_matrix.txt 3m1_Liver_gene_and_repeat_counts.txt 3m2_Liver_gene_and_repeat_counts.txt 3m3_Liver_gene_and_repeat_counts.txt 12m1_Liver_gene_and_repeat_counts.txt 12m2_Liver_gene_and_repeat_counts.txt 12m3_Liver_gene_and_repeat_counts.txt 29m1_Liver_gene_and_repeat_counts.txt 29m2_Liver_gene_and_repeat_counts.txt 29m3_Liver_gene_and_repeat_counts.txt
#parse_HOMER_mappings.pl 2017-12-12_heart_RNA_genes_and_repeats_matrix.txt 3m4_heart_gene_and_repeat_counts.txt 3m5_heart_gene_and_repeat_counts.txt 3m6_heart_gene_and_repeat_counts.txt 12m4_heart_gene_and_repeat_counts.txt 12m5_heart_gene_and_repeat_counts.txt 12m6_heart_gene_and_repeat_counts.txt 29m4_heart_gene_and_repeat_counts.txt 29m5_heart_gene_and_repeat_counts.txt 29m6_heart_gene_and_repeat_counts.txt
#parse_HOMER_mappings.pl 2017-12-12_Cereb_RNA_genes_and_repeats_matrix.txt Cereb_3m1_gene_and_repeat_counts.txt Cereb_3m2_gene_and_repeat_counts.txt Cereb_3m3_gene_and_repeat_counts.txt Cereb_12m1_gene_and_repeat_counts.txt Cereb_12m2_gene_and_repeat_counts.txt Cereb_12m3_gene_and_repeat_counts.txt Cereb_29m2_gene_and_repeat_counts.txt Cereb_29m3_gene_and_repeat_counts.txt Cereb_29m4_gene_and_repeat_counts.txt
#parse_HOMER_mappings.pl 2017-12-12_NPC_RNA_genes_and_repeats_matrix.txt NPC_3m5_gene_and_repeat_counts.txt NPC_3m6_gene_and_repeat_counts.txt NPC_12m5_gene_and_repeat_counts.txt NPC_12m6_gene_and_repeat_counts.txt NPC_29m5_gene_and_repeat_counts.txt NPC_29m6_gene_and_repeat_counts.txt
#parse_HOMER_mappings.pl 2017-12-12_OB_RNA_genes_and_repeats_matrix.txt OB_3m1_gene_and_repeat_counts.txt OB_3m2_gene_and_repeat_counts.txt OB_3m3_gene_and_repeat_counts.txt OB_12m1_gene_and_repeat_counts.txt OB_12m3_gene_and_repeat_counts.txt OB_29m1_gene_and_repeat_counts.txt OB_29m2_gene_and_repeat_counts.txt OB_29m3_gene_and_repeat_counts.txt


# 2017-12-14
# analyze fibroblast data
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_3m1_TAGs          > Cereb_3m1_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_3m2_TAGs          > Cereb_3m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_3m3_TAGs          > Cereb_3m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_12m1_TAGs         > Cereb_12m1_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_12m2_TAGs         > Cereb_12m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_12m3_TAGs         > Cereb_12m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_29m2_TAGs         > Cereb_29m2_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_29m3_TAGs         > Cereb_29m3_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d Fibro/Cereb_29m4_TAGs         > Cereb_29m4_RNA_genes.txt

               
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m1_cohort1_dermal_fibro_TAGs          > 3m1_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m2_cohort1_dermal_fibro_TAGs          > 3m2_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m3_cohort4_dermal_fibro_TAGs          > 3m3_cohort4_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m3c10R_R_TAGs                         > 3m3_cohort10R_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m4_cohort1_dermal_fibro_TAGs          > 3m4_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m5c11L_R_TAGs                         > 3m5_cohort11L_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m6_cohort4_dermal_fibro_TAGs          > 3m6_cohort4_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/3m7_cohort4_dermal_fibro_TAGs          > 3m7_cohort4_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/12m1_cohort1_dermal_fibro_TAGs         > 12m1_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/12m2_cohort1_dermal_fibro_TAGs         > 12m2_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/12m4_cohort1_dermal_fibro_TAGs         > 12m4_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29M1C10R_GOOD_R_TAGs                   > 29m1_cohort10R_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29m2_cohort1_dermal_fibro_TAGs         > 29m2_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29m2_cohort4_dermal_fibro_TAGs         > 29m2_cohort4_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29m3_cohort1_dermal_fibro_TAGs         > 29m3_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29M3C10L_GOOD_R_TAGs                   > 29m3_cohort10L_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29m4_cohort1_dermal_fibro_TAGs         > 29m4_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29m4_cohort4_dermal_fibro_TAGs         > 29m4_cohort4_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29M4C11R_GOOD_R_TAGs                   > 29m4_cohort11R_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29M5C10L_BAD_R_TAGs                    > 29m5_cohort10L_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29M5C11R_BAD_R_TAGs                    > 29m5_cohort11R_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29m6_cohort1_dermal_fibro_TAGs         > 29m6_cohort1_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29M6C11L_BAD_R_TAGs                    > 29m6_cohort11L_dermal_fibro_RNA_repeats.txt    
#analyzeRepeats.pl repeats mm9 -raw -d ./Fibro/29m7_cohort1_dermal_fibro_TAGs         > 29m7_cohort1_dermal_fibro_RNA_repeats.txt    
#
#
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m1_cohort1_dermal_fibro_TAGs          > 3m1_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m2_cohort1_dermal_fibro_TAGs          > 3m2_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m3_cohort4_dermal_fibro_TAGs          > 3m3_cohort4_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m3c10R_R_TAGs                         > 3m3_cohort10R_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m4_cohort1_dermal_fibro_TAGs          > 3m4_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m5c11L_R_TAGs                         > 3m5_cohort11L_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m6_cohort4_dermal_fibro_TAGs          > 3m6_cohort4_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/3m7_cohort4_dermal_fibro_TAGs          > 3m7_cohort4_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/12m1_cohort1_dermal_fibro_TAGs         > 12m1_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/12m2_cohort1_dermal_fibro_TAGs         > 12m2_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/12m4_cohort1_dermal_fibro_TAGs         > 12m4_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29M1C10R_GOOD_R_TAGs                   > 29m1_cohort10R_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29m2_cohort1_dermal_fibro_TAGs         > 29m2_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29m2_cohort4_dermal_fibro_TAGs         > 29m2_cohort4_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29m3_cohort1_dermal_fibro_TAGs         > 29m3_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29M3C10L_GOOD_R_TAGs                   > 29m3_cohort10L_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29m4_cohort1_dermal_fibro_TAGs         > 29m4_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29m4_cohort4_dermal_fibro_TAGs         > 29m4_cohort4_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29M4C11R_GOOD_R_TAGs                   > 29m4_cohort11R_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29M5C10L_BAD_R_TAGs                    > 29m5_cohort10L_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29M5C11R_BAD_R_TAGs                    > 29m5_cohort11R_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29m6_cohort1_dermal_fibro_TAGs         > 29m6_cohort1_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29M6C11L_BAD_R_TAGs                    > 29m6_cohort11L_dermal_fibro_RNA_genes.txt    
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./Fibro/29m7_cohort1_dermal_fibro_TAGs         > 29m7_cohort1_dermal_fibro_RNA_genes.txt    


#Parse_homer_counts.pl 3m1_cohort1_dermal_fibro_RNA_repeats.txt 3m1_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 3m2_cohort1_dermal_fibro_RNA_repeats.txt 3m2_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 3m3_cohort4_dermal_fibro_RNA_repeats.txt 3m3_cohort4_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 3m3_cohort10R_dermal_fibro_RNA_repeats.txt 3m3_cohort10R_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 3m4_cohort1_dermal_fibro_RNA_repeats.txt 3m4_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 3m5_cohort11L_dermal_fibro_RNA_repeats.txt 3m5_cohort11L_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 3m6_cohort4_dermal_fibro_RNA_repeats.txt 3m6_cohort4_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 3m7_cohort4_dermal_fibro_RNA_repeats.txt 3m7_cohort4_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 12m1_cohort1_dermal_fibro_RNA_repeats.txt 12m1_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 12m2_cohort1_dermal_fibro_RNA_repeats.txt 12m2_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 12m4_cohort1_dermal_fibro_RNA_repeats.txt 12m4_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m1_cohort10R_dermal_fibro_RNA_repeats.txt 29m1_cohort10R_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m2_cohort1_dermal_fibro_RNA_repeats.txt 29m2_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m2_cohort4_dermal_fibro_RNA_repeats.txt 29m2_cohort4_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m3_cohort1_dermal_fibro_RNA_repeats.txt 29m3_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m3_cohort10L_dermal_fibro_RNA_repeats.txt 29m3_cohort10L_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m4_cohort1_dermal_fibro_RNA_repeats.txt 29m4_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m4_cohort4_dermal_fibro_RNA_repeats.txt 29m4_cohort4_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m4_cohort11R_dermal_fibro_RNA_repeats.txt 29m4_cohort11R_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m5_cohort10L_dermal_fibro_RNA_repeats.txt 29m5_cohort10L_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m5_cohort11R_dermal_fibro_RNA_repeats.txt 29m5_cohort11R_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m6_cohort1_dermal_fibro_RNA_repeats.txt 29m6_cohort1_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m6_cohort11L_dermal_fibro_RNA_repeats.txt 29m6_cohort11L_dermal_fibro_RNA_genes.txt    
#Parse_homer_counts.pl 29m7_cohort1_dermal_fibro_RNA_repeats.txt 29m7_cohort1_dermal_fibro_RNA_genes.txt    
#
#parse_HOMER_mappings.pl 2017-12-14_Dermal_Fibro_RNA_genes_and_repeats_matrix.txt 3m1_cohort1_dermal_fibro_gene_and_repeat_counts.txt 3m2_cohort1_dermal_fibro_gene_and_repeat_counts.txt 3m3_cohort4_dermal_fibro_gene_and_repeat_counts.txt 3m3_cohort10R_dermal_fibro_gene_and_repeat_counts.txt 3m4_cohort1_dermal_fibro_gene_and_repeat_counts.txt 3m5_cohort11L_dermal_fibro_gene_and_repeat_counts.txt 3m6_cohort4_dermal_fibro_gene_and_repeat_counts.txt 3m7_cohort4_dermal_fibro_gene_and_repeat_counts.txt 12m1_cohort1_dermal_fibro_gene_and_repeat_counts.txt 12m2_cohort1_dermal_fibro_gene_and_repeat_counts.txt 12m4_cohort1_dermal_fibro_gene_and_repeat_counts.txt 29m1_cohort10R_dermal_fibro_gene_and_repeat_counts.txt 29m2_cohort1_dermal_fibro_gene_and_repeat_counts.txt 29m2_cohort4_dermal_fibro_gene_and_repeat_counts.txt 29m3_cohort1_dermal_fibro_gene_and_repeat_counts.txt 29m3_cohort10L_dermal_fibro_gene_and_repeat_counts.txt 29m4_cohort1_dermal_fibro_gene_and_repeat_counts.txt 29m4_cohort4_dermal_fibro_gene_and_repeat_counts.txt 29m4_cohort11R_dermal_fibro_gene_and_repeat_counts.txt 29m5_cohort10L_dermal_fibro_gene_and_repeat_counts.txt 29m5_cohort11R_dermal_fibro_gene_and_repeat_counts.txt 29m6_cohort1_dermal_fibro_gene_and_repeat_counts.txt 29m6_cohort11L_dermal_fibro_gene_and_repeat_counts.txt 29m7_cohort1_dermal_fibro_gene_and_repeat_counts.txt 


# 2017-12-15
# analyze CR RNAseq data

#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m1_Cerebellum_TAGs     > 3m1_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m6_Cerebellum_TAGs     > 3m6_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m7_Cerebellum_TAGs     > 3m7_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m10_Cereb_TAGs         > 3m10_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m1_Cerebellum_TAGs    > 24m1_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m6_Cerebellum_TAGs    > 24m6_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m7_Cerebellum_TAGs    > 24m7_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m8_Cereb_TAGs         > 24m8_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_1_Cerebellum_TAGs    > CR_1_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_2_Cerebellum_TAGs    > CR_2_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_3_Cerebellum_TAGs    > CR_3_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_4_Cerebellum_TAGs    > CR_4_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_5_Cerebellum_TAGs    > CR_5_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_7_Cerebellum_TAGs    > CR_7_Cereb_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_8_Cereb_TAGs         > CR_8_Cereb_RNA_repeats.txt
#
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m1_Liver_TAGs          > 3m1_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m6_Liver_TAGs          > 3m6_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m7_Liver_TAGs          > 3m7_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/3m10_Liver_TAGs         > 3m10_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m1_Liver_TAGs         > 24m1_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m6_Liver_TAGs         > 24m6_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m7_Liver_TAGs         > 24m7_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/24m8_Liver_TAGs         > 24m8_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_1_Liver_TAGs         > CR_1_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_2_Liver_TAGs         > CR_2_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_3_Liver_TAGs         > CR_3_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_4_Liver_TAGs         > CR_4_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_5_Liver_TAGs         > CR_5_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_6_Liver_TAGs         > CR_6_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_7_Liver_TAGs         > CR_7_Liver_RNA_repeats.txt
#analyzeRepeats.pl repeats mm9 -raw -d ./CR/CR_8_Liver_TAGs         > CR_8_Liver_RNA_repeats.txt
#
#
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m1_Cerebellum_TAGs     > 3m1_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m6_Cerebellum_TAGs     > 3m6_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m7_Cerebellum_TAGs     > 3m7_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m10_Cereb_TAGs         > 3m10_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m1_Cerebellum_TAGs    > 24m1_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m6_Cerebellum_TAGs    > 24m6_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m7_Cerebellum_TAGs    > 24m7_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m8_Cereb_TAGs         > 24m8_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_1_Cerebellum_TAGs    > CR_1_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_2_Cerebellum_TAGs    > CR_2_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_3_Cerebellum_TAGs    > CR_3_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_4_Cerebellum_TAGs    > CR_4_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_5_Cerebellum_TAGs    > CR_5_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_7_Cerebellum_TAGs    > CR_7_Cereb_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_8_Cereb_TAGs         > CR_8_Cereb_RNA_genes.txt
#
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m1_Liver_TAGs          > 3m1_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m6_Liver_TAGs          > 3m6_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m7_Liver_TAGs          > 3m7_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/3m10_Liver_TAGs         > 3m10_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m1_Liver_TAGs         > 24m1_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m6_Liver_TAGs         > 24m6_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m7_Liver_TAGs         > 24m7_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/24m8_Liver_TAGs         > 24m8_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_1_Liver_TAGs         > CR_1_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_2_Liver_TAGs         > CR_2_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_3_Liver_TAGs         > CR_3_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_4_Liver_TAGs         > CR_4_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_5_Liver_TAGs         > CR_5_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_6_Liver_TAGs         > CR_6_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_7_Liver_TAGs         > CR_7_Liver_RNA_genes.txt
#analyzeRepeats.pl rna mm9 -raw -condenseGenes -d ./CR/CR_8_Liver_TAGs         > CR_8_Liver_RNA_genes.txt
#
#Parse_homer_counts.pl 3m1_Cereb_RNA_repeats.txt 3m1_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 3m6_Cereb_RNA_repeats.txt 3m6_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 3m7_Cereb_RNA_repeats.txt 3m7_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 3m10_Cereb_RNA_repeats.txt 3m10_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 24m1_Cereb_RNA_repeats.txt 24m1_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 24m6_Cereb_RNA_repeats.txt 24m6_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 24m7_Cereb_RNA_repeats.txt 24m7_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 24m8_Cereb_RNA_repeats.txt 24m8_Cereb_RNA_genes.txt
#Parse_homer_counts.pl CR_1_Cereb_RNA_repeats.txt CR_1_Cereb_RNA_genes.txt
#Parse_homer_counts.pl CR_2_Cereb_RNA_repeats.txt CR_2_Cereb_RNA_genes.txt
#Parse_homer_counts.pl CR_3_Cereb_RNA_repeats.txt CR_3_Cereb_RNA_genes.txt
#Parse_homer_counts.pl CR_4_Cereb_RNA_repeats.txt CR_4_Cereb_RNA_genes.txt
#Parse_homer_counts.pl CR_5_Cereb_RNA_repeats.txt CR_5_Cereb_RNA_genes.txt
#Parse_homer_counts.pl CR_7_Cereb_RNA_repeats.txt CR_7_Cereb_RNA_genes.txt
#Parse_homer_counts.pl CR_8_Cereb_RNA_repeats.txt CR_8_Cereb_RNA_genes.txt
#Parse_homer_counts.pl 3m1_Liver_RNA_repeats.txt 3m1_Liver_RNA_genes.txt
#Parse_homer_counts.pl 3m6_Liver_RNA_repeats.txt 3m6_Liver_RNA_genes.txt
#Parse_homer_counts.pl 3m7_Liver_RNA_repeats.txt 3m7_Liver_RNA_genes.txt
#Parse_homer_counts.pl 3m10_Liver_RNA_repeats.txt 3m10_Liver_RNA_genes.txt
#Parse_homer_counts.pl 24m1_Liver_RNA_repeats.txt 24m1_Liver_RNA_genes.txt
#Parse_homer_counts.pl 24m6_Liver_RNA_repeats.txt 24m6_Liver_RNA_genes.txt
#Parse_homer_counts.pl 24m7_Liver_RNA_repeats.txt 24m7_Liver_RNA_genes.txt
#Parse_homer_counts.pl 24m8_Liver_RNA_repeats.txt 24m8_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_1_Liver_RNA_repeats.txt CR_1_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_2_Liver_RNA_repeats.txt CR_2_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_3_Liver_RNA_repeats.txt CR_3_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_4_Liver_RNA_repeats.txt CR_4_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_5_Liver_RNA_repeats.txt CR_5_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_6_Liver_RNA_repeats.txt CR_6_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_7_Liver_RNA_repeats.txt CR_7_Liver_RNA_genes.txt
#Parse_homer_counts.pl CR_8_Liver_RNA_repeats.txt CR_8_Liver_RNA_genes.txt


parse_HOMER_mappings.pl 2017-12-15_Cerebellum_CR_data_RNA_genes_and_repeats_matrix.txt 3m1_Cereb_gene_and_repeat_counts.txt 3m6_Cereb_gene_and_repeat_counts.txt 3m7_Cereb_gene_and_repeat_counts.txt 3m10_Cereb_gene_and_repeat_counts.txt 24m1_Cereb_gene_and_repeat_counts.txt 24m6_Cereb_gene_and_repeat_counts.txt 24m7_Cereb_gene_and_repeat_counts.txt 24m8_Cereb_gene_and_repeat_counts.txt CR_1_Cereb_gene_and_repeat_counts.txt CR_2_Cereb_gene_and_repeat_counts.txt CR_3_Cereb_gene_and_repeat_counts.txt CR_4_Cereb_gene_and_repeat_counts.txt CR_5_Cereb_gene_and_repeat_counts.txt CR_7_Cereb_gene_and_repeat_counts.txt CR_8_Cereb_gene_and_repeat_counts.txt
parse_HOMER_mappings.pl 2017-12-15_Liver_CR_data_RNA_genes_and_repeats_matrix.txt 3m1_Liver_gene_and_repeat_counts.txt 3m6_Liver_gene_and_repeat_counts.txt 3m7_Liver_gene_and_repeat_counts.txt 3m10_Liver_gene_and_repeat_counts.txt 24m1_Liver_gene_and_repeat_counts.txt 24m6_Liver_gene_and_repeat_counts.txt 24m7_Liver_gene_and_repeat_counts.txt 24m8_Liver_gene_and_repeat_counts.txt CR_1_Liver_gene_and_repeat_counts.txt CR_2_Liver_gene_and_repeat_counts.txt CR_3_Liver_gene_and_repeat_counts.txt CR_4_Liver_gene_and_repeat_counts.txt CR_5_Liver_gene_and_repeat_counts.txt CR_6_Liver_gene_and_repeat_counts.txt CR_7_Liver_gene_and_repeat_counts.txt CR_8_Liver_gene_and_repeat_counts.txt 





























