# setwd('/Volumes/BB_Backup_3/BD_aging_project/2018-09_revision_analyses/Other_aging_transcriptomes/From_Raw/Pathway_analysis/SUMMARY')
setwd('/Users/BB_2012/Dropbox/manuscripts_and_publications/2018_aging_epigenomics_data_description/aging_omics_paper/Genome_Research_submission/New_code_for_checking/Public_transcriptome_analysis/Pathway_analysis/summarizing_test//')

options(stringsAsFactors=F)

source('Process_gorilla_like_results_FUNCTIONS_vCOMP_v2.R')

# 2018-09-24
# compare our datasets to public
# needs "FDR5" folder in ../ folder to run (here, in INPUT)
# Stats table from previous run are in INPUT


############################################################################################################################################################################################################################
# read in stat table from mouse run
my.hallmark.mouse <- read.table('INPUT/2018-01-05_Enrichment_table_MSIgDB_Hallmark_Datasets_pathways_significant_in_4_or_more.txt',sep = "\t", header=T)
my.kegg.mouse2 <- read.csv('INPUT/2017-06-17_Enrichment_table_KEGG_2017_pathways_significant_in_4_or_more.txt',sep = "\t", header=T)


get_enrich_balloons_all_species("Hallmark", my.hallmark.mouse)
get_enrich_balloons_all_species("KEGG_2017", my.kegg.mouse2)
