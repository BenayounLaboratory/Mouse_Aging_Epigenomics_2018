# setwd('/Volumes/BB_Backup_3/BD_aging_project/Public_datasets/RNA_profile_other_species/Pathway_enrichment/Summary')
setwd('~/Downloads/Pathway_enrichment/')
options(stringsAsFactors=F)

source('Process_gorilla_like_results_FUNCTIONS_v3.R')

# need to create output folders MYDATA and Stats_tables

############################################################################################################################################################################################################################
# read in stat table from mouse run
my.hallmark.mouse <- read.table('./Mouse_Data/Mouse_tables/2017-05-22_Enrichment_table_MSIgDB_Hallmark_Datasets_pathways_significant_in_4_or_more.txt',sep = "\t", header=T)
my.kegg.mouse2 <- read.csv('./Mouse_Data/Mouse_tables/2017-06-17_Enrichment_table_KEGG_2017_pathways_significant_in_4_or_more.txt',sep = "\t", header=T)

get_enrich_balloons_all_species("MSIgDB_Hallmark_Datasets", my.hallmark.mouse)
get_enrich_balloons_all_species("KEGG_2017", my.kegg.mouse2)


