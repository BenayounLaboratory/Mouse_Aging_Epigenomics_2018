# setwd('/Volumes/BB_Backup_3/BD_aging_project/2018-09_revision_analyses/Other_aging_transcriptomes/From_Raw/Pathway_analysis/')
setwd('/Users/BB_2012/Dropbox/manuscripts_and_publications/2018_aging_epigenomics_data_description/aging_omics_paper/Genome_Research_submission/New_code_for_checking/Public_transcriptome_analysis/Pathway_analysis/Running_test/')
library('mHG')
options(stringsAsFactors=F)

source('GOrilla_statistics_functions.R')

# 2018-09-24
# Run on other datasets for revision
# Create folders ALL_PATHWAYS, FDR5percent, RData before running

my.gmt.sets <- c(paste("/Users/BB_2012/Dropbox/manuscripts_and_publications/2018_aging_epigenomics_data_description/aging_omics_paper/Genome_Research_submission/New_code_for_checking/Public_transcriptome_analysis/Pathway_analysis/Running_test/INPUT/GMT/Param/", list.files("/Users/BB_2012/Dropbox/manuscripts_and_publications/2018_aging_epigenomics_data_description/aging_omics_paper/Genome_Research_submission/New_code_for_checking/Public_transcriptome_analysis/Pathway_analysis/Running_test/INPUT/GMT/Param",pattern = "\\.gmt$"), sep="/"),
                 paste("/Users/BB_2012/Dropbox/manuscripts_and_publications/2018_aging_epigenomics_data_description/aging_omics_paper/Genome_Research_submission/New_code_for_checking/Public_transcriptome_analysis/Pathway_analysis/Running_test/INPUT/GMT/MSigDB", list.files("//Users/BB_2012/Dropbox/manuscripts_and_publications/2018_aging_epigenomics_data_description/aging_omics_paper/Genome_Research_submission/New_code_for_checking/Public_transcriptome_analysis/Pathway_analysis/Running_test/INPUT/GMT/MSigDB",pattern = "\\.gmt$"), sep="/"))

my.gmt.set.names <- c("KEGG_2017_no_diseases_UC",
                      "KEGG_2017_no_diseases",
                      "KEGG_2017_All",
                      
                      "C2_CGP",
                      "Biocarta",
                      "Kegg",
                      "reactome",
                      "C2cp",
                      "C3Mir",
                      "C3TF",
                      "C5BP",
                      "C5CC",
                      "C5MF",
                      "C7_All",
                      "MSigDB_Hallmarks"
)

cbind(my.gmt.sets, my.gmt.set.names)

# catalog pathways to run
my.run <- c(1,15)

source('GOrilla_statistics_functions.R')

##################################################################################### 
# Load RData RNA with age results
load("./INPUT/RData/2018-09-24_Boisvert_Cereb_astrocytes_RNAseq.RData")
load("./INPUT/RData/2018-09-24_Bochkis_Liver_RNAseq.RData")
load("./INPUT/RData/2018-09-24_White_Liver_RNAseq.RData")

# My data
run_pathway_enrich("Bochkis_Liver", my.bochkis.RNAseq.process[[1]])
run_pathway_enrich("White_Liver", my.white.RNAseq.process[[1]])
run_pathway_enrich("Boisvert_cereb_astrocytes", my.boisvert.RNAseq.process[[1]])

