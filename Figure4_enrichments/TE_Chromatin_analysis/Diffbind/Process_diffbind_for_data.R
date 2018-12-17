setwd('/Volumes/BB_Backup_3/BD_aging_project/2018-09_revision_analyses/TE_chromatin/RepeatMasker/Diffbind/')

source('diffbind_process_functions.R')

# 2018-09-12
# run diffbind to extract normalized intensity of ChIPs near repeats

###############################################################################
####################   Diffbind Analysis : H3K4me3 height  ####################
################################################################################
run_diffbind("Heart", "H3K4me3", "./Input/Heart_FIXSEQ_exp_K4me3_Repeats.csv")
run_diffbind("Liver", "H3K4me3", "./Input/Liver_FIXSEQ_exp_K4me3_Repeats.csv")
run_diffbind("Cerebellum", "H3K4me3", "./Input/Cereb_FIXSEQ_exp_K4me3_Repeats.csv")
run_diffbind("OB", "H3K4me3", "./Input/OB_FIXSEQ_exp_K4me3_Repeats.csv")
run_diffbind("NSPC", "H3K4me3", "./Input/NPC_FIXSEQ_exp_K4me3_Repeats.csv")

###############################################################################
####################   Diffbind Analysis : H3K27ac height  ####################
################################################################################
run_diffbind("Heart", "H3K27ac", "./Input/Heart_FIXSEQ_exp_K27ac_Repeats.csv")
run_diffbind("Liver", "H3K27ac", "./Input/Liver_FIXSEQ_exp_K27ac_Repeats.csv")
run_diffbind("Cerebellum", "H3K27ac", "./Input/Cereb_FIXSEQ_exp_K27ac_Repeats.csv")
run_diffbind("OB", "H3K27ac", "./Input/OB_FIXSEQ_exp_K27ac_Repeats.csv")
run_diffbind("NSPC", "H3K27ac", "./Input/NPC_FIXSEQ_exp_K27ac_Repeats.csv")
