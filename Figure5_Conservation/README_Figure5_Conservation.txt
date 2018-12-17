README - Figure6_Conservation
###########################

**************************
*** R Package versions ***
DESeq2 v1.6.3 (Most Data)
DESeq2 v1.16.1 (Human Both sex Data)
mHG v1.1
limma’ v 3.32.10
**************************
***   Other software   ***
Kallisto v0.43.0
homer	v4.10.3
**************************

* Comparisons
	- Rat
		process_rat_RNAseq.R: script to get DE genes (Kallisto tables are in "kallisto_output")
		RNAseq_process_functions.R: functions called by above script
		COmpare_to_Mouse_data.R: compare to the mouse RNAseq results
		cross_species_comparison_FUN.R: functions called by above script
		
	- GTex
		parse_GTex.R: parse meta data from GTEx
		process_GTex_v3.R: perform DE analysis for GTex data
		GTEX_analysis_functions_v3bis.R: functions called by above script
		Compare_to_Mouse_data_V2.R: compare to the mouse RNAseq results
		cross_species_comparison_FUN.R: functions called by above script

	- Human_microarray
		Process_microarray.R: process GEO array (inputs from GSE61260 files are provided for convenience)
		Compare_horvath_to_mouseaging.R: Microarray data analysis and comparison to mouse

	- Killifish_aging_RNAseq
		Brain
			parse_kallisto_mappings.pl: combine individual samples into matrix for R analysis (example input provided; ran through combine_kallisto_calls.sh)
			Process_RNAseq_forDEseq_brain.R: perform DE analysis for Brain data
			
		Liver
			parse_kallisto_mappings.pl: combine individual samples into matrix for R analysis (same code/procedure as above; ran through combine_kallisto_calls.sh)
			Process_RNAseq_forDEseq.R: perform DE analysis for Liver data
			
* Pathway_enrichment
	- run_GOrilla_like_for_nonMouse_v2.R: run pathway enrichment analysis (script results in "FDR5percent" folder). Data input from Comparisons scripts; uses orthology tables in Orthology.
	- GOrilla_statistics_functions_for_nonMouse.R: functions called by above script 
	
	- Process_gorilla_like_results_v3.R: get result summary (needs "FDR5percent" and Mouse_Data folders for input; Stats_tables and MYDATA need to be created for the script to run)
	- Process_gorilla_like_results_FUNCTIONS_v3.R: functions called by above script
