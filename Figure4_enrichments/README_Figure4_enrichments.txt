README - Figure4_enrichments
###########################

**************************
*** R Package versions ***
DESeq2 v1.6.3 (Normal Mouse Data)
DESeq2 v1.16.1 (TE Data)
Diffbind v2.4.8 (TE Data)
mHG v1.1
**************************
***   Other software   ***
homer	v4.10.3
TE-transcript 1.5.1
**************************

* Venn_different_molecules
	- Get_Venn_per_molecules.R: script to make Venn diagrams of overlaps (see README_Venn.txt)
	
* Functional_enrichment_RNA
	- run_GOrilla_like.R: enrichment script (pathway files are in "PATHWAYS")
	- GOrilla_statistics_functions.R: functions called by above script
	
	- Process_gorilla_like_resultsREV.R: summarize GOrilla results (use Results_of_1st_script content to run it)
	- Process_gorilla_like_results_FUNCTIONS_REV.R: functions called by above script
	
	
* Functional_enrichment_ChIP (uses a slightly modified input as above, so example input is not included)
	- run_GOrilla_like.R: enrichment script (pathway files are in "PATHWAYS")
	- GOrilla_statistics_functions.R:functions called by above script
	- Process_gorilla_like_results.R: summarize GOrilla results
	- Process_gorilla_like_results_FUNCTIONS.R: functions called by above script
	
* TE_RNA_analysis
	- TE-transcript
		tetrasncript_cluster_2017-03-13.sh: run TE-transcript (output count tables in TE-transcript_output)
		parse_TE_families.pl: parse TE family names from output (result: 2017-09-12_Parsed_TE_names_list.txt)
		process_RNAseq_datasets_TE_v2.R: run DE analysis for TEs (inputs in TE-transcript_output)
		TE_RNAseq_analysis_functions.R: functions called by above script

	- HOMER (different software to run the same analysis)
		quantify_repeat_expression.sh: bash script of pipeline
		Parse_homer_counts.pl: combine counts for genes and repeats so that DEseq2 can build the correct mean/variance model
		parse_HOMER_mappings.pl: combine different samples into a matrix for processing
		
		Process_Repeats_Exp_forAnne.R: DE analysis in R
		TE_DE_functions_forAnne.R: functions called by above script

* TE_Chromatin_analysis
	- Pipeline is fully described in README.txt


* IPA
	- Parse_IPA_v2.R: parse IPA results from Salah
	- Parse_IPA_functions_v2.R: functions called by above script
