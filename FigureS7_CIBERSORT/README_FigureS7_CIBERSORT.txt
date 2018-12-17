README - Figure5_CIBERSORT
###########################

* RNAseq_datasets_for_Deconvolution:
	- 2017-01-18
		process_read_counts_v10.R: process all public RNAseq datasets to generate reference RNAseq sets
			Feature count outputs for the >1000 samples are in "Feature_count_output"
	
			CIBERSORT was run on the webserver https://cibersort.stanford.edu/ with these files as inputs	
	
	- InSilicoMixtures
		Create_test_mixtures_v2.R: create in silico mixtures for validation of signature matrix
		
	- CiBERSORT_plotting
		get_CIBERSORT_results_v2.R: create beehive plots for figure 5
		(Input are in Figure5_CIBERSORT/RNAseq_datasets_for_Deconvolution/CIBERSORT)
		
* VST_transformed_coutns:
	- process_RNAseq_datasets.R: run VST transformation of DEseq2 counts prior to CIBERSORT treatment
	- RNAseq_analysis_functions.R: functions called by above script

* Disease_models:
	Alzheimer's models RNA-seq CIBERSORT results
	
* Immune_Gene_Expression
	Manual_Test.R: manual extraction of expression data for boxplots of expression
	2016-02-05_ALL_global_variance estimate_DESeq2_LINEAR_model_with_age _log2_counts_matrix.txt: input normalized count matrix