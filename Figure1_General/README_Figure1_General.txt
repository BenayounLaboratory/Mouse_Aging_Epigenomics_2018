README - Figure1_General
###########################

**************************
*** R Package versions ***
DESeq2 v1.6.3
Diffbind v1.12.3
**************************
***   Other software   ***
BedTools-2.16.0/2.16.1
MACS v2.08
Trimgalore v0.3.1 
STAR v2.4.0j
subread v1.4.5-p1
FIXSeq fixseq-e2cc14f19764
homer	v4.10.3
**************************

* Generic_NGS_scripts
	- process_FQ.sh: Clean up raw ChIP fastq files for downstream analysis
	- map_and_BAM.sh: map ChIP-seq data to genomes
	- clean_BAM_FS.sh: Run Fixseq on BAM aligned file to remove duplicated reads
	
* ChIP-seq_processing
	This code is separated in 2 sections: 
		(i) pre-process each tissue separately (Heart_example is given)
		(ii) process all datasets within a uniform pipeline (All_tissues)
		
	- Heart_example
		- H3K4me3_breadth
			merge_bed_call_calibration_heart.sh: 
				(i) get all merged reads for meta-peak calling
				(ii) downsample each H3 ChIP to the lowest read number across H3 ChIPs
				(iii) get coverage histograms of meta-peaks across downsamplings of H3K4me3 ChIPs by 5% increments  
				(iv) find correct adjustment level (see get_cov_cdf_heart.R code)
				(v) call individual peaks with callibrated reads
				(vi) create a file with all breadth coverage of each sample across metapeaks
				(vii) create bedgraphs and beds for visualization
				
			get_cov_cdf_heart.R (needs files in ./hists/ to run program): 
				uses the "Breadth_compare_bundle" code to determine appropriate downsampling levels
				
		- H3K4me3_height_heart
			K4me3_height_Diffbind.R: uses diffbind to get normalized read densities in each sample over metapeaks
			(Input are mapped reads - files are too big for GitHub.)
			
		- H3K27ac_height_heart
			H3K27ac_typical_enhancer_heaight_changes_Heart.R: uses diffbind to get normalized read densities in each sample over metapeaks
			(Input are mapped reads - files are too big for GitHub.)

		- SuperEnhancer
			get_Enhancers_heart.sh:
				(i) get merged H3K27ac reads for meta-peak calling
				(ii) get stiched peaks within 12.5kb as in Super Enhancer paper  
				(iii) get coverage in reads over stitched enhancers for each ChIP sampe
				(vi) create a file with all coverage of each sample across metapeaks/enhancers

			SuperEnhancer-Lite_analysis_UPDATED.R: call super enhancers across samples, and run diffbind to get normalized read densities in each sample over meta super enhancers

	- All_tissues
		- Breadth_modeling
			Merge_H3K4me3_breadth_files.pl: script called to combine breadth calling over tissues for general PCA/MDS analysis (input files in Breadth_Files; see README for function call)
			Get_global_clustering_v2.R: general MDS/PCA analysis using output of above
			
			process_breadth_datasets_linear_model_vPCA.R: process differential breadth over all tissues with uniform pipeline
			breadth_analysis_functions_linear_v4PCA.R: functions called by above script
			
		- H3K4me3_height
			Merge_H3K4me3_height_files.pl: combine height calling over tissues for general PCA/MDS analysis (see README for input command)
			Get_global_clustering.R: general MDS/PCA analysis
			
			Process_height_H3K4me3_linear_model_cor_vPCA.R: process differential H3K4me3 height over all tissues with uniform pipeline
			height_analysis_functions_linear_v3PCA.R : functions called by above script	
			
		- H3K27ac_height (the scripts are modified from H3K4me3_height to use the H3K27ac input, so example files are not included)
			Get_global_clustering.R: general MDS/PCA analysis (merged file generated with Merge_H3K4me3_height_files.pl script, as input files are in the same format, as are desired outputs)

			Process_height_H3K27ac_linear_modeling_with_reseq_cor.R: process differential H3K27ac height over all tissues with uniform pipeline
			H3K27ac_height_analysis_functions_linear_v2.R: functions called by above script
			
		- SuperEnhancer (the scripts are modified from H3K4me3_height to use the SE input, so example files are not included)
			Merge_SE_files.pl: combine Super Enahncer intensity calling over tissues for general PCA/MDS analysis (see README for function call)
			Get_global_clustering.R: general MDS/PCA analysis
			Process_height_in SE_H3K27ac_linear_modeling_with_reseq_PCA.R: process differential H3K27ac height at super enahncers over all tissues with uniform pipeline
			SE_H3K27ac_height_analysis_functions_linear_vPCA.R: functions called by above script


* Breadth_compare_bundle: collection of code used for differential H3K4me3 breadth calling (called in the above scripts)
	- breadth_ds_coverage_BED.sh: generate original and downsampled coverage histograms (Input files are mapped reads, to big to share on GitHub)
	- compare_adjusted_peaks.pl: compare depth adjusted peaks
	- Compare_ECDFs_Cov.R: Functions to maximize Kolmogoroff-Smirnov p-value on ECDFs of coverage over meta peaks (selecting the 'least different' coverage)
	- Down_sampling_bed_file.pl: Down sampling of aligned reads  (Input files are mapped reads, to big to share on GitHub)
	- get_densities_per_DS_BED.pl: Get histograms for a downsampled file  (Input files are mapped reads, to big to share on GitHub)
	- Gradual_down_sampling_bed_file.pl: Get histograms for each gradually downsampled file  (Input files are mapped reads, to big to share on GitHub)
	- merge_homer_intersect_files.pl: Get merged file with annotations for R processing


* RNA_Mapping_counting_DE
	- map_with_STAR.sh: run STAR mapping on SCG cluster (example script for Heart samples)
	- count_reads.sh: Count reads within gene models (example script with heart samples); featureCounts is a function of the subreads package

	- analyze_all_samples_RNA_nested.R:RNAseq processing for tissue-independent aging processing
	- analyze_all_samples_RNA_forMDS_PCA.R:RNAseq processing for global MDS/PCA

	- process_RNAseq_datasets_v4_PCA.R: RNAseq DE processing
	- RNAseq_analysis_functions_v4_forPCA.R:functions called by above script

* SuperEnhancer_scripts: collection of code used for Super Enhancer calling (called in the above scripts)
	- merge_coverage_files.pl: merging coverage fles of each sample into a table for R processing
	- SuperEnhancerScripts.R: functions modified from ROSE to call inflection of coverage for super enhancers[1,2]
	
[1] Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes 
Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young
Cell 153, 307-319, April 11, 2013

[2] Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers 
Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young
Cell 153, 320-334, April 11, 2013

* Correlation_computation: scripts to compute sample-to-sample correlations for ChIP-seq and RNA-seq