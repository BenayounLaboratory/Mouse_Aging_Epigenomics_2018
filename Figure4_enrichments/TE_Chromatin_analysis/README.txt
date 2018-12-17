#README

Analysis of chromatin changes at repeats: Similar approach to Liu N et al, Nature, 2018.

1. Use the mm9 RepeatMasker track from UCSC genome browser to get genomic position of repeats (Repeat_Masker_Reference/2018-09-12_Repeatmasker_mm9_ucsc.bed)

2. Intersect metapeaks with repeat library to reduce the number of regions in the genome that are under study (see Liu et al, 2018)
  Because we are using single end reads, mapping to repeats may be more delicate (since we have selected unique mappers)
  Thus, all of the elements that intersect metapeaks are extended 25bp upstream and downstream to improve capture of reads.
			Repeat_Masker_Reference/get_RM_intersections.sh

3. Generate a clean library of repeats to annotate diffbind with :
		* get one file of repeats with potential euchromatin:
			cat Repeat_Masker_Reference/Peaks_intersecting_repeats/*extend.bed | sort -u | sortBed -i - > RepeatMasker_intersecting_active_chromatin_elements.bed
		* Give unique names to each instance for downstream annotation
			cat RepeatMasker_intersecting_active_chromatin_elements.bed | perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[3]___$F[0]-$F[1]-$F[2]"' > RepeatMasker_intersecting_active_chromatin_elements.clean.bed

4. New diffbind (2.8.0) does not accept BED input anymore, and old diffbind could not be installed on new system.
  FIXseq BED reads are converted to bam using bedToBam (bedtools).

5. Use diffbind to generate normalized read counts at repeats:
		* Diffbind/Process_diffbind_for_data.R: code calling the processing functions
		* Diffbind/diffbind_process_functions.R: processing functions to obtain non NA peaks and normalized count matrix
		* Diffbind/Input and Diffbind/Output contain the sample sheet and output matrices and bed files.
		* Diffbind/Diffbind_annotate/annotate_diffbind_peaks.sh shows the code to extract annotations for DEseq2 (using the output of step 3)

6. Use DEseq2 to call differential histone mark intensity at repetitive elements
		* DEseq2/H3K4me3_TEs/Process_height_H3K4me3_linear_model_Repeats.R: code calling the processing functions
		* DEseq2/H3K4me3_TEs/height_analysis_functions_linear_vREPEATS.R: code to perform H3K4me3 DE at merged TE families
		
		* DEseq2/H3K27ac_TEs/Process_height_H3K27ac_linear_model_Repeats.R: code calling the processing functions
		* DEseq2/H3K27ac_TEs/height_analysis_functions_linear_K27ac_vREPEATS.R: code to perform H3K27ac DE at merged TE families
