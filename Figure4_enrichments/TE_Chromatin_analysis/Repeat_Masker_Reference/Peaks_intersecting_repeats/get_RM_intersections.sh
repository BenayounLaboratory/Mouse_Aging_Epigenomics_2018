# 2018-09-12
# intersect meta peaks to repeat masker
# similar method than 
# Selective silencing of euchromatic L1s revealed by genome-wide screens for L1 regulators
# Nian Liu, Cameron H. Lee, Tomek Swigut, Edward Grow, Bo Gu, Michael C. Bassik & Joanna Wysocka
# Nature volume 553, pages 228–232 (11 January 2018)

# intersecting elements
# require that 50% of the elements overlaps with a peak
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Cerebellum/H3K4me3_breadth/ALL_AGES_MERGED_Cerebellum_H3K4me3_peaks.bed                        > Repeatmasker_mm9_Cerebellum_H3K4me3_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Heart/H3K4me3_breadth/ALL_AGES_MERGED_Heart_H3K4me3_peaks.bed                                  > Repeatmasker_mm9_Heart_H3K4me3_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Liver/H3K4me3_breadth/ALL_AGES_MERGED_Liver_H3K4me3_broad_peaks.bed                            > Repeatmasker_mm9_Liver_H3K4me3_broad_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Olfactory_Bulb/H3K4me3_Breadth/ALL_AGES_MERGED_OB_H3K4me3_broad_peaks.bed                      > Repeatmasker_mm9_OB_H3K4me3_broad_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/NPC_cultures/New_rep_5_6/Combined_breadth/ALL_AGES_MERGED_CombinedNPC_H3K4me3_broad_peaks.bed  > Repeatmasker_mm9_CombinedNPC_H3K4me3_broad_peaks_intersecting_Repeats.bed

intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Cerebellum/SuperEnhancer/ALL_AGES_MERGED_Cerebellum_H3K27ac_broad_peaks.bed                    > Repeatmasker_mm9_Cerebellum_H3K27ac_broad_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Heart/SuperEnhancer/ALL_AGES_MERGED_Heart_H3K27ac_broad_peaks.bed                              > Repeatmasker_mm9_Heart_H3K27ac_broad_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Liver/SuperEnhancer/ALL_AGES_MERGED_Liver_H3K27ac_broad_peaks.bed                              > Repeatmasker_mm9_Liver_H3K27ac_broad_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/Olfactory_Bulb/SuperEnhancer/ALL_AGES_MERGED_OB_H3K27ac_PooledLanes_broad_peaks.bed            > Repeatmasker_mm9_OB_H3K27ac_PooledLanes_broad_peaks_intersecting_Repeats.bed
intersectBed -u -f 0.5 -a ../2018-09-12_Repeatmasker_mm9_ucsc.bed -b /Volumes/BB_Backup_3/BD_aging_project/ChIP-seq/NPC_cultures/SuperEnhancer/ALL_AGES_MERGED_NPCs_H3K27ac_PooledLanes_broad_peaks.bed            > Repeatmasker_mm9_NPCs_H3K27ac_PooledLanes_broad_peak_intersecting_Repeats.bed


# extend to get surrounding to get more unique mappers (surrounding bases pairs, 25 each way)
slopBed -b 100 -i Repeatmasker_mm9_Cerebellum_H3K4me3_peaks_intersecting_Repeats.bed             -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > Cereb_H3K4me3_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_Heart_H3K4me3_peaks_intersecting_Repeats.bed                  -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > Heart_H3K4me3_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_Liver_H3K4me3_broad_peaks_intersecting_Repeats.bed            -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > Liver_H3K4me3_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_OB_H3K4me3_broad_peaks_intersecting_Repeats.bed               -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > OB_H3K4me3_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_CombinedNPC_H3K4me3_broad_peaks_intersecting_Repeats.bed      -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > NSPCs_H3K4me3_Repeatmasker_mm9_25bp_extend.bed

slopBed -b 100 -i Repeatmasker_mm9_Cerebellum_H3K27ac_broad_peaks_intersecting_Repeats.bed       -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > Cereb_H3K27ac_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_Heart_H3K27ac_broad_peaks_intersecting_Repeats.bed            -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > Heart_H3K27ac_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_Liver_H3K27ac_broad_peaks_intersecting_Repeats.bed            -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > Liver_H3K27ac_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_OB_H3K27ac_PooledLanes_broad_peaks_intersecting_Repeats.bed   -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > OB_H3K27ac_Repeatmasker_mm9_25bp_extend.bed
slopBed -b 100 -i Repeatmasker_mm9_NPCs_H3K27ac_PooledLanes_broad_peak_intersecting_Repeats.bed  -g /Users/benayoun/Softwares/bedtools2.26/genomes/mouse.mm9.genome | sortBed -i - > NSPCs_H3K27ac_Repeatmasker_mm9_25bp_extend.bed





