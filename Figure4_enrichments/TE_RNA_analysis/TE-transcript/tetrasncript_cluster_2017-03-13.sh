#!/bin/bash

# set the name of the job
#$ -N TEtranscript
#
# set the maximum memory usage
#$ -l h_vmem=8G
#
#$ -q extended
#
#$ -pe shm 6
#
# set the maximum run time
#$ -l h_rt=150:00:00
#
#$ -l h_stack=15M
#
# send mail when job ends or aborts
#$ -m ea
#
# specify an email address
#$ -M benayoun@stanford.edu
#
# check for errors in the job submission options
#$ -w e
#
#$ -R y

source ~/.bash_profile;

cd /srv/gsfs0/projects/brunet/Berenice/BB/TE_mapping;

TETEXEC=/srv/gsfs0/projects/brunet/Berenice/BB/TE_mapping/TEToolkit-1.5.1/bin/
export TETEXEC

$TETEXEC/TEtranscripts --project Liver_12vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/12m1_Liver_TGACCA_RAligned.out.bam ./NEW_STAR_map/12m2_Liver_ACAGTG_RAligned.out.bam ./NEW_STAR_map/12m3_Liver_GCCAAT_RAligned.out.bam \
   -c ./NEW_STAR_map/3m1_Liver_ATCACG_RAligned.out.bam ./NEW_STAR_map/3m2_Liver_CGATGT_RAligned.out.bam ./NEW_STAR_map/3m3_Liver_TTAGGC_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi
$TETEXEC/TEtranscripts --project Liver_29vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/29m1_Liver_CAGATC_RAligned.out.bam ./NEW_STAR_map/29m2_Liver_ACTTGA_RAligned.out.bam ./NEW_STAR_map/29m3_Liver_GATCAG_RAligned.out.bam \
   -c ./NEW_STAR_map/3m1_Liver_ATCACG_RAligned.out.bam ./NEW_STAR_map/3m2_Liver_CGATGT_RAligned.out.bam ./NEW_STAR_map/3m3_Liver_TTAGGC_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi


$TETEXEC/TEtranscripts --project Heart_12vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/12m4_heart_TGACCA_RAligned.out.bam ./NEW_STAR_map/12m5_heart_ACAGTG_RAligned.out.bam ./NEW_STAR_map/12m6_heart_GCCAAT_RAligned.out.bam \
   -c ./NEW_STAR_map/3m4_heart_ATCACG_RAligned.out.bam ./NEW_STAR_map/3m5_heart_CGATGT_RAligned.out.bam ./NEW_STAR_map/3m6_heart_TTAGGC_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi
$TETEXEC/TEtranscripts --project Heart_29vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/29m4_heart_CAGATC_RAligned.out.bam ./NEW_STAR_map/29m5_heart_ACTTGA_RAligned.out.bam ./NEW_STAR_map/29m6_heart_GATCAG_RAligned.out.bam \
   -c ./NEW_STAR_map/3m4_heart_ATCACG_RAligned.out.bam ./NEW_STAR_map/3m5_heart_CGATGT_RAligned.out.bam ./NEW_STAR_map/3m6_heart_TTAGGC_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi


$TETEXEC/TEtranscripts --project OB_12vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/OB_12m1_TGACCA_RAligned.out.bam ./NEW_STAR_map/OB_12m3_GCCAAT_RAligned.out.bam \
   -c ./NEW_STAR_map/OB_3m1_ATCACG_RAligned.out.bam ./NEW_STAR_map/OB_3m2_CGATGT_RAligned.out.bam ./NEW_STAR_map/OB_3m3_TTAGGC_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi
$TETEXEC/TEtranscripts --project OB_29vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/OB_29m2_CAGATC_RAligned.out.bam ./NEW_STAR_map/OB_29m3_ACTTGA_RAligned.out.bam ./NEW_STAR_map/OB_29m4_GATCAG_RAligned.out.bam \
   -c ./NEW_STAR_map/OB_3m1_ATCACG_RAligned.out.bam ./NEW_STAR_map/OB_3m2_CGATGT_RAligned.out.bam ./NEW_STAR_map/OB_3m3_TTAGGC_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi


$TETEXEC/TEtranscripts --project NPCs_12vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/NPC_12m5-ACAGTG_RAligned.out.bam ./NEW_STAR_map/NPC_12m6-GCCAAT_RAligned.out.bam \
   -c ./NEW_STAR_map/NPC_3m5-CGATGT_RAligned.out.bam ./NEW_STAR_map/NPC_3m6-TGACCA_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi
$TETEXEC/TEtranscripts --project NPCs_29vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/NPC_29m5-CAGATC_RAligned.out.bam ./NEW_STAR_map/NPC_29m6-CTTGTA_RAligned.out.bam \
   -c ./NEW_STAR_map/NPC_3m5-CGATGT_RAligned.out.bam ./NEW_STAR_map/NPC_3m6-TGACCA_RAligned.out.bam \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi


$TETEXEC/TEtranscripts --project Cerebellum_12vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/Cereb_12m1_1st_run_TGACCA_RAligned.out.bam ./NEW_STAR_map/Cereb_12m2_1st_run_ACAGTG_RAligned.out.bam ./NEW_STAR_map/Cereb_12m3_1st_run_GCCAAT_RAligned.out.bam ./NEW_STAR_map/Cereb_12m1_2nd_run_TGACCA_RAligned.out.bam ./NEW_STAR_map/Cereb_12m2_2nd_run_ACAGTG_RAligned.out.bam ./NEW_STAR_map/Cereb_12m3_2nd_run_GCCAAT_RAligned.out.bam  \
   -c ./NEW_STAR_map/Cereb_3m1_1st_run_ATCACG_RAligned.out.bam ./NEW_STAR_map/Cereb_3m2_1st_run_CGATGT_RAligned.out.bam ./NEW_STAR_map/Cereb_3m3_1st_run_TTAGGC_RAligned.out.bam ./NEW_STAR_map/Cereb_3m1_2nd_run_ATCACG_RAligned.out.bam ./NEW_STAR_map/Cereb_3m2_2nd_run_CGATGT_RAligned.out.bam ./NEW_STAR_map/Cereb_3m3_2nd_run_TTAGGC_RAligned.out.bam  \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi
$TETEXEC/TEtranscripts --project Cerebellum_29vs3m_TEtranscripts --format BAM \
   -t ./NEW_STAR_map/Cereb_29m2_1st_run_CAGATC_RAligned.out.bam ./NEW_STAR_map/Cereb_29m3_1st_run_ACTTGA_RAligned.out.bam ./NEW_STAR_map/Cereb_29m4_1st_run_GATCAG_RAligned.out.bam ./NEW_STAR_map/Cereb_29m2_2nd_run_CAGATC_RAligned.out.bam ./NEW_STAR_map/Cereb_29m3_2nd_run_ACTTGA_RAligned.out.bam ./NEW_STAR_map/Cereb_29m4_2nd_run_GATCAG_RAligned.out.bam  \
   -c ./NEW_STAR_map/Cereb_3m1_1st_run_ATCACG_RAligned.out.bam ./NEW_STAR_map/Cereb_3m2_1st_run_CGATGT_RAligned.out.bam ./NEW_STAR_map/Cereb_3m3_1st_run_TTAGGC_RAligned.out.bam ./NEW_STAR_map/Cereb_3m1_2nd_run_ATCACG_RAligned.out.bam ./NEW_STAR_map/Cereb_3m2_2nd_run_CGATGT_RAligned.out.bam ./NEW_STAR_map/Cereb_3m3_2nd_run_TTAGGC_RAligned.out.bam  \
   --GTF mm9_with_ERCC.gtf --TE mm9_rmsk_TE.gtf --mode multi
