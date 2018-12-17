#!/bin/bash

# set the name of the job
#$ -N DF_star
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


cd /srv/gsfs0/projects/brunet/Berenice/BB/Nucleosomes/;

export DINUPEXEC=/srv/gsfs0/projects/brunet/Berenice/BB/Nucleosomes/dinup_1.3/bin

$DINUPEXEC/dinup --fdr=0.01 --treatment=29m_Cerebellum_H3_MERGED.FIXSEQ_CLEANED_READS.bed --control=3m_Cerebellum_H3_MERGED.FIXSEQ_CLEANED_READS.bed --name=3v29m_Cerebellum --feature --format='BED' --wig --region=100
$DINUPEXEC/dinup --fdr=0.01 --treatment=29m_Heart_H3_MERGED.FIXSEQ_CLEANED_READS.bed --control=3m_Heart_H3_MERGED.FIXSEQ_CLEANED_READS.bed --name=3v29m_Heart --feature --format='BED' --wig --region=100
$DINUPEXEC/dinup --fdr=0.01 --treatment=29m_Liver_H3_MERGED.FIXSEQ_CLEANED_READS.bed --control=3m_Liver_H3_MERGED.FIXSEQ_CLEANED_READS.bed --name=3v29m_Liver --feature --format='BED' --wig --region=100
$DINUPEXEC/dinup --fdr=0.01 --treatment=29m_NPCs_H3_MERGED.FIXSEQ_CLEANED_READS.bed --control=3m_NPCs_H3_MERGED.FIXSEQ_CLEANED_READS.bed --name=3v29m_NPCs --feature   --format='BED' --wig --region=100
$DINUPEXEC/dinup --fdr=0.01 --treatment=29m_OB_H3_MERGED.FIXSEQ_CLEANED_READS.bed --control=3m_OB_H3_MERGED.FIXSEQ_CLEANED_READS.bed --name=3v29m_OB --feature --format='BED' --wig --region=100
