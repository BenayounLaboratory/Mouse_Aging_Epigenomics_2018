
# 2017-01-04

# README file became empty after computer crash

# from memory:

run TE_transcript on cluster

# 2017-01-04
# download count tables to process in DEseq2

# to separate transposons and non transposons
cat Count_tables/Cereb_29vs3m_TEtranscripts.cntTable | cut -f 1 > gene_TE_names.txt

+ manual selection of the TE names and save as 


# 2017-03-13
rerun mapping with laxer multimapping data (cf TE Transcript).
allow 150 multimapping
tetrasncript_cluster_2017-03-13.sh

had to use untrimmed reads (otherwise didin't run)

# 2017-03-15
get TEtranscript results


# 2017-09-12
# parse TE names for grouped analysis
perl parse_TE_families.pl TE_names_list.txt > 2017-09-12_Parsed_TE_names_list.txt