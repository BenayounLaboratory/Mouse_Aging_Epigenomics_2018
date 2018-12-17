#!/bin/bash

if [[ "$#" -lt 1 ]]
then
    echo "$(basename $0) [BamDir] [oDir] "  1>&2
    exit 1
fi

BamDir=$(echo $1 | sed 's:/$::g')
oDir=$(echo $2 | sed 's:/$::g')

# make output directory if it doesnt exist
[[ ! -d "${oDir}" ]] && mkdir "${oDir}"

# run bowtie to map reads into SAM format
filePath="${BamDir}"
for f in $(find "$filePath" -name '*.bam')
do
	fileName=$(basename "${f}" | sed 's/\.bam/_counts_genes_mm9\.txt/g');
 	filePath="${oDir}"
    oFname="${filePath}/${fileName}"
    
    featureCounts -t exon -D 1500 -p --primary -T 3 -s 1 \
     -a /Volumes/MyBook_3/BD_aging_project/RNAseq/Cereb/1st_run/mm9_with_ERCC.gff \
     -o $oFname $f

done
echo "Finished read count\n"
