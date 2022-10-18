#!/bin/sh

# STAR alignment for Neanderthal experiments.

Index=$1
InDirectory=$2
OutDirectory=$3
Prefix=$4
Postfix=$5

echo $Index
echo $InDirectory
echo $OutDirectory
echo $Prefix
echo $Postfix

# star_align
STAR \
	--runThreadN 8 \
	--genomeDir $Index \
	--readFilesIn $InDirectory/$Prefix\_R1_001.fastq.gz $InDirectory/$Prefix\_R2_001.fastq.gz \
	--readFilesCommand zcat -c \
	--outSAMtype BAM Unsorted \
	--outFileNamePrefix $OutDirectory/star_align/$Postfix\_001. \
	--alignIntronMax 1 \
	--outFilterMismatchNmax 5 \
	--alignEndsType EndToEnd \
	--outFilterMultimapNmax 1
