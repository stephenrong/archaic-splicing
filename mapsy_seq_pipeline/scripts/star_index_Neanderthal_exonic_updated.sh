#!/bin/sh
# conda activate splicing-archaic-mapsy-v3
# conda env export > splicing-archaic-mapsy-v3.yaml

gunzip -fkv ../data/custom_reference/*.fasta.gz

# concat phix and neanderthal updated references
cat ../data/custom_reference/NC_001422.fna \
../data/custom_reference/neanderthal_updated_input.fasta >| \
../data/custom_reference/neanderthal_updated_phix_input.fasta
cat ../data/custom_reference/NC_001422.fna \
../data/custom_reference/neanderthal_updated_output.fasta >| \
../data/custom_reference/neanderthal_updated_phix_output.fasta

# create star index files updated
STAR \
	--runThreadN 32 \
	--runMode genomeGenerate \
	--genomeDir ../data/custom_star_index/neanderthal_updated_phix_input \
	--genomeFastaFiles ../data/custom_reference/neanderthal_updated_phix_input.fasta \
	--genomeSAindexNbases 8 \
	--genomeChrBinNbits 8

STAR \
	--runThreadN 32 \
	--runMode genomeGenerate \
	--genomeDir ../data/custom_star_index/neanderthal_updated_phix_output \
	--genomeFastaFiles ../data/custom_reference/neanderthal_updated_phix_output.fasta \
	--genomeSAindexNbases 8 \
	--genomeChrBinNbits 8

STAR \
	--runThreadN 32 \
	--runMode genomeGenerate \
	--genomeDir ../data/custom_star_index/neanderthal_updated_input \
	--genomeFastaFiles ../data/custom_reference/neanderthal_updated_input.fasta \
	--genomeSAindexNbases 8 \
	--genomeChrBinNbits 8

STAR \
	--runThreadN 32 \
	--runMode genomeGenerate \
	--genomeDir ../data/custom_star_index/neanderthal_updated_output \
	--genomeFastaFiles ../data/custom_reference/neanderthal_updated_output.fasta \
	--genomeSAindexNbases 8 \
	--genomeChrBinNbits 8
