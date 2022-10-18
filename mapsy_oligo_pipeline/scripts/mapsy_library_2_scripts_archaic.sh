#!/bin/sh
# Script file for generating Neanderthal Updated Library oligos.
Rscript mapsy_create_v1.R \
	-p mapsy_library_2_params.R \
	-f ../data/fasta/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz \
	-g ../data/gff/gencode.v32lift37.basic.annotation.gff3.gz \
	-v ../data/merge_variants/merge_variants_library_hub.vcf.gz \
	-o ../results/mapsy_libraries/library_2/Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2 \
	-d Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2
