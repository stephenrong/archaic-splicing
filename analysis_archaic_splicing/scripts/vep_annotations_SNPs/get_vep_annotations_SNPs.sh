#!/bin/sh

# Get VEP annotations from vcf files

# VCF to VEP consequences
vep -i ../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.vcf.gz --cache --vcf --force_overwrite --pick --assembly GRCh37 --compress_output bgzip -o ../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.vcf.gz
vep -i ../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hub.vcf.gz --cache --vcf --force_overwrite --pick --assembly GRCh37 --compress_output bgzip -o ../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hub_vep.vcf.gz

tabix -p vcf ../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.vcf.gz
tabix -p vcf ../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hub_vep.vcf.gz
