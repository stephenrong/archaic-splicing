#!/bin/sh
tabix -h ../../../final-splicing-archaic-variants/data/annotate_popgen_gnomAD/gnomAD_AF/gnomad.genomes.r2.1.1.sites.passFILTER.dropINFO.updated.vcf.gz -R ../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_archaic_hub.bed | bgzip -c >| ../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_archaic_gnomAD_file.vcf.gz
