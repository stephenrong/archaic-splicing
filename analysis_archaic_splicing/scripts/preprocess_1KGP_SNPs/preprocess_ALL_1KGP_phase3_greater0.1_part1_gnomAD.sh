#!/bin/sh
tabix -h ../../../final-splicing-archaic-variants/data/annotate_popgen_gnomAD/gnomAD_AF/gnomad.genomes.r2.1.1.sites.passFILTER.dropINFO.updated.vcf.gz -R ../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.1_hub.bed | bgzip -c >| ../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.1_gnomAD_file.vcf.gz
