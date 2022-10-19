#!/bin/sh
tabix -h ../../data/annotate_spliceai/spliceai_scores.raw.snv.hg19.vcf.gz -R ../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.9_hub.bed | bgzip >| ../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_greater0.9_spliceai_gs_raw.vcf.gz
tabix -h ../../data/annotate_spliceai/spliceai_scores.raw.snv.hg19.vcf.gz -R ../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hub.bed | bgzip >| ../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_spliceai_gs_raw.vcf.gz
