#!/bin/sh
tabix -h ../../data/annotate_spliceai/spliceai_scores.raw.snv.hg19.vcf.gz -R ../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.bed | bgzip >| ../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_spliceai_gs_raw.vcf.gz
