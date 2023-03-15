#!/bin/sh

RANDOM=111111

INPUT="../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_adaptive_hub.txt.gz"

# standard EUR_AC_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_weak.txt.gz" -d $i -o "../../results/enrichment_GTEx_QTLs_v2/adaptive_spliceai_weak/adaptive_spliceai_weak_matched_EUR_AC_bin_control_standard_SpliceAI" -n 1000 -s $RANDOM
done

# standard EUR_AC_binxEUR_LDtagN_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard_LDtagN.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_weak.txt.gz" -d $i -o "../../results/enrichment_GTEx_QTLs_v2/adaptive_spliceai_weak/adaptive_spliceai_weak_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_SpliceAI" -n 1000 -s $RANDOM
done
