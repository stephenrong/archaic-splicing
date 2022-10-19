#!/bin/sh

RANDOM=111111

INPUT="../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_adaptive_hub.txt.gz"

# standard EUR_AC_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_moderate.txt.gz" -d $i"_aQTL" -o "../../results/enrichment_GTEx_QTLs_v2/adaptive_spliceai_moderate/adaptive_spliceai_moderate_matched_EUR_AC_bin_control_standard_SpliceAI" -n 1000 -s $RANDOM
done

# standard EUR_AC_binxEUR_hapR2tag_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard_hapR2.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_moderate.txt.gz" -d $i"_aQTL" -o "../../results/enrichment_GTEx_QTLs_v2/adaptive_spliceai_moderate/adaptive_spliceai_moderate_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_SpliceAI" -n 1000 -s $RANDOM
done
