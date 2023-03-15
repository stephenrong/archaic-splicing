#!/bin/sh

RANDOM=111111

INPUT="../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_introgressed_hub.txt.gz"

# standard EUR_AC_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_strong.txt.gz" -d $i -o "../../results/enrichment_GTEx_QTLs_v2/introgressed_spliceai_strong/introgressed_spliceai_strong_matched_EUR_AC_bin_control_standard_SpliceAI" -n 1000 -s $RANDOM
done

# standard EUR_AC_binxEUR_LDtagN_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard_LDtagN.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_strong.txt.gz" -d $i -o "../../results/enrichment_GTEx_QTLs_v2/introgressed_spliceai_strong/introgressed_spliceai_strong_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_SpliceAI" -n 1000 -s $RANDOM
done
