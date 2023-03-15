#!/bin/sh

RANDOM=111111

INPUT="../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_introgressed_hub.txt.gz"

# standard EUR_AC_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_GTEx_sQTLs_hg19.txt.gz" -d $i"_sQTL" -o "../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_bin_control_standard_sQTLs_"$i -n 1000 -s $RANDOM
done

for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_GTEx_eQTLs_hg19.txt.gz" -d $i"_eQTL" -o "../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_bin_control_standard_eQTLs_"$i -n 1000 -s $RANDOM
done

# standard EUR_AC_binxEUR_LDtagN_bin
for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard_LDtagN.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_GTEx_sQTLs_hg19.txt.gz" -d $i"_sQTL" -o "../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_sQTLs_"$i -n 1000 -s $RANDOM
done

for i in "NULL"
do
	Rscript enrichment_GTEx_QTLs_standard_LDtagN.R -i $INPUT -a "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_GTEx_eQTLs_hg19.txt.gz" -d $i"_eQTL" -o "../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_eQTLs_"$i -n 1000 -s $RANDOM
done
