#!/bin/sh

RANDOM=111111

INPUT="../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_adaptive_hub.txt.gz"

# standard EUR_AC_bin
for i in "missense_variant" "synonymous_variant" "5_prime_UTR_variant" "3_prime_UTR_variant" "intron_variant" "upstream_gene_variant" "downstream_gene_variant" "intergenic_variant"
do
	Rscript enrichment_GTEx_QTLs_standard.R -i $INPUT -a "../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.txt.gz" -d $i -o "../../results/enrichment_GTEx_QTLs_v2/adaptively_introgressed/adaptively_introgressed_matched_EUR_AC_bin_control_standard_VEP_consequence_"$i -n 100 -s $RANDOM
done

# standard EUR_AC_binxEUR_LDtagN_bin
for i in "missense_variant" "synonymous_variant" "5_prime_UTR_variant" "3_prime_UTR_variant" "intron_variant" "upstream_gene_variant" "downstream_gene_variant" "intergenic_variant"
do
	Rscript enrichment_GTEx_QTLs_standard_LDtagN.R -i $INPUT -a "../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.txt.gz" -d $i -o "../../results/enrichment_GTEx_QTLs_v2/adaptively_introgressed/adaptively_introgressed_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_VEP_consequence_"$i -n 100 -s $RANDOM
done
