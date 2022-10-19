#!/bin/R

library(tidyverse)
library(data.table)
library(vcfR)

# load tables
process_spliceai <- function(vcf_spliceai) {
	vcf_spliceai <- vcf_spliceai %>% 
		mutate(hub_variant_ID = paste(as.character(CHROM), as.character(POS), paste(REF, ALT, sep="/"), sep="_"))
	vcf_spliceai <- vcf_spliceai %>% 
		separate(col=INFO, into=c("ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"), remove=F, sep="\\|")
	for (col in c("DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")) {
		vcf_spliceai[[col]] <- as.numeric(vcf_spliceai[[col]])
	}
	vcf_spliceai <- vcf_spliceai %>% 
		dplyr::select(-c(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO))
	# get spliceai max value 
	vcf_spliceai_temp <- vcf_spliceai %>% dplyr::select("DS_AG", "DS_AL", "DS_DG", "DS_DL")
	vcf_spliceai_temp[is.na(vcf_spliceai_temp)] <- 0
	vcf_spliceai$spliceai_max <- do.call(`pmax`, vcf_spliceai_temp)
	return(vcf_spliceai)
}

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.txt.gz"))

# 	raw version
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_spliceai_gs_raw <- as_tibble(read.vcfR("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_spliceai_gs_raw.vcf.gz")@fix) %>% process_spliceai()
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw <- left_join(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_spliceai_gs_raw)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw.txt.gz"))

# 	deduplicate
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw %>% 
	arrange(desc(spliceai_max)) %>% group_by(hub_variant_ID) %>% 
	slice(1) %>% ungroup()
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))

# # 	masked version
# final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_spliceai_gs_masked <- as_tibble(read.vcfR("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_spliceai_gs_masked.vcf.gz")@fix) %>% process_spliceai()
# final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_masked <- left_join(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_spliceai_gs_masked)
# write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_masked, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_masked.txt.gz"))

# # 	deduplicate
# final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_masked_dedup <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_masked %>% 
# 	arrange(desc(spliceai_max)) %>% group_by(hub_variant_ID) %>% 
# 	slice(1) %>% ungroup()
# write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_masked_dedup, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_masked_dedup.txt.gz"))
