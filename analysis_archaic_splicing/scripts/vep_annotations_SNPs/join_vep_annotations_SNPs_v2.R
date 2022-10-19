#!/bin/R

# Join VEP annotations to hub files

library(tidyverse)
library(data.table)
library(vcfR)
# source("../get_helper.R")

process_vep <- function(vcf_vep) {
	vcf_vep <- vcf_vep %>% 
		mutate(VEP = gsub("\\|.*", "", gsub("^.*?\\|", "", gsub(".*CSQ=", "", INFO))))
	for (i in c("stop_gained", "stop_lost", "missense_variant", "synonymous_variant", 
			"splice_donor_variant", "splice_acceptor_variant", "splice_region_variant", 
			"5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant", 
			"upstream_gene_variant", "downstream_gene_variant", 
			"intergenic_variant")) {
		vcf_vep[[i]] <- 
			ifelse(grepl(i, vcf_vep$VEP), TRUE, FALSE)
	}
	# vcf_vep <- vcf_vep %>% 
	# 	mutate(hub_variant_ID = ID)
	# vcf_vep <- vcf_vep %>% dplyr::select(-VEP) %>% 
	# 	dplyr::select(-c(CHROM, POS, REF, ALT, QUAL, FILTER, INFO))
	vcf_vep <- vcf_vep %>% dplyr::select(-VEP) %>% 
		dplyr::select(-c(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO))
	return(vcf_vep)
}

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- as_tibble(read.vcfR("../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.vcf.gz")@fix) %>% process_vep()
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_final <- bind_cols(as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz")), final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_final, "../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.txt.gz")

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep <- as_tibble(read.vcfR("../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.vcf.gz")@fix) %>% process_vep()
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_final <- bind_cols(as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub.txt.gz")), final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_final, "../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep.txt.gz")

ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep <- as_tibble(read.vcfR("../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hub_vep.vcf.gz")@fix) %>% process_vep()
ALL_1KGP_phase3_MAF0.01_hapR2_hub_final <- bind_cols(as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub.txt.gz")), ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep)
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_final, "../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.txt.gz")
