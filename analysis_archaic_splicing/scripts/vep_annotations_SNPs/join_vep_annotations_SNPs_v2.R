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
	vcf_vep <- vcf_vep %>% dplyr::select(-VEP) %>% 
		mutate(hub_variant_CHROM = CHROM, hub_variant_POS = POS, hub_variant_REF = REF, hub_variant_ALT = ALT) %>% 
		dplyr::select(-c(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO))
	return(vcf_vep)
}

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- as_tibble(read.vcfR("../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.vcf.gz")@fix) %>% process_vep() %>% 
	mutate(across(everything(), as.character))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz")) %>% 
	mutate(across(everything(), as.character))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_final <- full_join(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_final, "../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.txt.gz")

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep <- as_tibble(read.vcfR("../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.vcf.gz")@fix) %>% process_vep() %>% 
	mutate(across(everything(), as.character))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub.txt.gz")) %>% 
	mutate(across(everything(), as.character))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_final <- full_join(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_final, "../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep.txt.gz")

ALL_1KGP_phase3_MAF0.01_hub_vep <- as_tibble(read.vcfR("../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hub_vep.vcf.gz")@fix) %>% process_vep() %>% 
	mutate(across(everything(), as.character))
ALL_1KGP_phase3_MAF0.01_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hub.txt.gz")) %>% 
	mutate(across(everything(), as.character))
ALL_1KGP_phase3_MAF0.01_hub_final <- full_join(ALL_1KGP_phase3_MAF0.01_hub, ALL_1KGP_phase3_MAF0.01_hub_vep)
write_tsv(ALL_1KGP_phase3_MAF0.01_hub_final, "../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hub_vep.txt.gz")

ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep <- as_tibble(read.vcfR("../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.vcf.gz")@fix) %>% process_vep() %>% 
	mutate(across(everything(), as.character))
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub.txt.gz")) %>% 
	mutate(across(everything(), as.character))
ALL_1KGP_phase3_MAF0.01_hapR2_hub_final <- full_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub, ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep)
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_final, "../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.txt.gz")
