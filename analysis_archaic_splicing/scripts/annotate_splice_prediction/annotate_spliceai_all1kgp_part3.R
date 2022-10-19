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
	vcf_spliceai <- vcf_spliceai %>% 
		dplyr::select(-c(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO))
	for (col in c("DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")) {
		vcf_spliceai[[col]] <- as.numeric(vcf_spliceai[[col]])
	}
	# get spliceai max value 
	vcf_spliceai_temp <- vcf_spliceai %>% dplyr::select("DS_AG", "DS_AL", "DS_DG", "DS_DL")
	vcf_spliceai_temp[is.na(vcf_spliceai_temp)] <- 0
	vcf_spliceai$spliceai_max <- do.call(`pmax`, vcf_spliceai_temp)
	return(vcf_spliceai)
}

ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	filter(((EUR_AF>0.01)&(EUR_AF<0.99)) | ((EAS_AF>0.01)&(EAS_AF<0.99)) | ((SAS_AF>0.01)&(SAS_AF<0.99)) | ((AMR_AF>0.01)&(AMR_AF<0.99)) | ((AFR_AF>0.01)&(AFR_AF<0.99)))
write_tsv(ALL_1KGP_phase3_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_hub.txt.gz"))
ALL_1KGP_phase3_hub <- read_tsv("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_hub.txt.gz")

# 	raw version
ALL_1KGP_phase3_spliceai_gs_raw <- as_tibble(read.vcfR("../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw.vcf.gz")@fix) 
ALL_1KGP_phase3_spliceai_gs_raw <- ALL_1KGP_phase3_spliceai_gs_raw %>% filter(POS %in% ALL_1KGP_phase3_hub$hub_variant_POS)
write_tsv(ALL_1KGP_phase3_spliceai_gs_raw, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_spliceai_gs_raw.txt.gz"))
ALL_1KGP_phase3_spliceai_gs_raw <- read_tsv("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_spliceai_gs_raw.txt.gz")
ALL_1KGP_phase3_spliceai_gs_raw <- ALL_1KGP_phase3_spliceai_gs_raw %>% process_spliceai()
write_tsv(ALL_1KGP_phase3_spliceai_gs_raw, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_process_spliceai_gs_raw.txt.gz"))
ALL_1KGP_phase3_spliceai_gs_raw <- read_tsv("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_process_spliceai_gs_raw.txt.gz")

# 	join tables
ALL_1KGP_phase3_hub_spliceai_gs_raw <- left_join(ALL_1KGP_phase3_hub, ALL_1KGP_phase3_spliceai_gs_raw)
write_tsv(ALL_1KGP_phase3_hub_spliceai_gs_raw, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_hub_spliceai_gs_raw.txt.gz"))
ALL_1KGP_phase3_hub_spliceai_gs_raw <- read_tsv("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_hub_spliceai_gs_raw.txt.gz")

# 	deduplicate
ALL_1KGP_phase3_hub_spliceai_gs_raw_dedup <- ALL_1KGP_phase3_hub_spliceai_gs_raw %>% 
	group_by(hub_variant_ID) %>% 
	filter((is.na(spliceai_max))|(spliceai_max == max(spliceai_max, na.rm=T))) %>% 
	ungroup()
write_tsv(ALL_1KGP_phase3_hub_spliceai_gs_raw_dedup, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_hub_spliceai_gs_raw_dedup.txt.gz"))
