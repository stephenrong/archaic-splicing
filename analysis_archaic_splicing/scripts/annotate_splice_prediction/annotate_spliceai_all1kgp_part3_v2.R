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
		mutate(ALLELE = gsub("SpliceAI=", "", ALLELE))
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

ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))

# 	raw version
ALL_1KGP_phase3_spliceai_gs_raw <- as_tibble(read.vcfR("../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw.vcf.gz")@fix) %>% process_spliceai()
ALL_1KGP_phase3_hub_spliceai_gs_raw <- left_join(ALL_1KGP_phase3_hub, ALL_1KGP_phase3_spliceai_gs_raw)
write_tsv(ALL_1KGP_phase3_hub_spliceai_gs_raw, gzfile("../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub_spliceai_gs_raw.txt.gz"))

# 	deduplicate
ALL_1KGP_phase3_hub_spliceai_gs_raw_dedup <- ALL_1KGP_phase3_hub_spliceai_gs_raw %>% 
	arrange(desc(spliceai_max)) %>% group_by(hub_variant_ID) %>% 
	slice(1) %>% ungroup()
write_tsv(ALL_1KGP_phase3_hub_spliceai_gs_raw_dedup, gzfile("../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub_spliceai_gs_raw_dedup.txt.gz"))

# 	reload
ALL_1KGP_phase3_hub_spliceai_gs_raw <- as_tibble(fread("../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub_spliceai_gs_raw.txt.gz"))
ALL_1KGP_phase3_hub_spliceai_gs_raw_dedup <- as_tibble(fread("../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub_spliceai_gs_raw_dedup.txt.gz"))

# common
ALL_1KGP_phase3_common_hub_spliceai_gs_raw <- ALL_1KGP_phase3_hub_spliceai_gs_raw %>% 
	filter(((EUR_AF>0.01)&(EUR_AF<0.99)) | ((EAS_AF>0.01)&(EAS_AF<0.99)) | ((SAS_AF>0.01)&(SAS_AF<0.99)) | ((AMR_AF>0.01)&(AMR_AF<0.99)) | ((AFR_AF>0.01)&(AFR_AF<0.99)))
ALL_1KGP_phase3_common_hub_spliceai_gs_raw_dedup <- ALL_1KGP_phase3_common_hub_spliceai_gs_raw %>% 
	arrange(desc(spliceai_max)) %>% group_by(hub_variant_ID) %>% 
	slice(1) %>% ungroup()
write_tsv(ALL_1KGP_phase3_common_hub_spliceai_gs_raw, gzfile("../../results/annotate_splice_prediction/ALL_1KGP_phase3_common_hub_spliceai_gs_raw.txt.gz"))
write_tsv(ALL_1KGP_phase3_common_hub_spliceai_gs_raw_dedup, gzfile("../../results/annotate_splice_prediction/ALL_1KGP_phase3_common_hub_spliceai_gs_raw_dedup.txt.gz"))

# MAF0.01
ALL_1KGP_phase3_MAF0.01_hub_spliceai_gs_raw <- ALL_1KGP_phase3_hub_spliceai_gs_raw %>% 
	filter(EUR_AF >= 0.01, EUR_AF <= 0.99)
ALL_1KGP_phase3_MAF0.01_hub_spliceai_gs_raw_dedup <- ALL_1KGP_phase3_MAF0.01_hub_spliceai_gs_raw %>% 
	arrange(desc(spliceai_max)) %>% group_by(hub_variant_ID) %>% 
	slice(1) %>% ungroup()
write_tsv(ALL_1KGP_phase3_MAF0.01_hub_spliceai_gs_raw, gzfile("../../results/annotate_splice_prediction/ALL_1KGP_phase3_MAF0.01_hub_spliceai_gs_raw.txt.gz"))
write_tsv(ALL_1KGP_phase3_MAF0.01_hub_spliceai_gs_raw_dedup, gzfile("../../results/annotate_splice_prediction/ALL_1KGP_phase3_MAF0.01_hub_spliceai_gs_raw_dedup.txt.gz"))
