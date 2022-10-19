#!/bin/R

# Take collated QTL files, add 1KGP SNP info

library(tidyverse)
library(data.table)
library(vcfR)

# load 1KGP variants
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub.txt.gz"))

# merge SpliceAI scores
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

# # 	raw version
ALL_1KGP_phase3_MAF0.01_spliceai_gs_raw <- as_tibble(read.vcfR("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_spliceai_gs_raw.vcf.gz")@fix) %>% process_spliceai()
ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw <- left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub, ALL_1KGP_phase3_MAF0.01_spliceai_gs_raw)
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw, gzfile("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw.gz"))

# # 	deduplicated version
ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw %>% 
	mutate(index = row_number())
ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_temp <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw %>% 
	dplyr::select(index, hub_variant_ID, spliceai_max)
ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_temp <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_temp %>% 
	arrange(desc(spliceai_max)) %>% mutate(unique = !duplicated(hub_variant_ID)) %>% filter(unique) %>% 
	dplyr::select(-unique)
ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw %>% 
	filter(index %in% ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_temp$index)
ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw %>% 
	dplyr::select(-index)
ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup %>% 
	dplyr::select(-index)

# # # 	deduplicated version
# ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw %>% 
# 	arrange(desc(spliceai_max)) %>% group_by(hub_variant_ID) %>% 
# 	slice(1) %>% ungroup()

write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup, gzfile("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup.gz"))

# split SpliceAI scores
collate_SpliceAI_weak <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup %>% 
	filter((spliceai_max >= 0.01)&(spliceai_max < 0.2))
collate_SpliceAI_moderate <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup %>% 
	filter((spliceai_max >= 0.2)&(spliceai_max < 0.5))
collate_SpliceAI_strong <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup %>% 
	filter((spliceai_max >= 0.5))
collate_SpliceAI_sig <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup %>% 
	filter((spliceai_max >= 0.2))

# save SpliceAI scores
write_tsv(collate_SpliceAI_weak, gzfile("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_weak.txt.gz"))
write_tsv(collate_SpliceAI_moderate, gzfile("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_moderate.txt.gz"))
write_tsv(collate_SpliceAI_strong, gzfile("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_strong.txt.gz"))
write_tsv(collate_SpliceAI_sig, gzfile("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_sig.txt.gz"))

# load SpliceAI files
collate_SpliceAI_weak <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_weak.txt.gz")) %>% 
	mutate(hub_variant_CHROM = as.character(hub_variant_CHROM), CHROM = as.character(CHROM))
collate_SpliceAI_moderate <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_moderate.txt.gz")) %>% 
	mutate(hub_variant_CHROM = as.character(hub_variant_CHROM), CHROM = as.character(CHROM))
collate_SpliceAI_strong <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_strong.txt.gz")) %>% 
	mutate(hub_variant_CHROM = as.character(hub_variant_CHROM), CHROM = as.character(CHROM))
collate_SpliceAI_sig <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/collate_SpliceAI_sig.txt.gz")) %>% 
	mutate(hub_variant_CHROM = as.character(hub_variant_CHROM), CHROM = as.character(CHROM))

# join 1KGP variants
join_in1KGP_GTEx_weak <- collate_SpliceAI_weak %>% left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub)
join_in1KGP_GTEx_moderate <- collate_SpliceAI_moderate %>% left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub)
join_in1KGP_GTEx_strong <- collate_SpliceAI_strong %>% left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub)
join_in1KGP_GTEx_sig <- collate_SpliceAI_sig %>% left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub)

# filter 1KGP variants
filter_in1KGP_SpliceAI_bins_weak <- join_in1KGP_GTEx_weak %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only
filter_in1KGP_SpliceAI_bins_moderate <- join_in1KGP_GTEx_moderate %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only
filter_in1KGP_SpliceAI_bins_strong <- join_in1KGP_GTEx_strong %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only
filter_in1KGP_SpliceAI_bins_sig <- join_in1KGP_GTEx_sig %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only

# save filter 1KGP variants
write_tsv(filter_in1KGP_SpliceAI_bins_weak, "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_weak.txt.gz")
write_tsv(filter_in1KGP_SpliceAI_bins_moderate, "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_moderate.txt.gz")
write_tsv(filter_in1KGP_SpliceAI_bins_strong, "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_strong.txt.gz")
write_tsv(filter_in1KGP_SpliceAI_bins_sig, "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_sig.txt.gz")
