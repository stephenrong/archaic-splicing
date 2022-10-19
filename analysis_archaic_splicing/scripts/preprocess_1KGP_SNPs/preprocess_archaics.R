#!/bin/R

# Preprocess archaic data (get B, get mask, convert to hub format, join archaics together, add 1KGP SNP info, add gnomAD SNP info [convert to BED, tabix gnomAD, join gnomAD])

library(tidyverse)
library(data.table)
library(plyranges)
library(wrapr)
library(vcfR)
library(rtracklayer)
source("../get_helper.R")

# load
altai_denisovan <- read.vcfR("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/altai_denisovan/altai_denisovan_alt_bedfil_norm.vcf.gz")
altai_neanderthal <- read.vcfR("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/altai_neanderthal/altai_neanderthal_alt_bedfil_norm.vcf.gz")
chagyrskaya_neanderthal <- read.vcfR("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/chagyrskaya_neanderthal/chagyrskaya_neanderthal_alt_bedfil_norm.vcf.gz")
vindija_neanderthal <- read.vcfR("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/vindija_neanderthal/vindija_neanderthal_alt_bedfil_norm.vcf.gz")

# helper function
process_archaic_vcf <- function(archaic_vcf, archaic_vcf_output) {
	print(archaic_vcf_output)

	# clean up fix
	print("clean up fix")
	archaic_vcf_hub <- archaic_vcf@fix %>% as_tibble()
	archaic_vcf_hub <- archaic_vcf_hub %>% 
		separate(col=INFO, into=c("AC", "AN"), sep=";") %>% 
		mutate(AC = as.numeric(gsub(".*=", "", AC))) %>% 
		mutate(AN = as.numeric(gsub(".*=", "", AN)))

	# clean and bind gt
	print("clean and bind gt")
	archaic_vcf_gt <- archaic_vcf@gt %>% as_tibble()
	archaic_vcf_gt <- archaic_vcf_gt[,2]
	names(archaic_vcf_gt) <- c("GT")
	archaic_vcf_gt <- archaic_vcf_gt %>% 
		mutate(GT = gsub(":.*", "", GT))
	archaic_vcf_hub <- archaic_vcf_hub %>% 
		bind_cols(archaic_vcf_gt)

	# # filter to SNPs, unnecessary
	# print("filter to SNPs")
	# table(archaic_vcf_hub$ALT)
	# archaic_vcf_hub <- archaic_vcf_hub %>% 
	# 	filter(
	# 		REF %in% c("A", "C", "T", "G"),
	# 		ALT %in% c("A", "C", "T", "G")
	# 	)

	# convert to GRange
	print("convert to GRange")
	archaic_vcf_hub <- archaic_vcf_hub %>% 
		mutate(POS = as.numeric(POS))

	archaic_vcf_gr <- archaic_vcf_hub %>% 
		dplyr::rename(
			seqnames = CHROM
		) %>% 
		mutate(
			seqnames = paste("chr", seqnames, sep=""), 
			start = POS, 
			end = POS) %>% 
		GRanges()

	# get B statistics
	print("get B statistics")
	B_stats_hg19_gr <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/B_stats_hg19_gr.txt.gz")))
	# 	possible for liftOver to result in multiple overlapping regions
	# 	in this case, just pick the first such overlapping region
	archaic_overlap <- findOverlaps(archaic_vcf_gr, B_stats_hg19_gr) %>% 
		as_tibble() %>% filter(!duplicated(queryHits))
	archaic_overlap_index <- archaic_overlap$queryHits
	archaic_overlap_B_statistic <- B_stats_hg19_gr$B_statistic[archaic_overlap$subjectHits]
	archaic_overlap_B_statistic_bin <- B_stats_hg19_gr$B_statistic_bin[archaic_overlap$subjectHits]
	# 	now append B statistic information to table based on index
	archaic_vcf_hub$B_statistic <- NA
	archaic_vcf_hub$B_statistic_bin <- NA
	archaic_vcf_hub$B_statistic[archaic_overlap_index] <- archaic_overlap_B_statistic
	archaic_vcf_hub$B_statistic_bin[archaic_overlap_index] <- archaic_overlap_B_statistic_bin

	# get mask overlap
	print("get inter mask overlap")
	archaic_mask_inter <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_inter.bed.txt.gz")))
	seqlevels(archaic_mask_inter) <- paste("chr", seqlevels(archaic_mask_inter), sep="")
	mask_overlap <- findOverlaps(archaic_vcf_gr, archaic_mask_inter) %>% 
		as_tibble() %>% filter(!duplicated(queryHits))
	mask_index <- mask_overlap$queryHits
	archaic_vcf_hub$archaic_mask_inter <- FALSE 
	archaic_vcf_hub$archaic_mask_inter[mask_index] <- TRUE

	# get mask overlap
	print("get union mask overlap")
	archaic_mask_union <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_union.bed.txt.gz")))
	seqlevels(archaic_mask_union) <- paste("chr", seqlevels(archaic_mask_union), sep="")
	mask_overlap <- findOverlaps(archaic_vcf_gr, archaic_mask_union) %>% 
		as_tibble() %>% filter(!duplicated(queryHits))
	mask_index <- mask_overlap$queryHits
	archaic_vcf_hub$archaic_mask_union <- FALSE 
	archaic_vcf_hub$archaic_mask_union[mask_index] <- TRUE

	# convert to hub
	print("convert to hub")
	archaic_vcf_hub <- archaic_vcf_hub %>% 
		mutate(
			hub_reference_genome = "hg19", 
			hub_variant_CHROM = CHROM, 
			hub_variant_POS = POS, 
			hub_variant_ID = NA, 
			hub_variant_REF = REF, 
			hub_variant_ALT = ALT
		)

	# add ancestral
	print("add ancestral")
	archaic_vcf_hub <- archaic_vcf_hub %>% 
		mutate_hub_variant_ANC() %>% 
		mutate_hub_variant_DER() %>% 
		mutate_hub_variant_ID() %>% 
		split_ancestral() %>% 
		standard_sort()

	# save file
	print("save file")
	archaic_vcf_hub <- archaic_vcf_hub %>% 
		dplyr::select(starts_with("hub_"), !starts_with("hub_"))
	write_tsv(archaic_vcf_hub, archaic_vcf_output)
}

# iterate through archaic vcfs
process_archaic_vcf(altai_denisovan, gzfile("../../results/preprocess_1KGP_SNPs/altai_denisovan_hub.txt.gz"))
process_archaic_vcf(altai_neanderthal, gzfile("../../results/preprocess_1KGP_SNPs/altai_neanderthal_hub.txt.gz"))
process_archaic_vcf(chagyrskaya_neanderthal, gzfile("../../results/preprocess_1KGP_SNPs/chagyrskaya_neanderthal_hub.txt.gz"))
process_archaic_vcf(vindija_neanderthal, gzfile("../../results/preprocess_1KGP_SNPs/vindija_neanderthal_hub.txt.gz"))

# create filters
# strip_cols <- function(input) {
# 	input %>% dplyr::select(starts_with("hub_"))
# 	return(input)
# }

prefix_cols <- function(input, prefix, sep="_") {
	temp <- which(!(grepl("hub|B_statistic|archaic_mask", names(input))))
	names(input)[temp] <- paste(prefix, names(input)[temp], sep=sep)
	return(input)
}

# 	load archaics
altai_denisovan <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/altai_denisovan_hub.txt.gz")) %>% prefix_cols("altai_denisovan")
altai_neanderthal <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/altai_neanderthal_hub.txt.gz")) %>% prefix_cols("altai_neanderthal")
chagyrskaya_neanderthal <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/chagyrskaya_neanderthal_hub.txt.gz")) %>% prefix_cols("chagyrskaya_neanderthal")
vindija_neanderthal <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/vindija_neanderthal_hub.txt.gz")) %>% prefix_cols("vindija_neanderthal")

# 	join archaics
archaics_all_hub <- altai_denisovan %>% full_join(altai_neanderthal) %>% full_join(chagyrskaya_neanderthal) %>% full_join(vindija_neanderthal)
archaics_all_hub <- archaics_all_hub %>% rowwise() %>% 
	mutate(archaic_AC = sum(c(altai_denisovan_AC, altai_neanderthal_AC, chagyrskaya_neanderthal_AC, vindija_neanderthal_AC), na.rm=T)) %>% 
	mutate(archaic_AN = sum(c(altai_denisovan_AN, altai_neanderthal_AN, chagyrskaya_neanderthal_AN, vindija_neanderthal_AN), na.rm=T)) %>% 
	mutate(neanderthal_AC = sum(c(altai_neanderthal_AC, chagyrskaya_neanderthal_AC, vindija_neanderthal_AC), na.rm=T)) %>% 
	mutate(neanderthal_AN = sum(c(altai_neanderthal_AN, chagyrskaya_neanderthal_AN, vindija_neanderthal_AN), na.rm=T)) %>% 
	mutate(denisovan_AC = sum(c(altai_denisovan_AC), na.rm=T)) %>% 
	mutate(denisovan_AN = sum(c(altai_denisovan_AN), na.rm=T)) %>% 
	ungroup() %>% 
	dplyr::select(c(starts_with("hub_"), 
		altai_denisovan_AC, altai_denisovan_AN, 
		altai_neanderthal_AC, altai_neanderthal_AN, 
		chagyrskaya_neanderthal_AC, chagyrskaya_neanderthal_AN, 
		vindija_neanderthal_AC, vindija_neanderthal_AN, 
		archaic_AC, archaic_AN, 
		neanderthal_AC, neanderthal_AN, 
		denisovan_AC, denisovan_AN, 
		B_statistic, B_statistic_bin, 
		archaic_mask_inter, archaic_mask_union))
write_tsv(archaics_all_hub, gzfile("../../results/preprocess_1KGP_SNPs/archaics_all_hub.txt.gz"))

# # temp, remove
# archaics_all_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_hub.txt.gz"))

# add 1KG and archaic info
print("add 1KG info")
ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))
archaics_all_1KGP_archaic_hub <- ALL_1KGP_phase3_hub %>% right_join(archaics_all_hub)
write_tsv(archaics_all_1KGP_archaic_hub, gzfile("../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_hub.txt.gz"))

# convert to bed
# 	temp, remove
archaics_all_1KGP_archaic_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_hub.txt.gz"))

# 	table to bed
print("gnomAD table to bed")
archaics_all_1KGP_archaic_bed <- "../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_hub.bed"
hub_to_bed(archaics_all_1KGP_archaic_hub, archaics_all_1KGP_archaic_bed)

# # 	bed to vcf
# print("gnomAD bed to vcf")
# gnomAD_file <- "../../../final-splicing-archaic-variants/data/annotate_popgen_gnomAD/gnomAD_AF/gnomad.genomes.r2.1.1.sites.passFILTER.dropINFO.updated.vcf.gz"
# archaics_all_1KGP_archaic_gnomAD_file <- "../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_gnomAD_file.vcf.gz"
# # 	https://stackoverflow.com/questions/40845452/specify-which-shell-to-use-in-r/40845667
# system(paste("echo ", shQuote(paste("#!/bin/sh\ntabix -h ", gnomAD_file, " -R ", archaics_all_1KGP_archaic_bed, " | bcftools sort | bgzip -c >| ", archaics_all_1KGP_archaic_gnomAD_file, 
# 	"\n", "tabix -p vcf ", archaics_all_1KGP_archaic_gnomAD_file, sep="")), " >| preprocess_archaics_part1_gnomAD.sh", sep=""))
# # system("sh preprocess_archaics_part1_gnomAD.sh")
