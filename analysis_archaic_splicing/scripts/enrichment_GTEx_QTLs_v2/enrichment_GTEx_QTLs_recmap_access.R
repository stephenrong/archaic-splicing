#!/bin/R

library(tidyverse)
library(data.table)
library(optparse)
# library(wrapr)

# inputs
print("inputs")
option_list = list(
	make_option(
		c("-i", "--input_file"), type="character", default="../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz", help="input file of variants as hub file", metavar="character"), 
	make_option(
		c("-a", "--anno_file"), type="character", default=NULL, help="annotation file of functional variants as hub file", metavar="character"), 
	make_option(
		c("-d", "--anno_col"), type="character", default=NULL, help="col of interest in annotation file", metavar="character"), 
	make_option(
		c("-b", "--background_file"), type="character", default="../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_recmap_access_hub.txt.gz", help="background file of variants as hub file", metavar="character"), 
	make_option(
		c("-x", "--background_index"), type="character", default="../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_recmap_access_hub_EUR_AC_bin_index.rds", help="index of binned rows in background file", metavar="character"), 
	make_option(
		c("-o", "--output_prefix"), type="character", default=NULL, help="output report location", metavar="character"),
	make_option(
		c("-n", "--number_samples"), type="character", default=NULL, help="number of sampling iterations", metavar="character"),
	make_option(
		c("-s", "--set_seed"), type="character", default=NULL, help="set random seed", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input_file
anno_file <- opt$anno_file
anno_col <- opt$anno_col
background_file <- opt$background_file
background_index <- opt$background_index
output_prefix <- opt$output_prefix
number_samples <- opt$number_samples
set_seed <- opt$set_seed

# load files
print("load files")
print(input_file)
load_input_file <- as_tibble(fread(input_file))
print(anno_file)
load_anno_file <- as_tibble(fread(anno_file))

# filter relevant
print("filter relevant")
if (is.null(anno_col) | (grepl("NULL", anno_col))) {
	filter_anno_file <- load_anno_file[1:nrow(load_anno_file),][,"hub_variant_ID"]
} else {
	filter_anno_file <- load_anno_file[which(load_anno_file[[anno_col]] >= 1),][,"hub_variant_ID"]	
}

filter_input_file <- load_input_file %>% 
	dplyr::select(hub_variant_ID, EUR_AC_bin)

# load background
print("load background")
load_background_file <- as_tibble(fread(background_file))
# 	is there alternative that doesn't load into memory?

# sample matched rows
print("sample matched rows")
if (!is.null(set_seed)) {
	set.seed(set_seed)
}

load_background_index <- readRDS(background_index)

filter_input_file_tally <- filter_input_file %>% 
	group_by(EUR_AC_bin) %>% tally()

sample_input_file_index_list <- as.list(rep(NA, number_samples))
for (i in 1:number_samples) {
	print(i)

	# join and sample controls
	sample_input_file_tally <- filter_input_file_tally %>% 
		left_join(load_background_index) %>% 
		rowwise() %>% mutate(index_n = length(index)) %>% ungroup() %>% 
		rowwise() %>% mutate(sample_n = min(n, index_n)) %>% ungroup() %>%  
		rowwise() %>% mutate(sample_extra_n = n-min(n, index_n)) %>% ungroup() %>%  
		rowwise() %>% mutate(sample = list(sample(index, sample_n, replace=TRUE))) %>% ungroup() %>% 
		rowwise() %>% mutate(sample_extra = list(sample(index, sample_extra_n, replace=TRUE))) %>% ungroup()

	# # unmatched sampled controls
	# sample_input_file_tally <- sort(sample(1:nrow(load_background_file), nrow(filter_input_file), replace=FALSE))

	# get row indices
	sample_input_file_index <- sort(c(unlist(sample_input_file_tally$sample), unlist(sample_input_file_tally$sample_extra)))
	sample_input_file_index_list[[i]] <- sample_input_file_index
}

# enrichments
print("enrichments")
table_input <- table(filter_input_file$hub_variant_ID %in% filter_anno_file$hub_variant_ID)
# table_input <- table(gsub("_.*", "", filter_input_file$hub_variant_ID) %in% c("1"))

pseudocount <- function(x) {
	sum(x, 1, na.rm=T)
}

enrichment_score_list <- as.list(rep(NA, number_samples))
for (i in 1:number_samples) {
	print(i)

	# retrieve matched rows
	retrieved_input_file <- load_background_file[sample_input_file_index_list[[i]],]

	# compute enrichments
	table_background <- table(retrieved_input_file$hub_variant_ID %in% filter_anno_file$hub_variant_ID)

	# log enrichment, with pseudocount
	enrichment_score <- log2(pseudocount(table_input[2])) - log2(pseudocount(table_input[1])) - log2(pseudocount(table_background[2])) + log2(pseudocount(table_background[1]))
	enrichment_score_list[[i]] <- enrichment_score
}

# 	checking
# write_tsv(filter_input_file, paste(output_prefix, "-filter_input_file_temp.txt", sep=""))
# write_tsv(retrieved_input_file, paste(output_prefix, "-retrieved_input_file_temp.txt", sep=""))

# 	summarise
enrichment_score_table <- enframe(as.vector(unlist(enrichment_score_list)))
enrichment_score_summ <- enrichment_score_table %>% 
	summarise(
		mean_enrichment = mean(value, na.rm=T),
		sd_enrichment = sd(value, na.rm=T)
	)

# 	save
write_tsv(enrichment_score_table, paste(output_prefix, "-enrichment_score_table.txt", sep=""))
write_tsv(enrichment_score_summ, paste(output_prefix, "-enrichment_score_summ.txt", sep=""))
