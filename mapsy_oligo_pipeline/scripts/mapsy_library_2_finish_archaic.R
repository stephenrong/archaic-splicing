#!/bin/R
library(tidyverse)
source("mapsy_library_2_params.R")

# Add primer and common sequences to create final oligos.

# Preamble
# load commons
final_commons_assignments <- read_tsv(
	"../data/commons/commons_verified_exonWT_master.txt")

final_commons_assignments_five <- final_commons_assignments %>% 
	dplyr::select(Common_id, Sequence_5_master_55_15) %>% 
	mutate(Common_id = gsub("Common_", "Common_five_", Common_id)) %>% 
	dplyr::rename(Common_id_five = Common_id, 
		Common_five = Sequence_5_master_55_15)

final_commons_assignments_three <- final_commons_assignments %>% 
	dplyr::select(Common_id, Sequence_3_master_55_15) %>% 
	mutate(Common_id = gsub("Common_", "Common_three_", Common_id)) %>% 
	dplyr::rename(Common_id_three = Common_id, 
		Common_three = Sequence_3_master_55_15)

# load primers
primer_table <- read_tsv(
	"../data/primers/primers_to_sequences_library2.txt")
primer_to_sequence_assignments <- primer_table %>% 
	filter(!is.na(Primer_pair_new)) %>% 
	dplyr::select(Primer_pair_new, Primer_side, Sequence_forward_5_to_3) %>% 
	pivot_wider(names_from=Primer_side, values_from=Sequence_forward_5_to_3)

# check primers
table(primer_to_sequence_assignments$Primer_5_prime)
table(primer_to_sequence_assignments$Primer_3_prime)

Biostrings::stringDist(primer_to_sequence_assignments$Primer_5_prime)
Biostrings::stringDist(primer_to_sequence_assignments$Primer_3_prime)

# load sublibraries
sublibrary_to_primer_assignments <- read_tsv(
	"../data/primers/sublibrary_to_primers_library2.txt")

# downsample short halfs
# down_sample_short_half <- 1000
down_sample_short_half <- Inf

# New Neanderthal library exonic
Neanderthal_updated_library_exon <- NULL
Neanderthal_updated_library_exon$short_full_exon_variant_three <- read_tsv(
	"../results/mapsy_libraries/library_2/Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2_short_full_exon_variant_three.txt.gz")
Neanderthal_updated_library_exon$short_full_exon_variant_five <- read_tsv(
	"../results/mapsy_libraries/library_2/Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2_short_full_exon_variant_five.txt.gz")
Neanderthal_updated_library_exon$short_half_exon_variant_three <- read_tsv(
	"../results/mapsy_libraries/library_2/Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2_short_half_exon_variant_three.txt.gz")
Neanderthal_updated_library_exon$short_half_exon_variant_five <- read_tsv(
	"../results/mapsy_libraries/library_2/Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2_short_half_exon_variant_five.txt.gz")
Neanderthal_updated_library_exon$long_half_exon_variant_three <- read_tsv(
	"../results/mapsy_libraries/library_2/Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2_long_half_exon_variant_three.txt.gz")
Neanderthal_updated_library_exon$long_half_exon_variant_five <- read_tsv(
	"../results/mapsy_libraries/library_2/Neanderthal_GRCh37_GENCODE_32_merge_variants_library_2_long_half_exon_variant_five.txt.gz")

# down sample
set.seed(20200626)
Neanderthal_updated_library_exon$short_half_exon_variant_three <- 
	Neanderthal_updated_library_exon$short_half_exon_variant_three[sort(sample(
		c(1:nrow(Neanderthal_updated_library_exon$short_half_exon_variant_three)), 
		min(down_sample_short_half, nrow(Neanderthal_updated_library_exon$short_half_exon_variant_three)))),]
Neanderthal_updated_library_exon$short_half_exon_variant_five <- 
	Neanderthal_updated_library_exon$short_half_exon_variant_five[sort(sample(
		c(1:nrow(Neanderthal_updated_library_exon$short_half_exon_variant_five)), 
		min(down_sample_short_half, nrow(Neanderthal_updated_library_exon$short_half_exon_variant_five)))),]

# bind tables
Neanderthal_updated_library_exon_variable <- bind_rows(
	Neanderthal_updated_library_exon, .id="construct_type")

# add commons
# 	create commons mapping
Neanderthal_updated_library_exon_common_mapping <- tibble(
	construct_type = c(
		"short_half_exon_variant_five", 
		"long_half_exon_variant_five",
		"short_half_exon_variant_three", 
		"long_half_exon_variant_three"
	),
	Common_id_three = as.character(c(
		"Common_three_05", "Common_three_05",
		NA, NA
	)),
	Common_id_five = as.character(c(
		NA, NA,
		"Common_five_05", "Common_five_05"
	)))

# 	join to main table
Neanderthal_updated_library_exon_common <- left_join(
	Neanderthal_updated_library_exon_variable, 
	Neanderthal_updated_library_exon_common_mapping)

# 	pull sequences
Neanderthal_updated_library_exon_common <- Neanderthal_updated_library_exon_common %>% 
	left_join(final_commons_assignments_five) %>% 
	left_join(final_commons_assignments_three) %>% 
	mutate(Common_five = ifelse(is.na(Common_five), "", Common_five)) %>% 
	mutate(Common_three = ifelse(is.na(Common_three), "", Common_three))

# manually assign primers
# 	check
Neanderthal_updated_library_exon_common %>% group_by(construct_type) %>% 
	summarise(length(final_wt_seq) + length(final_mt_seq))
Neanderthal_updated_library_exon_common %>% group_by(construct_type) %>% 
	summarise(length(unique(c(toupper(final_wt_seq), toupper(final_mt_seq)))))

# 	assign
Neanderthal_updated_library_exon_primer <- Neanderthal_updated_library_exon_common %>% 
	mutate(Sub_library = NA) %>% 
	mutate(Sub_library = ifelse(construct_type %in% 
		c("short_full_exon_variant_three"), 
		"Sub_2.01", Sub_library)) %>% 
	mutate(Sub_library = ifelse(construct_type %in% 
		c("short_half_exon_variant_three"), 
		"Sub_2.02", Sub_library)) %>% 
	mutate(Sub_library = ifelse((construct_type %in% 
		c("long_half_exon_variant_three")) & 
		(variant_start %% 2 == 0), 
		"Sub_2.03a", Sub_library)) %>% 
	mutate(Sub_library = ifelse((construct_type %in% 
		c("long_half_exon_variant_three")) & 
		(variant_start %% 2 == 1), 
		"Sub_2.03b", Sub_library)) %>% 
	mutate(Sub_library = ifelse(construct_type %in% 
		c("short_full_exon_variant_five"), 
		"Sub_2.04", Sub_library)) %>% 
	mutate(Sub_library = ifelse(construct_type %in% 
		c("short_half_exon_variant_five"), 
		"Sub_2.05", Sub_library)) %>% 
	mutate(Sub_library = ifelse((construct_type %in% 
		c("long_half_exon_variant_five")) & 
		(variant_start %% 2 == 0), 
		"Sub_2.06a", Sub_library)) %>% 
	mutate(Sub_library = ifelse((construct_type %in% 
		c("long_half_exon_variant_five")) & 
		(variant_start %% 2 == 1), 
		"Sub_2.06b", Sub_library))

Neanderthal_updated_library_exon_primer <- Neanderthal_updated_library_exon_primer %>% 
	left_join(sublibrary_to_primer_assignments) %>% 
	left_join(primer_to_sequence_assignments)

# 	check
Neanderthal_updated_library_exon_primer %>% group_by(Sub_library) %>% 
	summarise(length(final_wt_seq) + length(final_mt_seq))
Neanderthal_updated_library_exon_primer %>% group_by(Sub_library) %>% 
	summarise(length(unique(c(toupper(final_wt_seq), toupper(final_mt_seq)))))

# get final ids
Neanderthal_updated_library_exon_final <- Neanderthal_updated_library_exon_primer %>% 
	mutate(Order_temp_id = paste(
		description, design, variant_id, variable_alleles, 
		exon_gene_name, exon_transcript_name, exon_exon_rank, exon_strand,
		Common_id_three, Common_id_five, Sub_library, Order_assign, Primer_pair_new, sep="|")) %>% 
	mutate(Order_final_id_wt = paste(Order_temp_id, "wt", sep="|")) %>% 
	mutate(Order_final_id_mt = paste(Order_temp_id, "mt", sep="|")) %>% 
	dplyr::select(-Order_temp_id)

# get final sequences
Neanderthal_updated_library_exon_final <- Neanderthal_updated_library_exon_final %>% 
	mutate(Order_debug_sequence_wt = paste(
		Primer_5_prime, chartr("ACGTacgt", "acgtACGT", Common_three), final_wt_seq, 
		chartr("ACGTacgt", "acgtACGT", Common_five), Primer_3_prime, sep="_")) %>% 
	mutate(Order_debug_sequence_mt = paste(
		Primer_5_prime, chartr("ACGTacgt", "acgtACGT", Common_three), final_mt_seq, 
		chartr("ACGTacgt", "acgtACGT", Common_five), Primer_3_prime, sep="_")) %>% 
	mutate(Order_final_sequence_wt = gsub("_", "", Order_debug_sequence_wt)) %>% 
	mutate(Order_final_sequence_mt = gsub("_", "", Order_debug_sequence_mt)) %>% 
	mutate(Order_length_wt = nchar(Order_final_sequence_wt)) %>% 
	mutate(Order_length_mt = nchar(Order_final_sequence_mt))

# 	check sequence characters
table(grepl("N|n", Neanderthal_updated_library_exon_final$Order_final_sequence_wt))
table(grepl("N|n", Neanderthal_updated_library_exon_final$Order_final_sequence_mt))

# 	check sequence lengths
table(Neanderthal_updated_library_exon_final$Order_length_wt)
table(Neanderthal_updated_library_exon_final$Order_length_mt)

# save above table
write_tsv(Neanderthal_updated_library_exon_final, 
	gzfile("../results/mapsy_orders/library_2/Neanderthal_updated_library_exon_final.txt.gz"))

# 	check
Neanderthal_updated_library_exon_final %>% group_by(Sub_library) %>% 
	summarise(length(Order_final_sequence_wt) + length(Order_final_sequence_mt))
Neanderthal_updated_library_exon_final %>% group_by(Sub_library) %>% 
	summarise(length(unique(c(toupper(Order_final_sequence_wt), toupper(Order_final_sequence_mt)))))

# 	check
table(table(toupper(Neanderthal_updated_library_exon_final$Order_final_sequence_wt)))
table(table(toupper(Neanderthal_updated_library_exon_final$Order_final_sequence_mt)))

# save order format
Neanderthal_updated_library_exon_order_wt <- Neanderthal_updated_library_exon_final %>% 
	dplyr::select(Order_final_id_wt, Order_final_sequence_wt)
Neanderthal_updated_library_exon_order_mt <- Neanderthal_updated_library_exon_final %>% 
	dplyr::select(Order_final_id_mt, Order_final_sequence_mt)
names(Neanderthal_updated_library_exon_order_wt) <- c("id", "seq")
names(Neanderthal_updated_library_exon_order_mt) <- c("id", "seq")

# 	with duplicates
Neanderthal_updated_library_exon_order_temp <- rbind(
	Neanderthal_updated_library_exon_order_wt, Neanderthal_updated_library_exon_order_mt) %>% 
	arrange(id)
write_tsv(Neanderthal_updated_library_exon_order_temp, 
	gzfile("../results/mapsy_orders/library_2/Neanderthal_updated_library_exon_order_temp.txt.gz"))

# 	without duplicates
Neanderthal_updated_library_exon_order_final <- Neanderthal_updated_library_exon_order_temp[
	which(!duplicated(toupper(Neanderthal_updated_library_exon_order_temp$seq))),] %>% 
	arrange(id)
write_tsv(Neanderthal_updated_library_exon_order_final, 
	gzfile("../results/mapsy_orders/library_2/Neanderthal_updated_library_exon_order_final.txt.gz"))

# 	check
nrow(Neanderthal_updated_library_exon_order_temp)
nrow(Neanderthal_updated_library_exon_order_final)

# 	check
Neanderthal_updated_library_exon_diagnostic <- Neanderthal_updated_library_exon_final %>% 
	group_by(Sub_library) %>% 
	summarise(
		variants_wt_mt = length(Order_final_id_wt), 
		length_unique_wt = length(unique(
			toupper(Order_final_sequence_wt))),
		length_unique_mt = length(unique(
			toupper(Order_final_sequence_mt))),
		length_unique_wt_mt = length(unique(c(
			toupper(Order_final_sequence_wt), 
			toupper(Order_final_sequence_mt)))))
Neanderthal_updated_library_exon_diagnostic
write_tsv(Neanderthal_updated_library_exon_diagnostic, 
	"../results/mapsy_orders/library_2/Neanderthal_updated_library_exon_diagnostic_sub_library.txt")

Neanderthal_updated_library_exon_diagnostic <- Neanderthal_updated_library_exon_final %>% 
	group_by(construct_type) %>% 
	summarise(
		variants_wt_mt = length(Order_final_id_wt), 
		length_unique_wt = length(unique(
			toupper(Order_final_sequence_wt))),
		length_unique_mt = length(unique(
			toupper(Order_final_sequence_mt))),
		length_unique_wt_mt = length(unique(c(
			toupper(Order_final_sequence_wt), 
			toupper(Order_final_sequence_mt)))))
Neanderthal_updated_library_exon_diagnostic
write_tsv(Neanderthal_updated_library_exon_diagnostic, 
	"../results/mapsy_orders/library_2/Neanderthal_updated_library_exon_diagnostic_construct_type.txt")

# check size of sub libraries
Neanderthal_updated_library_exon_primer %>% group_by(Sub_library, construct_type) %>% 
	summarise(count_seqs = length(unique(c(toupper(final_wt_seq), toupper(final_mt_seq)))))

# check size of primer pairs
Neanderthal_updated_library_exon_primer %>% group_by(Primer_pair_new, construct_type) %>% 
	summarise(count_seqs = length(unique(c(toupper(final_wt_seq), toupper(final_mt_seq)))))

# check size of orders
Neanderthal_updated_library_exon_primer %>% group_by(Order_assign, Primer_pair_new) %>% 
	summarise(count_seqs = length(unique(c(toupper(final_wt_seq), toupper(final_mt_seq)))))

# check for not assigned primer pair
sum(is.na(Neanderthal_updated_library_exon_primer$Primer_pair_new))

# check for not assigned primer 5'
sum(is.na(Neanderthal_updated_library_exon_primer$Primer_5_prime))

# check for not assigned primer 3'
sum(is.na(Neanderthal_updated_library_exon_primer$Primer_3_prime))

# check length of wt sequence
sum(Neanderthal_updated_library_exon_final$Order_length_wt != 230)

# check length of mt sequence
sum(Neanderthal_updated_library_exon_final$Order_length_mt != 230)

# check wt sequence characters
sum(grepl("N|n", Neanderthal_updated_library_exon_final$Order_final_sequence_wt))

# check mt sequence characters
sum(grepl("N|n", Neanderthal_updated_library_exon_final$Order_final_sequence_wt))
