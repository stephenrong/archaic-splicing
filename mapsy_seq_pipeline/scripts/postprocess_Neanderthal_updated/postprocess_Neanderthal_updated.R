#!/bin/R

# MaPSy analysis on Neanderthal updated library.
library(tidyverse)
library(data.table)

# load table
Neanderthal_updated_ref <- as_tibble(read.table(
	"../../data/custom_reference/neanderthal_updated_table_ref.txt.gz", sep="\t", header=T))

# load idxstats
idxstats_cleanup <- function(input) {
	input <- input %>% 
		dplyr::select(V1, V3) %>% 
		dplyr::rename(id=V1, readc=V3)
	if (all((input$readc %% 2) == 0)) {
		input <- input %>% 
			mutate(readc = readc/2) %>%  # paired read
			filter(!grepl("exon_skipped", id)) %>% 
			filter(id != "*")
	} else {
		input <- input %>% 
			mutate(readc = readc) %>%  # single read
			filter(!grepl("exon_skipped", id)) %>% 
			filter(id != "*")
	}
	return(input)
}

idxstats_dir <- "../../data/Neanderthal-HEK_293T-final_expanded_exonic/star_idxstats/"

idxstats_ls <- list()
idxstats_ls$N2A_input_a <- as_tibble(read.table(
	paste(idxstats_dir, "N2A-input-a_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2A_input_b <- as_tibble(read.table(
	paste(idxstats_dir, "N2A-input-b_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2A_input_c <- as_tibble(read.table(
	paste(idxstats_dir, "N2A-input-c_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2A_output_a <- as_tibble(read.table(
	paste(idxstats_dir, "N2A-output-a_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2A_output_b <- as_tibble(read.table(
	paste(idxstats_dir, "N2A-output-b_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2A_output_c <- as_tibble(read.table(
	paste(idxstats_dir, "N2A-output-c_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2B_input_a <- as_tibble(read.table(
	paste(idxstats_dir, "N2B-input-a_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2B_input_b <- as_tibble(read.table(
	paste(idxstats_dir, "N2B-input-b_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2B_input_c <- as_tibble(read.table(
	paste(idxstats_dir, "N2B-input-c_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2B_output_a <- as_tibble(read.table(
	paste(idxstats_dir, "N2B-output-a_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2B_output_b <- as_tibble(read.table(
	paste(idxstats_dir, "N2B-output-b_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()
idxstats_ls$N2B_output_c <- as_tibble(read.table(
	paste(idxstats_dir, "N2B-output-c_001.Aligned.unique.sortedByCoord.idxstats.txt", sep=""), sep="\t")) %>% idxstats_cleanup()

idxstats_tb <- idxstats_ls %>% 
	bind_rows(.id="NGS_library") %>% 
	pivot_wider(names_from=NGS_library, values_from=readc)

idxstats_tb_temp <- idxstats_tb %>% rowwise() %>% 
	mutate(input_readc_a = sum(N2A_input_a, N2B_input_a, na.rm=T)) %>% 
	mutate(input_readc_b = sum(N2A_input_b, N2B_input_b, na.rm=T)) %>% 
	mutate(input_readc_c = sum(N2A_input_c, N2B_input_c, na.rm=T)) %>% 
	mutate(output_readc_a = sum(N2A_output_a, N2B_output_a, na.rm=T)) %>% 
	mutate(output_readc_b = sum(N2A_output_b, N2B_output_b, na.rm=T)) %>% 
	mutate(output_readc_c = sum(N2A_output_c, N2B_output_c, na.rm=T))

idxstats_tb_final <- idxstats_tb_temp %>% rowwise() %>% 
	mutate(readc_a = sum(input_readc_a, output_readc_a, na.rm=T)) %>% 
	mutate(readc_b = sum(input_readc_b, output_readc_b, na.rm=T)) %>% 
	mutate(readc_c = sum(input_readc_c, output_readc_c, na.rm=T)) %>% 
	ungroup() %>% dplyr::select(id, readc_a, readc_b, readc_c)

# visualize read readc corrs
library(GGally)
# correlation plots
# https://stackoverflow.com/questions/45873483/
# ggpairs-plot-with-heatmap-of-correlation-values
my_fn <- function(data, mapping, method="p", use="pairwise", ...){
	x <- eval_data_col(data, mapping$x)
	y <- eval_data_col(data, mapping$y)
	corr <- cor(x, y, method=method, use=use)
	colFn <- colorRampPalette(c("#ffffff", "#ffffff", "#ffffff"), interpolate ='spline')
	# colFn <- colorRampPalette(c("#0571b0", "#f7f7f7", "#ca0020"), interpolate ='spline')
	fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
	ggally_cor(data = data, mapping = mapping, ...) + 
		theme_void() + theme(panel.background = element_rect(fill=fill))
}

# change limits
# https://stackoverflow.com/questions/53277656/
# how-to-define-facet-axis-limits-in-ggpairs-function
limitRange <- function(data, mapping, ...) { 
	ggplot(data = data, mapping = mapping, ...) + # theme_classic() + 
		geom_point(..., size=0.2, alpha=0.05) + 
		scale_y_continuous(limits = c(min(data), max(data))) +
		scale_x_continuous(limits = c(min(data), max(data))) 
}

idxstats_tb_final_input <- idxstats_tb_final %>% 
	filter(grepl("input", id)) %>% 
	filter(rowSums(.[,2:4])!=0)

names(idxstats_tb_final_input) <- 
	paste("input", names(idxstats_tb_final_input), sep="_")
names(idxstats_tb_final_input) <- 
	gsub("readc_", "", names(idxstats_tb_final_input))
idxstats_tb_final_input <- 
	log10(idxstats_tb_final_input[,2:4]+1)
ggpairs(idxstats_tb_final_input, 
	title="Input log10(read counts+1)", 
	upper = list(continuous = my_fn), 
	lower = list(continuous = limitRange))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_input_corr.png", scale=0.6)

idxstats_tb_final_output <- idxstats_tb_final %>% 
	filter(grepl("output", id)) %>% 
	filter(rowSums(.[,2:4])!=0)

names(idxstats_tb_final_output) <- 
	paste("output", names(idxstats_tb_final_output), sep="_")
names(idxstats_tb_final_output) <- 
	gsub("readc_", "", names(idxstats_tb_final_output))
idxstats_tb_final_output <- 
	log10(idxstats_tb_final_output[,2:4]+1)
ggpairs(idxstats_tb_final_output, 
	title="Output log10(read counts+1)", 
	upper = list(continuous = my_fn), 
	lower = list(continuous = limitRange))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_output_corr.png", scale=0.6)

# join tables
idxstats_temp <- function(input, pre, post) {
	names(input) <- paste(pre, names(input), post, sep="_")
	input
}

Neanderthal_updated_variant <- Neanderthal_updated_ref %>% 
	left_join(idxstats_temp(idxstats_tb_final, "ref_input", "wt")) %>% 
	left_join(idxstats_temp(idxstats_tb_final, "ref_input", "mt")) %>% 
	left_join(idxstats_temp(idxstats_tb_final, "ref_output", "wt")) %>% 
	left_join(idxstats_temp(idxstats_tb_final, "ref_output", "mt")) %>% 
	filter(exon_exon_internal)

# fix error in repeating short full
Neanderthal_updated_variant <- Neanderthal_updated_variant %>% 
	mutate(temp = !(construct_type=="short_full_exon_variant_five")) %>% 
	mutate(construct_type = ifelse(construct_type=="short_full_exon_variant_five", "short_full_exon_variant_three", construct_type)) %>% 
	group_by(source, construct_type, variant_id, exon_exon_id, exon_gene_id, exon_transcript_id, Common_id_five, Common_id_three) %>% 
	mutate(  # because the output seqs are repeated
		ref_input_readc_a_wt = sum(ref_input_readc_a_wt, na.rm=T), 
		ref_input_readc_b_wt = sum(ref_input_readc_b_wt, na.rm=T), 
		ref_input_readc_c_wt = sum(ref_input_readc_c_wt, na.rm=T), 
		ref_input_readc_a_mt = sum(ref_input_readc_a_mt, na.rm=T), 
		ref_input_readc_b_mt = sum(ref_input_readc_b_mt, na.rm=T), 
		ref_input_readc_c_mt = sum(ref_input_readc_c_mt, na.rm=T)) %>% 
	ungroup() %>% filter(temp) %>% dplyr::select(-temp)

# filter by known canonical
gencode_v32_basic_canonical_exons <- read_tsv("../../data/known_canonical/gencode_v32_lift37_basic_canonical_exons.txt.gz")
names(gencode_v32_basic_canonical_exons) <- paste("exon", names(gencode_v32_basic_canonical_exons), sep="_")

Neanderthal_updated_variant <- Neanderthal_updated_variant %>% 
	left_join(gencode_v32_basic_canonical_exons) %>% filter(exon_canonical_exon)

# remove gene read through
Neanderthal_updated_variant <- Neanderthal_updated_variant %>% 
	rowwise() %>% filter(ifelse(grepl("HLA-", exon_gene_name), TRUE, ifelse(grepl("-", exon_gene_name), ifelse(is.na(as.numeric(gsub(".*-", "", exon_gene_name))), FALSE, TRUE), TRUE))) %>% ungroup()

# remove really long exons
Neanderthal_updated_variant <- Neanderthal_updated_variant %>% 
	filter(exon_width <= 500)

# quality control
qc_readc_filters <- function(idxstats) {
	# minimum wt read counts
	idxstats <- idxstats %>% 
		rowwise() %>% 
		mutate(filter_input_readcs = 
			((ref_input_readc_a_wt>=20) & (ref_input_readc_b_wt>=20) & (ref_input_readc_c_wt>=20)) & 
			((ref_input_readc_a_mt>=20) & (ref_input_readc_b_mt>=20) & (ref_input_readc_c_mt>=20))) %>% 
		mutate(filter_output_readcs = 
			((ref_output_readc_a_mt>=1) | (ref_output_readc_b_mt>=1) | (ref_output_readc_c_mt>=1)) | 
			((ref_output_readc_a_mt>=1) | (ref_output_readc_b_mt>=1) | (ref_output_readc_c_mt>=1))) %>% 
		ungroup()
}

# qc_readc_filters <- function(idxstats) {
# 	# minimum wt read counts
# 	idxstats <- idxstats %>% 
# 		rowwise() %>% 
# 		mutate(filter_input_readcs = 
# 			(min(ref_input_readc_a_wt, ref_input_readc_b_wt, ref_input_readc_c_wt)>=20) &  # both
# 			(min(ref_input_readc_a_mt, ref_input_readc_b_mt, ref_input_readc_c_mt)>=20)) %>% 
# 		mutate(filter_output_readcs = # TRUE)
# 			(min(ref_output_readc_a_wt, ref_output_readc_b_wt, ref_output_readc_c_wt)>=20) |  # either
# 			(min(ref_output_readc_a_mt, ref_output_readc_b_mt, ref_output_readc_c_mt)>=20)) %>% ungroup()
# }

# qc_dupli_seq_filters <- function(idxstats) {
# 	# 	duplicated
# 	idxstats <- idxstats %>% 
# 		mutate(filter_unique_seqs = !duplicated(paste(
# 			ref_input_seq_mt, ref_output_seq_mt, 
# 			ref_input_seq_wt, ref_output_seq_wt, sep="")))
# }

Neanderthal_updated_variant <- 
	Neanderthal_updated_variant %>% 
	qc_readc_filters() # %>% 
	# qc_dupli_seq_filters()

# add construct type x construct id column
Neanderthal_updated_variant <- 
	Neanderthal_updated_variant %>% 
	mutate(construct_type_id = as.factor(paste(construct_type, Common_id_three, Common_id_five, sep="_")))

# naive splicing scores
naive_log2_splicing_scores_helper <- function(out_mt, in_mt, out_wt, in_wt, filter=T) {
	if (filter) {
		log2(out_mt+1) - log2(in_mt+1) - log2(out_wt+1) + log2(in_wt+1)
	} else {
		NA
	}
}

naive_log2_splicing_scores <- function(idxstats) {
	idxstats <- idxstats %>% rowwise() %>% 
		mutate(naive_log2_splicing_score_a = naive_log2_splicing_scores_helper(
			ref_output_readc_a_mt, ref_input_readc_a_mt, 
			ref_output_readc_a_wt, ref_input_readc_a_wt, 
			filter_input_readcs*filter_output_readcs)) %>% 
		mutate(naive_log2_splicing_score_b = naive_log2_splicing_scores_helper(
			ref_output_readc_b_mt, ref_input_readc_b_mt, 
			ref_output_readc_b_wt, ref_input_readc_b_wt, 
			filter_input_readcs*filter_output_readcs)) %>% 
		mutate(naive_log2_splicing_score_c = naive_log2_splicing_scores_helper(
			ref_output_readc_c_mt, ref_input_readc_c_mt, 
			ref_output_readc_c_wt, ref_input_readc_c_wt, 
			filter_input_readcs*filter_output_readcs)) %>% ungroup()
	return(idxstats)
}

Neanderthal_updated_variant <- Neanderthal_updated_variant %>% 
	naive_log2_splicing_scores()

# visualize distribution of naive splicing scores
pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_check_a.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_score_a), ylab="naive_log2_splicing_score_a")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_check_b.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_score_b), ylab="naive_log2_splicing_score_b")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_check_c.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_score_c), ylab="naive_log2_splicing_score_c")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

# visualize density of naive splicing scores
ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_score_a)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_facet_a.pdf", scale=2)

ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_score_b)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_facet_b.pdf", scale=2)

ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_score_c)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_facet_c.pdf", scale=2)

# visualize correlation of naive splicing scores
ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_score_a, naive_log2_splicing_score_b), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_scatter_a_b.pdf")

ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_score_a, naive_log2_splicing_score_c), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_scatter_a_c.pdf")

ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_score_b, naive_log2_splicing_score_c), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_score_scatter_b_c.pdf")

# naive splicing efficiency
naive_log2_splicing_efficiency_helper <- function(out_mt, in_mt, out_all, in_all, filter=T) {
	if (filter) {
		log2(out_mt+1) - log2(in_mt+1) - log2(out_all) + log2(in_all)
	} else {
		NA
	}
}

naive_log2_splicing_efficiency <- function(idxstats) {
	# wildtype
	ref_input_readc_a_wt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_input_id_wt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_input_readc_a_wt %>% sum(.+1)
	ref_input_readc_b_wt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_input_id_wt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_input_readc_b_wt %>% sum(.+1)
	ref_input_readc_c_wt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_input_id_wt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_input_readc_c_wt %>% sum(.+1)
	ref_output_readc_a_wt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_output_id_wt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_output_readc_a_wt %>% sum(.+1)
	ref_output_readc_b_wt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_output_id_wt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_output_readc_b_wt %>% sum(.+1)
	ref_output_readc_c_wt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_output_id_wt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_output_readc_c_wt %>% sum(.+1)

	idxstats <- idxstats %>% rowwise() %>% 
		mutate(naive_log2_splicing_efficiency_wt_a = naive_log2_splicing_efficiency_helper(
			ref_output_readc_a_wt, ref_input_readc_a_wt, 
			ref_output_readc_a_wt_all, ref_input_readc_a_wt_all, 
			filter_input_readcs*filter_output_readcs)) %>% 
		mutate(naive_log2_splicing_efficiency_wt_b = naive_log2_splicing_efficiency_helper(
			ref_output_readc_b_wt, ref_input_readc_b_wt, 
			ref_output_readc_b_wt_all, ref_input_readc_b_wt_all, 
			filter_input_readcs*filter_output_readcs)) %>% 
		mutate(naive_log2_splicing_efficiency_wt_c = naive_log2_splicing_efficiency_helper(
			ref_output_readc_c_wt, ref_input_readc_c_wt, 
			ref_output_readc_c_wt_all, ref_input_readc_c_wt_all, 
			filter_input_readcs*filter_output_readcs)) %>% ungroup()

	# mutant
	ref_input_readc_a_mt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_input_id_mt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_input_readc_a_mt %>% sum(.+1)
	ref_input_readc_b_mt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_input_id_mt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_input_readc_b_mt %>% sum(.+1)
	ref_input_readc_c_mt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_input_id_mt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_input_readc_c_mt %>% sum(.+1)
	ref_output_readc_a_mt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_output_id_mt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_output_readc_a_mt %>% sum(.+1)
	ref_output_readc_b_mt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_output_id_mt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_output_readc_b_mt %>% sum(.+1)
	ref_output_readc_c_mt_all <- Neanderthal_updated_variant %>% 
		filter(!duplicated(ref_output_id_mt) & filter_input_readcs*filter_output_readcs) %>% 
		.$ref_output_readc_c_mt %>% sum(.+1)

	idxstats <- idxstats %>% rowwise() %>% 
		mutate(naive_log2_splicing_efficiency_mt_a = naive_log2_splicing_efficiency_helper(
			ref_output_readc_a_mt, ref_input_readc_a_mt, 
			ref_output_readc_a_mt_all, ref_input_readc_a_mt_all, 
			filter_input_readcs*filter_output_readcs)) %>% 
		mutate(naive_log2_splicing_efficiency_mt_b = naive_log2_splicing_efficiency_helper(
			ref_output_readc_b_mt, ref_input_readc_b_mt, 
			ref_output_readc_b_mt_all, ref_input_readc_b_mt_all, 
			filter_input_readcs*filter_output_readcs)) %>% 
		mutate(naive_log2_splicing_efficiency_mt_c = naive_log2_splicing_efficiency_helper(
			ref_output_readc_c_mt, ref_input_readc_c_mt, 
			ref_output_readc_c_mt_all, ref_input_readc_c_mt_all, 
			filter_input_readcs*filter_output_readcs)) %>% ungroup()

	return(idxstats)
}

Neanderthal_updated_variant <- Neanderthal_updated_variant %>% 
	naive_log2_splicing_efficiency()

# visualize density of naive splicing efficiency
ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_efficiency_wt_a)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_facet_a.pdf", scale=2)

ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_efficiency_wt_b)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_facet_b.pdf", scale=2)

ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_efficiency_wt_c)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_facet_c.pdf", scale=2)

# visualize correlation of naive splicing efficiency
ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_efficiency_wt_a, naive_log2_splicing_efficiency_wt_b), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_scatter_a_b.pdf")

ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_efficiency_wt_a, naive_log2_splicing_efficiency_wt_c), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_scatter_a_c.pdf")

ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_efficiency_wt_b, naive_log2_splicing_efficiency_wt_c), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_scatter_b_c.pdf")

# visualize distribution of naive splicing efficiency
pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_check_a.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_efficiency_wt_a), ylab="naive_log2_splicing_efficiency_wt_a")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_check_b.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_efficiency_wt_b), ylab="naive_log2_splicing_efficiency_wt_b")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_wt_check_c.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_efficiency_wt_c), ylab="naive_log2_splicing_efficiency_wt_c")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

# visualize density of naive splicing efficiency
ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_efficiency_mt_a)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_facet_a.pdf", scale=2)

ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_efficiency_mt_b)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_facet_b.pdf", scale=2)

ggplot(Neanderthal_updated_variant) + 
	geom_density(aes(naive_log2_splicing_efficiency_mt_c)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_facet_c.pdf", scale=2)

# visualize correlation of naive splicing efficiency
ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_efficiency_mt_a, naive_log2_splicing_efficiency_mt_b), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_scatter_a_b.pdf")

ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_efficiency_mt_a, naive_log2_splicing_efficiency_mt_c), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_scatter_a_c.pdf")

ggplot(Neanderthal_updated_variant) + 
	geom_point(aes(naive_log2_splicing_efficiency_mt_b, naive_log2_splicing_efficiency_mt_c), pch=1, alpha=0.1) + 
	facet_wrap(~construct_type) + theme(aspect.ratio=1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_scatter_b_c.pdf")

# visualize distribution of naive splicing efficiency
pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_check_a.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_efficiency_mt_a), ylab="naive_log2_splicing_efficiency_mt_a")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_check_b.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_efficiency_mt_b), ylab="naive_log2_splicing_efficiency_mt_b")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

pdf("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_efficiency_mt_check_c.pdf")
plot(sort(Neanderthal_updated_variant$naive_log2_splicing_efficiency_mt_c), ylab="naive_log2_splicing_efficiency_mt_c")
abline(h=log2(2), col="red"); abline(h=-log2(2), col="red"); 
dev.off()

# save idxstats thus far
write_tsv(Neanderthal_updated_variant, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant.txt.gz"))

# evaluate bias
library(ggpubr)

Neanderthal_updated_variant_bias <- 
	Neanderthal_updated_variant %>% 
	# get family size before filtering
	group_by(
		construct_type_id,
		ref_input_id_wt, ref_output_id_wt
		) %>% 
	mutate(
		bias_family_index = cur_group_id(),
		bias_family_size = length(unique(variant_id))
		) %>% 
	ungroup() %>% 
	# filter by read counts
	filter(
		filter_input_readcs, 
		filter_output_readcs
		) %>% 
	# get bias stats after filtering
	group_by(# source,
		construct_type_id, 
		ref_input_id_wt, ref_output_id_wt, 
		bias_family_index, bias_family_size
		) %>% 
	summarise(
		ref_input_readc_a_wt = mean(ref_input_readc_a_wt, na.rm=T), 
		ref_input_readc_b_wt = mean(ref_input_readc_b_wt, na.rm=T), 
		ref_input_readc_c_wt = mean(ref_input_readc_c_wt, na.rm=T), 
		ref_output_readc_a_wt = mean(ref_output_readc_a_wt, na.rm=T), 
		ref_output_readc_b_wt = mean(ref_output_readc_b_wt, na.rm=T), 
		ref_output_readc_c_wt = mean(ref_output_readc_c_wt, na.rm=T), 
		bias_naive_log2_splicing_score_mean_a = mean(naive_log2_splicing_score_a, na.rm=T), 
		bias_naive_log2_splicing_score_mean_b = mean(naive_log2_splicing_score_b, na.rm=T), 
		bias_naive_log2_splicing_score_mean_c = mean(naive_log2_splicing_score_c, na.rm=T),
		bias_naive_log2_splicing_score_sd_a = sd(naive_log2_splicing_score_a, na.rm=T), 
		bias_naive_log2_splicing_score_sd_b = sd(naive_log2_splicing_score_b, na.rm=T), 
		bias_naive_log2_splicing_score_sd_c = sd(naive_log2_splicing_score_c, na.rm=T)
		) %>% 
	ungroup()

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_old_a.pdf")

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_old_b.pdf")

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_old_c.pdf")

# faceted version
ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_old_facet_a.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_old_facet_b.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_old_facet_c.pdf", scale=2)

# bias on a variant level
Neanderthal_updated_variant_biasvar <- 
	Neanderthal_updated_variant %>% 
	# get family size before filtering
	group_by(
		construct_type_id,
		ref_input_id_wt, ref_output_id_wt
		) %>% 
	mutate(
		bias_family_index = cur_group_id(),
		bias_family_size = length(unique(variant_id))
		) %>% 
	ungroup() %>% 
	# filter by read counts
	filter(
		filter_input_readcs, 
		filter_output_readcs
		)

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_old_a.pdf")

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_old_b.pdf")

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_old_c.pdf")

# faceted version
ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_old_facet_a.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_old_facet_b.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_old_facet_c.pdf", scale=2)

# resampling experiment
Neanderthal_updated_variant_bias_resample <- 
	Neanderthal_updated_variant_bias

Neanderthal_updated_variant_bias_null <- 
	Neanderthal_updated_variant_bias %>% 
	filter(bias_family_size == 1)

resample_helper <- function(idxstats_bias_null, bias_family_size) {
	mean(sample(c(
		idxstats_bias_null$bias_naive_log2_splicing_score_mean_a,
		idxstats_bias_null$bias_naive_log2_splicing_score_mean_b,
		idxstats_bias_null$bias_naive_log2_splicing_score_mean_c), 
		bias_family_size, replace=TRUE), na.rm=T)
}

Neanderthal_updated_variant_bias_resample <- 
	Neanderthal_updated_variant_bias_resample %>% 
	rowwise() %>% 
	mutate(bias_naive_log2_splicing_score_mean_a = 
		resample_helper(., bias_family_size)) %>% 
	mutate(bias_naive_log2_splicing_score_mean_b = 
		resample_helper(., bias_family_size)) %>% 
	mutate(bias_naive_log2_splicing_score_mean_c = 
		resample_helper(., bias_family_size)) %>% 
	ungroup()

ggscatter(Neanderthal_updated_variant_bias_resample, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_resample_a.pdf")

ggscatter(Neanderthal_updated_variant_bias_resample, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_resample_b.pdf")

ggscatter(Neanderthal_updated_variant_bias_resample, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_resample_c.pdf")

# faceted version
ggscatter(Neanderthal_updated_variant_bias_resample, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_resample_facet_a.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_bias_resample, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_resample_facet_b.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_bias_resample, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_resample_facet_c.pdf", scale=2)

# fit and apply correction
# bias_fit_a <- lm(bias_naive_log2_splicing_score_mean_a~bias_family_size,  # *construct_type_id, 
# 	data=Neanderthal_updated_variant_bias)
# bias_fit_b <- lm(bias_naive_log2_splicing_score_mean_b~bias_family_size,  # *construct_type_id, 
# 	data=Neanderthal_updated_variant_bias)
# bias_fit_c <- lm(bias_naive_log2_splicing_score_mean_c~bias_family_size,  # *construct_type_id, 
# 	data=Neanderthal_updated_variant_bias)
# bias_fit_a$coefficients[which(!grepl("bias_family_size", names(bias_fit_a$coefficients)))] <- 0
# bias_fit_b$coefficients[which(!grepl("bias_family_size", names(bias_fit_b$coefficients)))] <- 0
# bias_fit_c$coefficients[which(!grepl("bias_family_size", names(bias_fit_c$coefficients)))] <- 0

# Neanderthal_updated_variant_bias$bias_correction_a <- 1
	# 2**(-predict(bias_fit_a, newadata=Neanderthal_updated_variant_bias$bias_family_size))
# Neanderthal_updated_variant_bias$bias_correction_b <- 1
	# 2**(-predict(bias_fit_b, newadata=Neanderthal_updated_variant_bias$bias_family_size))
# Neanderthal_updated_variant_bias$bias_correction_c <- 1
	# 2**(-predict(bias_fit_c, newadata=Neanderthal_updated_variant_bias$bias_family_size))

# join bias and correct
Neanderthal_updated_variant <- 
	Neanderthal_updated_variant %>% 
	left_join(Neanderthal_updated_variant_bias)

Neanderthal_updated_variant <-
	Neanderthal_updated_variant %>% 
	mutate(
		# ref_input_readc_orig_a_wt=ref_input_readc_a_wt,
		# ref_input_readc_orig_b_wt=ref_input_readc_b_wt,
		# ref_input_readc_orig_c_wt=ref_input_readc_c_wt,
		# ref_output_readc_orig_a_wt=ref_output_readc_a_wt,
		# ref_output_readc_orig_b_wt=ref_output_readc_b_wt,
		# ref_output_readc_orig_c_wt=ref_output_readc_c_wt,
		ref_input_readc_a_wt=ref_input_readc_a_wt,
		ref_input_readc_b_wt=ref_input_readc_b_wt,
		ref_input_readc_c_wt=ref_input_readc_c_wt,
		ref_output_readc_a_wt=ref_output_readc_a_wt,  # /(bias_correction_a),
		ref_output_readc_b_wt=ref_output_readc_b_wt,  # /(bias_correction_b),
		ref_output_readc_c_wt=ref_output_readc_c_wt  # /(bias_correction_c)
	) %>% 
	naive_log2_splicing_scores()

Neanderthal_updated_variant_bias <- 
	Neanderthal_updated_variant %>% 
	# get family size before filtering
	group_by(
		construct_type_id,
		ref_input_id_wt, ref_output_id_wt
		) %>% 
	mutate(
		bias_family_index = cur_group_id(),
		bias_family_size = length(unique(variant_id))
		) %>% 
	ungroup() %>% 
	# filter by read counts
	filter(
		filter_input_readcs, 
		filter_output_readcs
		) %>% 
	# get bias stats after filtering
	group_by(# source,
		construct_type_id, 
		ref_input_id_wt, ref_output_id_wt, 
		bias_family_index, bias_family_size
		) %>% 
	summarise(
		ref_input_readc_a_wt = mean(ref_input_readc_a_wt, na.rm=T), 
		ref_input_readc_b_wt = mean(ref_input_readc_b_wt, na.rm=T), 
		ref_input_readc_c_wt = mean(ref_input_readc_c_wt, na.rm=T), 
		ref_output_readc_a_wt = mean(ref_output_readc_a_wt, na.rm=T), 
		ref_output_readc_b_wt = mean(ref_output_readc_b_wt, na.rm=T), 
		ref_output_readc_c_wt = mean(ref_output_readc_c_wt, na.rm=T), 
		bias_naive_log2_splicing_score_mean_a = mean(naive_log2_splicing_score_a, na.rm=T), 
		bias_naive_log2_splicing_score_mean_b = mean(naive_log2_splicing_score_b, na.rm=T), 
		bias_naive_log2_splicing_score_mean_c = mean(naive_log2_splicing_score_c, na.rm=T),
		bias_naive_log2_splicing_score_sd_a = sd(naive_log2_splicing_score_a, na.rm=T), 
		bias_naive_log2_splicing_score_sd_b = sd(naive_log2_splicing_score_b, na.rm=T), 
		bias_naive_log2_splicing_score_sd_c = sd(naive_log2_splicing_score_c, na.rm=T)
		) %>% 
	ungroup()

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 10) +
	stat_regline_equation(label.x = 5, label.y = 9) + 
	ylim(c(-5, 10)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_new_bias_a.pdf")

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 10) +
	stat_regline_equation(label.x = 5, label.y = 9) + 
	ylim(c(-5, 10)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_new_bias_b.pdf")

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 10) +
	stat_regline_equation(label.x = 5, label.y = 9) + 
	ylim(c(-5, 10)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_new_bias_c.pdf")

# faceted version
ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_new_facet_a.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_new_facet_b.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_bias, 
	x = "bias_family_size", 
	y = "bias_naive_log2_splicing_score_mean_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias_new_facet_c.pdf", scale=2)

write_tsv(Neanderthal_updated_variant_bias, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_bias.txt.gz"))

# bias on a variant level
Neanderthal_updated_variant_biasvar <- 
	Neanderthal_updated_variant %>% 
	# get family size before filtering
	group_by(
		construct_type_id,
		ref_input_id_wt, ref_output_id_wt
		) %>% 
	mutate(
		bias_family_index = cur_group_id(),
		bias_family_size = length(unique(variant_id))
		) %>% 
	ungroup() %>% 
	# filter by read counts
	filter(
		filter_input_readcs, 
		filter_output_readcs
		)

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_new_a.pdf")

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_new_b.pdf")

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA")
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_new_c.pdf")

# faceted version
ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_a", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_new_facet_a.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_b", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_new_facet_b.pdf", scale=2)

ggscatter(Neanderthal_updated_variant_biasvar, 
	x = "bias_family_size", 
	y = "naive_log2_splicing_score_c", 
	add = "reg.line", alpha=0.1, 
	add.params = list(color = "#377eb8", fill = "lightgray")) +
	stat_cor(label.x = 5, label.y = 5) +
	stat_regline_equation(label.x = 5, label.y = 4.5) + 
	ylim(c(-5, 5)) + 
	geom_hline(yintercept=0, color="#AAAAAA") + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_biasvar_new_facet_c.pdf", scale=2)

# mpralm splicing scores
source("../mpralm_splicing_scores/mpralm_analyze_exp.R")
mpralm_log2_splicing_scores <- function(idxstats, correlate=FALSE) {
	# only use unique pairs for mpralm
	idxstats <- idxstats %>%   # need for later!
		mutate(index = row_number())

	idxstats_temp_tb <- 
		idxstats %>% filter(filter_input_readcs, filter_output_readcs)

	if (correlate==FALSE) {
		# mpralm without correlate
		# create experiment object
		temp_experiment <- list()
		temp_experiment$input_counts_allelic <- 
			idxstats_temp_tb %>% 
				dplyr::select(
					ref_input_readc_a_mt, ref_input_readc_a_wt, 
					ref_input_readc_b_mt, ref_input_readc_b_wt, 
					ref_input_readc_c_mt, ref_input_readc_c_wt)
		temp_experiment$output_spliced_counts_allelic <- 
			idxstats_temp_tb %>% 
				dplyr::select(
					ref_output_readc_a_mt, ref_output_readc_a_wt, 
					ref_output_readc_b_mt, ref_output_readc_b_wt, 
					ref_output_readc_c_mt, ref_output_readc_c_wt)

		# calculate mpralm and fet scores
		temp_analyze_experiment <- analyze_experiment(temp_experiment)
		temp_analyze_experiment <- temp_analyze_experiment %>% 
			bind_cols(idxstats_temp_tb["index"])  # need for later!
	} else {
		# mpralm with correlate
		# create experiment object
		temp_analyze_experiment <- NULL
		for (i in unique(idxstats_temp_tb$construct_type_id)) {
			idxstats_temp_tb_i <- idxstats_temp_tb %>% 
				filter(construct_type_id==i)

			temp_experiment <- list()
			temp_experiment$input_counts_allelic <- 
				idxstats_temp_tb_i %>% 
					dplyr::select(
						ref_input_readc_a_mt, ref_input_readc_a_wt, 
						ref_input_readc_b_mt, ref_input_readc_b_wt, 
						ref_input_readc_c_mt, ref_input_readc_c_wt)
			temp_experiment$output_spliced_counts_allelic <- 
				idxstats_temp_tb_i %>% 
					dplyr::select(
						ref_output_readc_a_mt, ref_output_readc_a_wt, 
						ref_output_readc_b_mt, ref_output_readc_b_wt, 
						ref_output_readc_c_mt, ref_output_readc_c_wt)

			# calculate mpralm and fet scores
			temp_analyze_experiment_temp <- analyze_experiment(temp_experiment)
			temp_analyze_experiment_temp <- temp_analyze_experiment_temp %>% 
				bind_cols(idxstats_temp_tb_i["index"])  # need for later!
			if (is.null(temp_analyze_experiment)) {
				temp_analyze_experiment <- temp_analyze_experiment_temp
			} else {
				temp_analyze_experiment <- bind_rows(
					temp_analyze_experiment, temp_analyze_experiment_temp)
			}
		}
	}

	# combine back into original table
	idxstats_temp_tb_final <- 
		idxstats_temp_tb %>% 
		dplyr::select(index,
			ref_input_seq_wt, ref_input_seq_mt, 
			ref_output_seq_wt, ref_output_seq_mt) %>% 
		left_join(temp_analyze_experiment) %>% 
		dplyr::select(-contains("allele"))

	# return revised table
	idxstats %>% left_join(idxstats_temp_tb_final) %>% 
		dplyr::select(!c("mpralm.AveExpr", "mpralm.adj.P.Val", "mpralm.B", "chisq.adj.P.Val"))
}

Neanderthal_updated_variant_mpralm <- 
	Neanderthal_updated_variant %>% 
	mpralm_log2_splicing_scores()  # correlate=TRUE)

# get significant variants
mpralm_signif_variants <- function(idxstats_mpralm) {
	idxstats_mpralm %>% 
		mutate(mpralm.adj.Pval = p.adjust(mpralm.P.Value, method="fdr")) %>% 
		mutate(mpralm.sigvar = (mpralm.adj.Pval < 0.05)&(abs(mpralm.logFC) > log2(1.5))) %>% 
		mutate(mpralm.sigclass = ifelse(
			(mpralm.adj.Pval < 0.05)&(abs(mpralm.logFC) > log2(3)), 
			"strong", ifelse((mpralm.adj.Pval < 0.05)&(abs(mpralm.logFC) > log2(1.5)),
			"weak", "ns"))) %>%  
		mutate(chisq.adj.Pval = p.adjust(chisq.P.Value, method="fdr")) %>% 
		mutate(chisq.sigvar = (chisq.adj.Pval < 0.05)&(abs(chisq.logFC) > log2(1.5))) %>% 
		mutate(chisq.sigclass = ifelse(
			(chisq.adj.Pval < 0.05)&(abs(chisq.logFC) > log2(3)), 
			"strong", ifelse((chisq.adj.Pval < 0.05)&(abs(chisq.logFC) > log2(1.5)),
			"weak", "ns"))) %>%  
		dplyr::select(!contains(c("mpralm", "chisq")), contains("mpralm"), contains("chisq"))
}

Neanderthal_updated_variant_mpralm <- 
	Neanderthal_updated_variant_mpralm %>% 
	mpralm_signif_variants()

ggplot(Neanderthal_updated_variant_mpralm) + 
	geom_density(aes(mpralm.logFC)) + 
	facet_wrap(~(construct_type_id))
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_facet_mpralm.pdf", scale=2)

# clean up, save, and visualize
write_tsv(Neanderthal_updated_variant_mpralm, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_mpralm.txt.gz"))

# # save short variant-exon summary
# # 	old version
# # Neanderthal_updated_variant_mpralm_short <- 
# # 	Neanderthal_updated_variant_mpralm %>% 
# # 	dplyr::select(source, construct_type, construct_type_id, 
# # 		starts_with("variant"), starts_with("exon"), 
# # 		starts_with("Common"), starts_with("Primer"), 
# # 		Sub_library, Order_assign, 
# # 		starts_with("ref_input_id"), 
# # 		starts_with("ref_output_id"), 
# # 		contains("readc"), 
# # 		starts_with("filter"), starts_with("naive"), 
# # 		starts_with("bias"), 
# # 		starts_with("mpralm"), starts_with("chisq")) %>% 
# # 	unique()

# 	new version
Neanderthal_updated_variant_mpralm_short <- 
	Neanderthal_updated_variant_mpralm %>% 
	dplyr::select(source, construct_type, 
		starts_with("exon"), 
		starts_with("variant"), 
		starts_with("variable"), 
		starts_with("buffer"), 
		starts_with("Common"), 
		starts_with("Primer"), 
		starts_with("Order_debug"), 
		starts_with("Order_final"), 
		Order_assign, 
		Sub_library, 
		contains("readc"), 
		starts_with("filter"), 
		starts_with("naive_log2_splicing_score"), 
		starts_with("mpralm"), 
		starts_with("chisq")
	) %>% unique() %>% 
	dplyr::select(!(c("variable_wt_allele", "variable_mt_allele", "variable_alleles")))

write_tsv(Neanderthal_updated_variant_mpralm_short, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_mpralm_short.txt.gz"))

# output versions
Neanderthal_updated_variant_mpralm_fil_short_full <- Neanderthal_updated_variant_mpralm_short %>% 
	filter(construct_type %in% "short_full_exon_variant_three") %>% 
	mutate(var_pos_in_exon = ifelse(exon_strand == "+", variant_start-exon_start+1, (exon_end-variant_end+1))) %>% 
	mutate(var_prop_in_exon = var_pos_in_exon/exon_width)

Neanderthal_updated_variant_mpralm_fil_long_half_three <- Neanderthal_updated_variant_mpralm_short %>% 
	filter(construct_type %in% c("long_half_exon_variant_three", "short_half_exon_variant_three")) %>% 
	# filter(!(variant_id %in% Neanderthal_updated_variant_mpralm_fil_short_full$variant_id)) %>% 
	mutate(var_pos_in_exon = ifelse(exon_strand == "+", variant_start-exon_start+1, (exon_end-variant_end+1))) %>% 
	mutate(var_prop_in_exon = var_pos_in_exon/exon_width) # %>% filter(var_pos_in_exon <= floor(exon_width/2))

Neanderthal_updated_variant_mpralm_fil_long_half_five <- Neanderthal_updated_variant_mpralm_short %>% 
	filter(construct_type %in% c("long_half_exon_variant_five", "short_half_exon_variant_five")) %>% 
	# filter(!(variant_id %in% Neanderthal_updated_variant_mpralm_fil_short_full$variant_id)) %>% 
	mutate(var_pos_in_exon = ifelse(exon_strand == "+", variant_start-exon_start+1, (exon_end-variant_end+1))) %>% 
	mutate(var_prop_in_exon = var_pos_in_exon/exon_width) # %>% filter(!(var_pos_in_exon <= floor(exon_width/2)))

Neanderthal_updated_variant_mpralm_fil <- bind_rows(
	Neanderthal_updated_variant_mpralm_fil_short_full,
	Neanderthal_updated_variant_mpralm_fil_long_half_three,
	Neanderthal_updated_variant_mpralm_fil_long_half_five
)

Neanderthal_updated_variant_ind <- 
	Neanderthal_updated_variant_mpralm_fil %>% 
	filter(construct_type %in% c(
		"short_full_exon_variant_three", 
		"long_half_exon_variant_three"  # ,
		# "long_half_exon_variant_five"
	)) %>% 
	unique()
write_tsv(Neanderthal_updated_variant_ind, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_ind.txt.gz"))

Neanderthal_updated_variant_val <- 
	Neanderthal_updated_variant_mpralm_fil %>% 
	filter(construct_type %in% c(
		"short_half_exon_variant_three",
		"short_half_exon_variant_five"
	)) %>% 
	unique()
write_tsv(Neanderthal_updated_variant_val, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_val.txt.gz"))

# take min p-value
Neanderthal_updated_variant_ind_minp <- Neanderthal_updated_variant_ind %>% 
	group_by(source, 
		variant_seqnames, variant_start, variant_end, variant_width, 
		variant_strand, variant_REF, variant_ALT, variant_QUAL, 
		variant_FILTER, variant_name, variant_id, variant_type) %>% 
	arrange(mpralm.P.Value) %>% 
	mutate(mpralm.P.Value = p.adjust(mpralm.P.Value, method="bonferroni")) %>% 
	mutate(chisq.P.Value = p.adjust(chisq.P.Value, method="bonferroni")) %>% 
	filter(row_number()==1) %>% 
	ungroup() %>% 
	mpralm_signif_variants()
write_tsv(Neanderthal_updated_variant_ind, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_ind_minp.txt.gz"))

# recalculate on per variant*gene basis
Neanderthal_updated_variant_input <- 
	Neanderthal_updated_variant_ind %>% 
	dplyr::select(source, exon_strand, exon_gene_id, exon_gene_name, exon_gene_type, 
		starts_with("variant"), starts_with("ref_input_id"), 
		starts_with("ref_input_readc")) %>% unique()
Neanderthal_updated_variant_output <- 
	Neanderthal_updated_variant_ind %>% 
	dplyr::select(source, exon_strand, exon_gene_id, exon_gene_name, exon_gene_type, 
		starts_with("variant"), starts_with("ref_output_id"),
		starts_with("ref_output_readc")) %>% unique()

write_tsv(Neanderthal_updated_variant_input, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_input_sum_gene.txt.gz"))
write_tsv(Neanderthal_updated_variant_output, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_output_sum_gene.txt.gz"))

table(table(paste(Neanderthal_updated_variant_input$ref_input_id_wt, 
	Neanderthal_updated_variant_input$ref_input_id_mt, sep="_")))
table(table(paste(Neanderthal_updated_variant_output$ref_output_id_wt, 
	Neanderthal_updated_variant_output$ref_output_id_mt, sep="_")))

Neanderthal_updated_variant_input_sum <- 
	Neanderthal_updated_variant_input %>% 
	group_by(source, exon_strand, exon_gene_id, exon_gene_name, exon_gene_type, 
		variant_seqnames, variant_start, variant_end, variant_width, 
		variant_strand, variant_REF, variant_ALT, variant_QUAL, 
		variant_FILTER, variant_name, variant_id, variant_type) %>% 
	summarise(
		ref_input_readc_a_wt=sum(ref_input_readc_a_wt, na.rm=T), 
		ref_input_readc_b_wt=sum(ref_input_readc_b_wt, na.rm=T), 
		ref_input_readc_c_wt=sum(ref_input_readc_c_wt, na.rm=T),
		ref_input_readc_a_mt=sum(ref_input_readc_a_mt, na.rm=T), 
		ref_input_readc_b_mt=sum(ref_input_readc_b_mt, na.rm=T), 
		ref_input_readc_c_mt=sum(ref_input_readc_c_mt, na.rm=T)  # ,
		# ref_input_readc_orig_a_wt=sum(ref_input_readc_orig_a_wt, na.rm=T), 
		# ref_input_readc_orig_b_wt=sum(ref_input_readc_orig_b_wt, na.rm=T), 
		# ref_input_readc_orig_c_wt=sum(ref_input_readc_orig_c_wt, na.rm=T)
	) %>% ungroup()
Neanderthal_updated_variant_output_sum <- 
	Neanderthal_updated_variant_output %>% 
	group_by(source, exon_strand, exon_gene_id, exon_gene_name, exon_gene_type, 
		variant_seqnames, variant_start, variant_end, variant_width, 
		variant_strand, variant_REF, variant_ALT, variant_QUAL, 
		variant_FILTER, variant_name, variant_id, variant_type) %>% 
	summarise(
		ref_output_readc_a_wt=sum(ref_output_readc_a_wt, na.rm=T), 
		ref_output_readc_b_wt=sum(ref_output_readc_b_wt, na.rm=T), 
		ref_output_readc_c_wt=sum(ref_output_readc_c_wt, na.rm=T),
		ref_output_readc_a_mt=sum(ref_output_readc_a_mt, na.rm=T), 
		ref_output_readc_b_mt=sum(ref_output_readc_b_mt, na.rm=T), 
		ref_output_readc_c_mt=sum(ref_output_readc_c_mt, na.rm=T)  # ,
		# ref_output_readc_orig_a_wt=sum(ref_output_readc_orig_a_wt, na.rm=T), 
		# ref_output_readc_orig_b_wt=sum(ref_output_readc_orig_b_wt, na.rm=T), 
		# ref_output_readc_orig_c_wt=sum(ref_output_readc_orig_c_wt, na.rm=T)
	) %>% ungroup()

Neanderthal_updated_variant_sum <- full_join(
		Neanderthal_updated_variant_input_sum, 
		Neanderthal_updated_variant_output_sum) %>% 
	mutate(  # faked, so it passes through functions
		ref_input_seq_wt=as.character(row_number()), 
		ref_input_seq_mt=as.character(row_number()),
		ref_output_seq_wt=as.character(row_number()), 
		ref_output_seq_mt=as.character(row_number())) %>% 
	qc_readc_filters() %>% 
	# qc_dupli_seq_filters() %>% 
	naive_log2_splicing_scores() %>% 
	mpralm_log2_splicing_scores() %>% 
	mpralm_signif_variants()

Neanderthal_updated_variant_sum <- Neanderthal_updated_variant_sum %>% 
	dplyr::select(-c(index, experiment_eid, # filter_unique_seqs, 
		ref_input_seq_wt, ref_input_seq_mt, ref_output_seq_wt, ref_output_seq_mt))

write_tsv(Neanderthal_updated_variant_sum, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_sum_gene.txt.gz"))

ggplot(Neanderthal_updated_variant_sum) + 
	geom_point(aes(mpralm.logFC, -log10(mpralm.adj.Pval)), pch=1, alpha=0.1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_sum_gene.pdf")

# check retention of variants
length(unique(Neanderthal_updated_variant_mpralm_fil$variant_id))
length(unique(Neanderthal_updated_variant_sum$variant_id))

# recalculate on per variant*exon basis
Neanderthal_updated_variant_input <- 
	Neanderthal_updated_variant_ind %>% 
	dplyr::select(source, exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		starts_with("variant"), starts_with("ref_input_id"), 
		starts_with("ref_input_readc")) %>% unique()
Neanderthal_updated_variant_output <- 
	Neanderthal_updated_variant_ind %>% 
	dplyr::select(source, exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		starts_with("variant"), starts_with("ref_output_id"),
		starts_with("ref_output_readc")) %>% unique()

write_tsv(Neanderthal_updated_variant_input, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_input_sum_exon.txt.gz"))
write_tsv(Neanderthal_updated_variant_output, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_output_sum_exon.txt.gz"))

table(table(paste(Neanderthal_updated_variant_input$ref_input_id_wt, 
	Neanderthal_updated_variant_input$ref_input_id_mt, sep="_")))
table(table(paste(Neanderthal_updated_variant_output$ref_output_id_wt, 
	Neanderthal_updated_variant_output$ref_output_id_mt, sep="_")))

Neanderthal_updated_variant_input_sum <- 
	Neanderthal_updated_variant_input %>% 
	group_by(source, exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		variant_seqnames, variant_start, variant_end, variant_width, 
		variant_strand, variant_REF, variant_ALT, variant_QUAL, 
		variant_FILTER, variant_name, variant_id, variant_type) %>% 
	summarise(
		ref_input_readc_a_wt=sum(ref_input_readc_a_wt, na.rm=T), 
		ref_input_readc_b_wt=sum(ref_input_readc_b_wt, na.rm=T), 
		ref_input_readc_c_wt=sum(ref_input_readc_c_wt, na.rm=T),
		ref_input_readc_a_mt=sum(ref_input_readc_a_mt, na.rm=T), 
		ref_input_readc_b_mt=sum(ref_input_readc_b_mt, na.rm=T), 
		ref_input_readc_c_mt=sum(ref_input_readc_c_mt, na.rm=T)  # ,
		# ref_input_readc_orig_a_wt=sum(ref_input_readc_orig_a_wt, na.rm=T), 
		# ref_input_readc_orig_b_wt=sum(ref_input_readc_orig_b_wt, na.rm=T), 
		# ref_input_readc_orig_c_wt=sum(ref_input_readc_orig_c_wt, na.rm=T)
	) %>% ungroup()
Neanderthal_updated_variant_output_sum <- 
	Neanderthal_updated_variant_output %>% 
	group_by(source, exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		variant_seqnames, variant_start, variant_end, variant_width, 
		variant_strand, variant_REF, variant_ALT, variant_QUAL, 
		variant_FILTER, variant_name, variant_id, variant_type) %>% 
	summarise(
		ref_output_readc_a_wt=sum(ref_output_readc_a_wt, na.rm=T), 
		ref_output_readc_b_wt=sum(ref_output_readc_b_wt, na.rm=T), 
		ref_output_readc_c_wt=sum(ref_output_readc_c_wt, na.rm=T),
		ref_output_readc_a_mt=sum(ref_output_readc_a_mt, na.rm=T), 
		ref_output_readc_b_mt=sum(ref_output_readc_b_mt, na.rm=T), 
		ref_output_readc_c_mt=sum(ref_output_readc_c_mt, na.rm=T)  # ,
		# ref_output_readc_orig_a_wt=sum(ref_output_readc_orig_a_wt, na.rm=T), 
		# ref_output_readc_orig_b_wt=sum(ref_output_readc_orig_b_wt, na.rm=T), 
		# ref_output_readc_orig_c_wt=sum(ref_output_readc_orig_c_wt, na.rm=T)
	) %>% ungroup()

Neanderthal_updated_variant_sum <- full_join(
		Neanderthal_updated_variant_input_sum, 
		Neanderthal_updated_variant_output_sum) %>% 
	mutate(  # faked, so it passes through functions
		ref_input_seq_wt=as.character(row_number()), 
		ref_input_seq_mt=as.character(row_number()),
		ref_output_seq_wt=as.character(row_number()), 
		ref_output_seq_mt=as.character(row_number())) %>% 
	qc_readc_filters() %>% 
	# qc_dupli_seq_filters() %>% 
	naive_log2_splicing_scores() %>% 
	mpralm_log2_splicing_scores() %>% 
	mpralm_signif_variants()

Neanderthal_updated_variant_sum <- Neanderthal_updated_variant_sum %>% 
	dplyr::select(-c(index, experiment_eid, # filter_unique_seqs, 
		ref_input_seq_wt, ref_input_seq_mt, ref_output_seq_wt, ref_output_seq_mt))

write_tsv(Neanderthal_updated_variant_sum, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_sum_exon.txt.gz"))

ggplot(Neanderthal_updated_variant_sum) + 
	geom_point(aes(mpralm.logFC, -log10(mpralm.adj.Pval)), pch=1, alpha=0.1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_sum_exon.pdf")

# check retention of variants
length(unique(Neanderthal_updated_variant_mpralm_fil$variant_id))
length(unique(Neanderthal_updated_variant_sum$variant_id))

# recalculate on per variant*exon*construct basis
Neanderthal_updated_variant_input <- 
	Neanderthal_updated_variant_ind %>% 
	dplyr::select(source, construct_type, 
		exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		starts_with("variant"), starts_with("ref_input_id"), 
		starts_with("ref_input_readc")) %>% unique()
Neanderthal_updated_variant_output <- 
	Neanderthal_updated_variant_ind %>% 
	dplyr::select(source, construct_type,
		exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		starts_with("variant"), starts_with("ref_output_id"),
		starts_with("ref_output_readc")) %>% unique()

write_tsv(Neanderthal_updated_variant_input, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_input_sum_exon_construct.txt.gz"))
write_tsv(Neanderthal_updated_variant_output, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_output_sum_exon_construct.txt.gz"))

table(table(paste(Neanderthal_updated_variant_input$ref_input_id_wt, 
	Neanderthal_updated_variant_input$ref_input_id_mt, sep="_")))
table(table(paste(Neanderthal_updated_variant_output$ref_output_id_wt, 
	Neanderthal_updated_variant_output$ref_output_id_mt, sep="_")))

Neanderthal_updated_variant_input_sum <- 
	Neanderthal_updated_variant_input %>% 
	group_by(source, construct_type, 
		exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		variant_seqnames, variant_start, variant_end, variant_width, 
		variant_strand, variant_REF, variant_ALT, variant_QUAL, 
		variant_FILTER, variant_name, variant_id, variant_type) %>% 
	summarise(
		ref_input_readc_a_wt=sum(ref_input_readc_a_wt, na.rm=T), 
		ref_input_readc_b_wt=sum(ref_input_readc_b_wt, na.rm=T), 
		ref_input_readc_c_wt=sum(ref_input_readc_c_wt, na.rm=T),
		ref_input_readc_a_mt=sum(ref_input_readc_a_mt, na.rm=T), 
		ref_input_readc_b_mt=sum(ref_input_readc_b_mt, na.rm=T), 
		ref_input_readc_c_mt=sum(ref_input_readc_c_mt, na.rm=T)  # ,
		# ref_input_readc_orig_a_wt=sum(ref_input_readc_orig_a_wt, na.rm=T), 
		# ref_input_readc_orig_b_wt=sum(ref_input_readc_orig_b_wt, na.rm=T), 
		# ref_input_readc_orig_c_wt=sum(ref_input_readc_orig_c_wt, na.rm=T)
	) %>% ungroup()
Neanderthal_updated_variant_output_sum <- 
	Neanderthal_updated_variant_output %>% 
	group_by(source, construct_type, 
		exon_seqnames, exon_start, exon_end, 
		exon_width, exon_strand, exon_exon_id, 
		exon_exon_name, exon_exon_rank, exon_exon_internal, 
		exon_gene_id, exon_gene_type, exon_gene_name, 
		exon_transcript_id, exon_transcript_type, exon_transcript_name, 
		variant_seqnames, variant_start, variant_end, variant_width, 
		variant_strand, variant_REF, variant_ALT, variant_QUAL, 
		variant_FILTER, variant_name, variant_id, variant_type) %>% 
	summarise(
		ref_output_readc_a_wt=sum(ref_output_readc_a_wt, na.rm=T), 
		ref_output_readc_b_wt=sum(ref_output_readc_b_wt, na.rm=T), 
		ref_output_readc_c_wt=sum(ref_output_readc_c_wt, na.rm=T),
		ref_output_readc_a_mt=sum(ref_output_readc_a_mt, na.rm=T), 
		ref_output_readc_b_mt=sum(ref_output_readc_b_mt, na.rm=T), 
		ref_output_readc_c_mt=sum(ref_output_readc_c_mt, na.rm=T)  # ,
		# ref_output_readc_orig_a_wt=sum(ref_output_readc_orig_a_wt, na.rm=T), 
		# ref_output_readc_orig_b_wt=sum(ref_output_readc_orig_b_wt, na.rm=T), 
		# ref_output_readc_orig_c_wt=sum(ref_output_readc_orig_c_wt, na.rm=T)
	) %>% ungroup()

Neanderthal_updated_variant_sum <- full_join(
		Neanderthal_updated_variant_input_sum, 
		Neanderthal_updated_variant_output_sum) %>% 
	mutate(  # faked, so it passes through functions
		ref_input_seq_wt=as.character(row_number()), 
		ref_input_seq_mt=as.character(row_number()),
		ref_output_seq_wt=as.character(row_number()), 
		ref_output_seq_mt=as.character(row_number())) %>% 
	qc_readc_filters() %>% 
	# qc_dupli_seq_filters() %>% 
	naive_log2_splicing_scores() %>% 
	mpralm_log2_splicing_scores() %>% 
	mpralm_signif_variants()

Neanderthal_updated_variant_sum <- Neanderthal_updated_variant_sum %>% 
	dplyr::select(-c(index, experiment_eid, # filter_unique_seqs, 
		ref_input_seq_wt, ref_input_seq_mt, ref_output_seq_wt, ref_output_seq_mt))

write_tsv(Neanderthal_updated_variant_sum, 
	gzfile("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_sum_exon_construct.txt.gz"))

ggplot(Neanderthal_updated_variant_sum) + 
	geom_point(aes(mpralm.logFC, -log10(mpralm.adj.Pval)), pch=1, alpha=0.1)
ggsave("../../results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_sum_exon_construct.pdf")

# check retention of variants
length(unique(Neanderthal_updated_variant_mpralm_fil$variant_id))
length(unique(Neanderthal_updated_variant_sum$variant_id))
