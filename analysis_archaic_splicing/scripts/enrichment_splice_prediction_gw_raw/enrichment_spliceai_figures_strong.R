#!/bin/R

library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(Hmisc)

# load table
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs %>% filter(!(hub_variant_CHROM %in% c("X")))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs %>% filter(!is.na(spliceai_max))  # remove sites where spliceai does not make prediction

# spliceai effect bin
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs %>% 
	mutate(spliceai_sig = (spliceai_max>=0.5))

# get proportion
spliceai_proportion <- function(table, column) {
	table_prop <- table[which(table[[column]]),] %>% 
		summarise(
			spliceai_sig = sum(spliceai_sig), 
			spliceai_nonsig = n()-sum(spliceai_sig), 
			spliceai_total = n()
		)
	table_prop <- table_prop %>% 
		mutate(spliceai_prop = binconf(spliceai_sig, spliceai_total)[1]) %>% 
		mutate(spliceai_propLCI = binconf(spliceai_sig, spliceai_total)[2]) %>% 
		mutate(spliceai_propUCI = binconf(spliceai_sig, spliceai_total)[3])
	table_prop$source <- column
	table_prop <- table_prop %>% dplyr::select(source, !source)
	return(table_prop)
}

# create plotting function
enrichment_spliceai_figures <- function(list_cols, list_names, output_file, spacer, width) {
	# get proportions
	spliceai_proportion_tab <- NULL
	for (col in list_cols){
		temp <- spliceai_proportion(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs, col)
		if (is.null(spliceai_proportion_tab)) {
			spliceai_proportion_tab <- temp
		} else {
			spliceai_proportion_tab <- bind_rows(spliceai_proportion_tab, temp)
		}
	}
	spliceai_proportion_tab$description <- list_names
	spliceai_proportion_tab <- spliceai_proportion_tab %>% 
		mutate(description = factor(description, levels=list_names, labels=paste0(description, " (n = ", scales::label_comma()(spliceai_total), ")")))

	# add statistical comparison
	stat_test_pair <- pairwise_fisher_test(spliceai_proportion_tab %>% 
		dplyr::select(description, spliceai_sig, spliceai_nonsig) %>% 
		column_to_rownames("description"), detailed=TRUE)
	stat_test_pair_sig <- stat_test_pair %>% 
		mutate(p = signif(p, 2)) %>% 
		filter(p < 0.10)
		# filter((p < 0.05) | (p == min(p)))

	# does have significant pair
	if (min(stat_test_pair$p) < 0.05) {
		stat_test_pair_sig$y.position <- 
			max(spliceai_proportion_tab$spliceai_propUCI) + 
			c(1:(nrow(stat_test_pair_sig)))*spacer

		p <- ggbarplot(spliceai_proportion_tab %>% arrange(spliceai_prop), 
			x = "description", y = "spliceai_prop", fill="#386cb0") + 
			geom_errorbar(aes(description, 
				ymin=spliceai_propLCI, 
				ymax=spliceai_propUCI), stat="identity", width=width) + 
			stat_pvalue_manual(stat_test_pair_sig, label = "p", 
				tip.length = 0.01, label.size=3) + 
			theme(legend.position = "none") + ylab("Proportion SpliceAI >= 0.5") + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
			theme(aspect.ratio=1.5) + 
			theme(axis.title.x=element_blank()) + 
			ylim(c(0, max(stat_test_pair_sig$y.position)+spacer)) + 
			theme(rect = element_rect(fill = "transparent"))
	# does not have significant pair
	} else {
		p <- ggbarplot(spliceai_proportion_tab %>% arrange(spliceai_prop), 
			x = "description", y = "spliceai_prop", fill="#386cb0") + 
			geom_errorbar(aes(description, 
				ymin=spliceai_propLCI, 
				ymax=spliceai_propUCI), stat="identity", width=width) + 
			theme(legend.position = "none") + ylab("Proportion SpliceAI >= 0.5") + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
			theme(aspect.ratio=1.5) + 
			theme(axis.title.x=element_blank()) + 
			theme(rect = element_rect(fill = "transparent"))
	}
	return(p)
}

# Lineage specific

# plot - lineage specific
list_cols <- c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean")
list_names <- c("Modern", "Archaic", "Neanderthal")
output_file <- "../../results/enrichment_splice_prediction_gw_raw/Neanderthal_updated_spliceai_prop_bar_plot_stat_strong-lineage_specific_nodeni.pdf"
spacer <- 0.0003
scale <- 0.55
width <- 0.2
enrichment_spliceai_figures(list_cols, list_names, output_file, spacer, width)
ggsave(output_file, scale=scale)

# plot - lineage specific
list_cols <- c("hub_in_final_study_modern_nearly_0.5", "hub_in_final_study_modern_nearly_0.6", "hub_in_final_study_modern_nearly_0.7", "hub_in_final_study_modern_nearly_0.8", "hub_in_final_study_modern_nearly_0.9", "hub_in_final_study_modern_nearly_0.99", "hub_in_final_study_archaic", "hub_in_final_study_archaic_fixed", "hub_in_final_study_nean", "hub_in_final_study_nean_fixed")
list_names <- c("Modern (AF >= 0.5)", "Modern (AF >= 0.6)", "Modern (AF >= 0.7)", "Modern (AF >= 0.8)", "Modern (AF >= 0.9)", "Modern (AF >= 0.99)", "Archaic (AF >= 7/8)", "Archaic (AF = 8/8)", "Neanderthal (AF >= 5/6)", "Neanderthal (AF = 8/8)")
output_file <- "../../results/enrichment_splice_prediction_gw_raw/Neanderthal_updated_spliceai_prop_bar_plot_stat_strong-lineage_specific_nodeni-range.pdf"
spacer <- 0.0003
scale <- 0.8  # 0.55
width <- 0.2
enrichment_spliceai_figures(list_cols, list_names, output_file, spacer, width) + 
	theme(aspect.ratio=1.5)
ggsave(output_file, scale=scale)

# Introgressed

# plot - introgressed vs adaptively
list_cols <- c("hub_in_final_study_introgressed", "hub_in_final_study_adaptive")
list_names <- c("Introgressed", "Adaptively Introgressed")
output_file <- "../../results/enrichment_splice_prediction_gw_raw/Neanderthal_updated_spliceai_prop_bar_plot_stat_strong-introgressed.pdf"
spacer <- 0.0003
scale <- 0.65
width <- 0.2
enrichment_spliceai_figures(list_cols, list_names, output_file, spacer, width)
ggsave(output_file, scale=scale)
