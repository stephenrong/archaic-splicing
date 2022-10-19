#!/bin/R

library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(Hmisc)

# load table
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup_intlin.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- as_tibble(fread("../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep[,c(4,130:142)]
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- right_join(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep)

# vep_nonsense effect bin
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep %>% 
	mutate(vep_nonsense_sig = (stop_gained))

# get proportion
vep_nonsense_proportion <- function(table, column) {
	table_prop <- table[which(table[[column]]),] %>% 
		filter(!is.na(vep_nonsense_sig)) %>% 
		summarise(
			vep_nonsense_sig = sum(vep_nonsense_sig), 
			vep_nonsense_nonsig = n()-sum(vep_nonsense_sig), 
			vep_nonsense_total = n()
		)
	table_prop <- table_prop %>% 
		mutate(vep_nonsense_prop = binconf(vep_nonsense_sig, vep_nonsense_total)[1]) %>% 
		mutate(vep_nonsense_propLCI = binconf(vep_nonsense_sig, vep_nonsense_total)[2]) %>% 
		mutate(vep_nonsense_propUCI = binconf(vep_nonsense_sig, vep_nonsense_total)[3])
	table_prop$source <- column
	table_prop <- table_prop %>% dplyr::select(source, !source)
	return(table_prop)
}

# create plotting function
enrichment_vep_nonsense_figures <- function(list_cols, list_names, output_file, spacer, width) {
	# get proportions
	vep_nonsense_proportion_tab <- NULL
	for (col in list_cols){
		temp <- vep_nonsense_proportion(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep, col)
		if (is.null(vep_nonsense_proportion_tab)) {
			vep_nonsense_proportion_tab <- temp
		} else {
			vep_nonsense_proportion_tab <- bind_rows(vep_nonsense_proportion_tab, temp)
		}
	}
	vep_nonsense_proportion_tab$description <- list_names

	vep_nonsense_proportion_tab <- vep_nonsense_proportion_tab %>% 
		mutate(description = factor(description, levels=list_names))

	# add statistical comparison
	stat_test_pair <- pairwise_fisher_test(vep_nonsense_proportion_tab %>% 
		dplyr::select(description, vep_nonsense_sig, vep_nonsense_nonsig) %>% 
		column_to_rownames("description"), detailed=TRUE)
	stat_test_pair_sig <- stat_test_pair %>% 
		filter(p < 0.10)
		# filter((p < 0.05) | (p == min(p)))

	# does have significant pair
	if (min(stat_test_pair$p) < 0.1) {
		stat_test_pair_sig$y.position <- 
			max(vep_nonsense_proportion_tab$vep_nonsense_propUCI) + 
			c(1:(nrow(stat_test_pair_sig)))*spacer

		p <- ggbarplot(vep_nonsense_proportion_tab %>% arrange(vep_nonsense_prop), 
			x = "description", y = "vep_nonsense_prop", fill="#2ca25f") + 
			geom_errorbar(aes(description, 
				ymin=vep_nonsense_propLCI, 
				ymax=vep_nonsense_propUCI), stat="identity", width=width) + 
			stat_pvalue_manual(stat_test_pair_sig, label = "p", 
				tip.length = 0.01, label.size=3) + 
			theme(legend.position = "none") + ylab("Proportion nonsense variant") + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
			theme(aspect.ratio=1.5) + 
			theme(axis.title.x=element_blank()) + 
			ylim(c(0, max(stat_test_pair_sig$y.position)+spacer)) + 
			theme(rect = element_rect(fill = "transparent"))
	# does not have significant pair
	} else {
		p <- ggbarplot(vep_nonsense_proportion_tab %>% arrange(vep_nonsense_prop), 
			x = "description", y = "vep_nonsense_prop", fill="#2ca25f") + 
			geom_errorbar(aes(description, 
				ymin=vep_nonsense_propLCI, 
				ymax=vep_nonsense_propUCI), stat="identity", width=width) + 
			theme(legend.position = "none") + ylab("Proportion nonsense variant") + 
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
output_file <- "../../results/additional_analyses_enrichment_vep_prediction/Neanderthal_updated_vep_nonsense_prop_bar_plot_stat_nonsense-lineage_specific_nodeni.pdf"
spacer <- 0.0005
scale <- 0.55
width <- 0.2
enrichment_vep_nonsense_figures(list_cols, list_names, output_file, spacer, width)
ggsave(output_file, scale=scale)
