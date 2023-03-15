#!/bin/R
library(tidyverse)
library(data.table)
library(ggrepel)
library(Hmisc)
library(gtools)
library(rstatix)
library(ggpubr)

# load tables
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
mapsy_variant_table <- mapsy_variant_table %>% filter(!(exon_seqnames %in% c("chrX")))

# get proportion
mapsy_proportion <- function(table, column) {
	table_prop <- table[which(table[[column]]),] %>% 
		summarise(
			mapsy_sig = sum((mpralm.sigvar)&((mpralm.ANCDER.logFC)< -log2(3)), na.rm=T), 
			mapsy_nonsig = n()-sum((mpralm.sigvar)&((mpralm.ANCDER.logFC)< -log2(3)), na.rm=T), 
			mapsy_total = n()
		)
	table_prop <- table_prop %>% 
		mutate(mapsy_prop = binconf(mapsy_sig, mapsy_total, method="wilson")[1]) %>% 
		mutate(mapsy_propLCI = binconf(mapsy_sig, mapsy_total, method="wilson")[2]) %>% 
		mutate(mapsy_propUCI = binconf(mapsy_sig, mapsy_total, method="wilson")[3])
	table_prop$source <- column
	table_prop <- table_prop %>% dplyr::select(source, !source)
	return(table_prop)
}

# create plotting function
enrichment_mapsy_figures <- function(list_cols, list_names, output_file, spacer, width) {
	# get proportions
	mapsy_proportion_tab <- NULL
	for (col in list_cols){
		temp <- mapsy_proportion(mapsy_variant_table, col)
		if (is.null(mapsy_proportion_tab)) {
			mapsy_proportion_tab <- temp
		} else {
			mapsy_proportion_tab <- bind_rows(mapsy_proportion_tab, temp)
		}
	}
	mapsy_proportion_tab$description <- list_names
	mapsy_proportion_tab <- mapsy_proportion_tab %>% mutate(description = paste0(description, " (n = ", scales::label_comma()(mapsy_total), ")"))

	# add statistical comparison
	stat_test_pair <- pairwise_fisher_test(mapsy_proportion_tab %>% 
		dplyr::select(description, mapsy_sig, mapsy_nonsig) %>% 
		column_to_rownames("description"), detailed=TRUE)
	stat_test_pair_sig <- stat_test_pair %>% 
		mutate(p = signif(p, 2)) %>% 
		filter(p < 0.10)
		# filter((p < 0.05) | (p == min(p)))

	# does have significant pair
	if (min(stat_test_pair$p) < 0.1) {
		stat_test_pair_sig$y.position <- 
			max(mapsy_proportion_tab$mapsy_propUCI) + 
			c(1:(nrow(stat_test_pair_sig)))*spacer+spacer

		p <- ggbarplot(mapsy_proportion_tab, 
			x = "description", y = "mapsy_prop", fill="#D6604D") + 
			geom_errorbar(aes(description, 
				ymin=mapsy_propLCI, 
				ymax=mapsy_propUCI), stat="identity", width=width) + 
			stat_pvalue_manual(stat_test_pair_sig, label = "p", 
				tip.length = 0.01, label.size=3) + 
			theme(legend.position = "none") + ylab("Proportion MaPSy strong neg ESM") + 
			theme(aspect.ratio=1.5) + 
			theme(axis.title.x=element_blank()) + 
			ylim(c(0, max(stat_test_pair_sig$y.position)+spacer)) + 
			theme(rect = element_rect(fill = "transparent")) + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	# does not have significant pair
	} else {
		p <- ggbarplot(mapsy_proportion_tab, 
			x = "description", y = "mapsy_prop", fill="#D6604D") + 
			geom_errorbar(aes(description, 
				ymin=mapsy_propLCI, 
				ymax=mapsy_propUCI), stat="identity", width=width) + 
			theme(legend.position = "none") + ylab("Proportion MaPSy strong neg ESM") + 
			theme(aspect.ratio=1.5) + 
			theme(axis.title.x=element_blank()) + 
			theme(rect = element_rect(fill = "transparent")) + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
	}
	return(p)
}

# Lineage specific

# plot - lineage specific
list_cols <- c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean")
list_names <- c("Modern", "Archaic", "Neanderthal")
output_file <- "../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_prop_bar_plot_stat_strong_neg-lineage_specific_nodeni.pdf"
spacer <- 0.01
scale <- 0.55
width <- 0.2
enrichment_mapsy_figures(list_cols, list_names, output_file, spacer, width)
ggsave(output_file, scale=scale)

# plot - lineage specific
list_cols <- c("hub_in_final_study_modern", "hub_in_final_study_modern_nearly_0.99", "hub_in_final_study_archaic", "hub_in_final_study_archaic_fixed", "hub_in_final_study_nean", "hub_in_final_study_nean_fixed")
list_names <- c("Modern (AF >= 0.9)", "Modern (AF >= 0.99)", "Archaic (AF >= 7/8)", "Archaic (AF = 8/8)", "Neanderthal (AF >= 5/6)", "Neanderthal (AF = 6/6)")
output_file <- "../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_prop_bar_plot_stat_strong_neg-lineage_specific_nodeni-fixed.pdf"
spacer <- 0.01
scale <- 0.65  # 0.55
width <- 0.2
enrichment_mapsy_figures(list_cols, list_names, output_file, spacer, width)
ggsave(output_file, scale=scale)

# Introgressed

# plot - introgressed vs adaptively
list_cols <- c("hub_in_final_study_introgressed", "hub_in_final_study_adaptive")
list_names <- c("Introgressed", "Adaptively Introgressed")
output_file <- "../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_prop_bar_plot_stat_strong_neg-introgressed.pdf"
spacer <- 0.01
scale <- 0.65
width <- 0.2
enrichment_mapsy_figures(list_cols, list_names, output_file, spacer, width)
ggsave(output_file, scale=scale)
