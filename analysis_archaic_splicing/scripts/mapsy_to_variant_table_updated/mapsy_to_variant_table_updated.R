#!/bin/R
library(tidyverse)
library(data.table)
library(ggrepel)
library(Hmisc)
library(gtools)
library(rstatix)
library(ggpubr)
library(viridis)

# merge variant data
variant_table <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.txt.gz"))

# mapsy data
mapsy_table <- as_tibble(fread("../../../final-postmapsy-archaic-EndToEnd2/results/postprocess_Neanderthal_updated/Neanderthal_updated_variant_sum_exon_construct.txt.gz"))

# join tables
mapsy_variant_table <- mapsy_table %>% 
	mutate(hub_variant_ID=gsub("chr", "", variant_id)) %>% 
	mutate(hub_variant_ID=gsub(":", "_", hub_variant_ID)) %>% 
	left_join(variant_table)

# filter exon width
mapsy_variant_table <- mapsy_variant_table %>% 
	filter(exon_width <= 500)

# anc/der polarize
# mapsy_variant_table <- mapsy_variant_table %>% 
# 	mutate(mpralm.ANCDER.logFC = ifelse(
# 		hub_variant_ALT==hub_variant_DER, mpralm.logFC, -mpralm.logFC))
mapsy_variant_table <- mapsy_variant_table %>% 
	mutate(mpralm.ANCDER.logFC = 
		ifelse(hub_variant_ALT==hub_variant_DER, mpralm.logFC, 
			ifelse(hub_variant_ALT==hub_variant_ANC, -mpralm.logFC, NA)))

mapsy_variant_table <- mapsy_variant_table %>% 
	dplyr::select(!starts_with("hub_"), starts_with("hub_"))

# remove variants with ambiguous anc/der polarization
mapsy_variant_table <- mapsy_variant_table %>% 
	filter(!is.na(mpralm.ANCDER.logFC))

# save tables
write_tsv(mapsy_variant_table, 
	gzfile("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))

# volcano plots
plot_save_volcano <- function(head_n, table, hub_in, title_pdf, out_pdf) {
	temp <- table %>%  # https://stackoverflow.com/questions/34219912/how-to-use-a-variable-in-dplyrfilter
		filter(!!rlang::sym(hub_in))  # https://stackoverflow.com/questions/49786597/r-dplyr-filter-with-a-dynamic-variable-name
	# head_sum <- nrow(temp %>% filter(mpralm.sigvar)6

	if(!is.na(head_n)) {
		temp_text_left <- temp %>% 
			filter(is.finite(mpralm.ANCDER.logFC)) %>% 
			arrange(mpralm.ANCDER.logFC) %>% 
			head(head_n)
		temp_text_right <- temp %>% 
			filter(is.finite(mpralm.ANCDER.logFC)) %>% 
			arrange(mpralm.ANCDER.logFC) %>% 
			tail(head_n)
		# temp_text_up <- temp %>% 
		# 	filter(mpralm.sigvar) %>% 
		# 	arrange(-log10(mpralm.adj.Pval)) %>% 
		# 	tail(head_n)
		temp_text <- bind_rows(temp_text_left, temp_text_right) %>% unique()  # , temp_text_up) %>% unique()
	} else {
		temp_text <- temp %>% 
			filter(mpralm.sigvar) %>% unique()
	}

	ggplot(temp) + 
		ggtitle(title_pdf) + theme_classic() + 
		geom_vline(xintercept = 0, color="grey") + 
		geom_hline(yintercept = -log10(0.05), color="grey") + 
		geom_point(aes(
			x=mpralm.ANCDER.logFC, 
			y=-log10(mpralm.adj.Pval), 
			color=mpralm.sigclass), alpha=0.8, size=1) + 
		geom_text_repel(data=temp_text, mapping=aes(
			x=mpralm.ANCDER.logFC, 
			y=-log10(mpralm.adj.Pval), 
			label=exon_gene_name), size=3, min.segment.length=0, max.overlaps=60) + 
		# scale_color_viridis_d() + 
		scale_color_manual(values=c("#dcdcdc", "#d6604d", "#92c5de")) + 
		xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + 
		ylab(expression(paste("-log"[10], " FDR"))) + 
		theme(aspect.ratio=1) + theme(legend.position = "none") + 
		theme(plot.title = element_text(hjust = 0.5)) + 
		xlim(c(
			min(table$mpralm.ANCDER.logFC, na.rm=T)-0.5, 
			max(table$mpralm.ANCDER.logFC, na.rm=T)+0.5)) + 
		ylim(c(0, 
			max(-log10(table$mpralm.adj.Pval), na.rm=T)+1))
}

plot_save_volcano(5, 
	mapsy_variant_table,
	"hub_in_final_study_modern", 
	"Modern specific")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_final_study_modern.pdf", scale=0.4)

plot_save_volcano(5, 
	mapsy_variant_table,
	"hub_in_final_study_archaic", 
	"Archaic specific")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_final_study_archaic.pdf", scale=0.4)

plot_save_volcano(5, 
	mapsy_variant_table,
	"hub_in_final_study_nean", 
	"Neanderthal specific")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_final_study_nean.pdf", scale=0.4)

plot_save_volcano(5, 
	mapsy_variant_table,
	"hub_in_final_study_deni", 
	"Denisovan specific")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_final_study_deni.pdf", scale=0.4)

plot_save_volcano(5, 
	mapsy_variant_table,
	"hub_in_final_study_introgressed", 
	"Introgressed")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_final_study_introgressed.pdf", scale=0.4)

plot_save_volcano(NA, 
	mapsy_variant_table,
	"hub_in_final_study_adaptive", 
	"Adaptively introgressed")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_final_study_adaptive.pdf", scale=0.4)

plot_save_volcano(5, 
	mapsy_variant_table,
	"hub_in_ldvernot_akey_2016", 
	"Introgressed (Vernot 2016)")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_vernot_akey_2016.pdf", scale=0.4)

plot_save_volcano(5, 
	mapsy_variant_table,
	"hub_in_browning_2018", 
	"Introgressed (Browning 2018)")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_browning_2018.pdf", scale=0.4)

plot_save_volcano(NA, 
	mapsy_variant_table,
	"hub_in_gittelman_2016", 
	"Adaptively (Gittelman 2016)")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_gittelman_2016.pdf", scale=0.4)

plot_save_volcano(NA, 
	mapsy_variant_table,
	"hub_in_racimo_2017_AI", 
	"Adaptively (Racimo 2017)")
ggsave("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_volcano_hub_in_racimo_2017_AI.pdf", scale=0.4)
