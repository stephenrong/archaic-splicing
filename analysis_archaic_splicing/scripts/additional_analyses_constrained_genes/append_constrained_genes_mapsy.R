#!/bin/R

library(tidyverse)
library(data.table)
library(ggpubr)
library(Hmisc)

mapsy_updated <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
gnomad_v2_constraint <- as_tibble(fread("../../data/gnomAD_v2_constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz"))

gnomad_v2_constraint_small <- gnomad_v2_constraint %>% 
	dplyr::select(gene, oe_lof_upper_bin, oe_lof_upper_bin_6) %>% 
	dplyr::rename(exon_gene_name = gene, LOEUF_bin_10=oe_lof_upper_bin, LOEUF_bin_6=oe_lof_upper_bin_6) %>% 
	mutate(LOEUF_bin_10 = as.factor(LOEUF_bin_10), LOEUF_bin_6 = as.factor(LOEUF_bin_6))

mapsy_updated_bin <- mapsy_updated %>% 
	mutate(
		mapsy_bin_ns = (mpralm.sigclass=="ns"),
		mapsy_bin_weak_neg = (mpralm.sigclass=="weak")&(mpralm.ANCDER.logFC<0),
		mapsy_bin_weak_pos = (mpralm.sigclass=="weak")&(mpralm.ANCDER.logFC>0),
		mapsy_bin_strong_neg = (mpralm.sigclass=="strong")&(mpralm.ANCDER.logFC<0),
		mapsy_bin_strong_pos = (mpralm.sigclass=="strong")&(mpralm.ANCDER.logFC>0)
	)

mapsy_updated_constraint <- mapsy_updated_bin %>% 
	left_join(gnomad_v2_constraint_small)

mapsy_updated_constraint_summ <- mapsy_updated_constraint %>% 
	group_by(LOEUF_bin_6) %>% 
	summarise(
		mapsy_bin_ns_n = sum(mapsy_bin_ns, na.rm=T), 
		mapsy_bin_weak_neg_n = sum(mapsy_bin_weak_neg, na.rm=T), 
		mapsy_bin_weak_pos_n = sum(mapsy_bin_weak_pos, na.rm=T), 
		mapsy_bin_strong_neg_n = sum(mapsy_bin_strong_neg, na.rm=T),
		mapsy_bin_strong_pos_n = sum(mapsy_bin_strong_pos, na.rm=T)
	) %>% 
	mutate(
		mapsy_bin_n = mapsy_bin_ns_n + mapsy_bin_weak_neg_n + mapsy_bin_weak_pos_n + mapsy_bin_strong_neg_n + mapsy_bin_strong_pos_n
	) %>% 
	rowwise() %>% 
	mutate(
		mapsy_bin_ns_p = binconf(mapsy_bin_ns_n, mapsy_bin_n)[[1]], 
		mapsy_bin_ns_l = binconf(mapsy_bin_ns_n, mapsy_bin_n)[[2]], 
		mapsy_bin_ns_u = binconf(mapsy_bin_ns_n, mapsy_bin_n)[[3]], 

		mapsy_bin_weak_neg_p = binconf(mapsy_bin_weak_neg_n, mapsy_bin_n)[[1]], 
		mapsy_bin_weak_neg_l = binconf(mapsy_bin_weak_neg_n, mapsy_bin_n)[[2]], 
		mapsy_bin_weak_neg_u = binconf(mapsy_bin_weak_neg_n, mapsy_bin_n)[[3]], 

		mapsy_bin_weak_pos_p = binconf(mapsy_bin_weak_pos_n, mapsy_bin_n)[[1]], 
		mapsy_bin_weak_pos_l = binconf(mapsy_bin_weak_pos_n, mapsy_bin_n)[[2]], 
		mapsy_bin_weak_pos_u = binconf(mapsy_bin_weak_pos_n, mapsy_bin_n)[[3]], 

		mapsy_bin_strong_neg_p = binconf(mapsy_bin_strong_neg_n, mapsy_bin_n)[[1]], 
		mapsy_bin_strong_neg_l = binconf(mapsy_bin_strong_neg_n, mapsy_bin_n)[[2]], 
		mapsy_bin_strong_neg_u = binconf(mapsy_bin_strong_neg_n, mapsy_bin_n)[[3]],

		mapsy_bin_strong_pos_p = binconf(mapsy_bin_strong_pos_n, mapsy_bin_n)[[1]], 
		mapsy_bin_strong_pos_l = binconf(mapsy_bin_strong_pos_n, mapsy_bin_n)[[2]], 
		mapsy_bin_strong_pos_u = binconf(mapsy_bin_strong_pos_n, mapsy_bin_n)[[3]]
	) %>% 
	ungroup()

ggbarplot(mapsy_updated_constraint_summ %>% filter(!is.na(LOEUF_bin_6)), 
	x = "LOEUF_bin_6", y = "mapsy_bin_ns_p", fill="LOEUF_bin_6") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_6, 
		ymin=mapsy_bin_ns_l, 
		ymax=mapsy_bin_ns_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF sextile") + ylab("Proportion MaPSy non-ESM") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/mapsy_updated_constrained_genes_mapsy_ns.pdf", scale=0.45)

ggbarplot(mapsy_updated_constraint_summ %>% filter(!is.na(LOEUF_bin_6)), 
	x = "LOEUF_bin_6", y = "mapsy_bin_weak_neg_p", fill="LOEUF_bin_6") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_6, 
		ymin=mapsy_bin_weak_neg_l, 
		ymax=mapsy_bin_weak_neg_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF sextile") + ylab("Proportion MaPSY weak neg ESM") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/mapsy_updated_constrained_genes_mapsy_weak_neg.pdf", scale=0.45)

ggbarplot(mapsy_updated_constraint_summ %>% filter(!is.na(LOEUF_bin_6)), 
	x = "LOEUF_bin_6", y = "mapsy_bin_weak_pos_p", fill="LOEUF_bin_6") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_6, 
		ymin=mapsy_bin_weak_pos_l, 
		ymax=mapsy_bin_weak_pos_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF sextile") + ylab("Proportion MaPSy weak pos ESM") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/mapsy_updated_constrained_genes_mapsy_weak_pos.pdf", scale=0.45)

ggbarplot(mapsy_updated_constraint_summ %>% filter(!is.na(LOEUF_bin_6)), 
	x = "LOEUF_bin_6", y = "mapsy_bin_strong_neg_p", fill="LOEUF_bin_6") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_6, 
		ymin=mapsy_bin_strong_neg_l, 
		ymax=mapsy_bin_strong_neg_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF sextile") + ylab("Proportion MaPSy strong neg ESM") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/mapsy_updated_constrained_genes_mapsy_strong_neg.pdf", scale=0.45)

ggbarplot(mapsy_updated_constraint_summ %>% filter(!is.na(LOEUF_bin_6)), 
	x = "LOEUF_bin_6", y = "mapsy_bin_strong_pos_p", fill="LOEUF_bin_6") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_6, 
		ymin=mapsy_bin_strong_pos_l, 
		ymax=mapsy_bin_strong_pos_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF sextile") + ylab("Proportion MaPSy strong pos ESM") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/mapsy_updated_constrained_genes_mapsy_strong_pos.pdf", scale=0.45)
