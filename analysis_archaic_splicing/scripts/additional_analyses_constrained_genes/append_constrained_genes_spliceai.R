#!/bin/R

library(tidyverse)
library(data.table)
library(ggpubr)
library(Hmisc)

final_v2_spliceai <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
gnomad_v2_constraint <- as_tibble(fread("../../data/gnomAD_v2_constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz"))

gnomad_v2_constraint_small <- gnomad_v2_constraint %>% 
	dplyr::select(gene, oe_lof_upper_bin, oe_lof_upper_bin_6) %>% 
	dplyr::rename(SYMBOL = gene, LOEUF_bin_10=oe_lof_upper_bin, LOEUF_bin_6=oe_lof_upper_bin_6) %>% 
	mutate(LOEUF_bin_10 = as.factor(LOEUF_bin_10), LOEUF_bin_6 = as.factor(LOEUF_bin_6))

final_v2_spliceai_bin <- final_v2_spliceai %>% 
	mutate(
		spliceai_bin_ns = (spliceai_max < 0.01),
		spliceai_bin_weak = ((spliceai_max >= 0.01)&(spliceai_max < 0.2)), 
		spliceai_bin_moderate = ((spliceai_max >= 0.2)&(spliceai_max < 0.5)), 
		spliceai_bin_strong = (spliceai_max >= 0.5)
	)

final_v2_spliceai_constraint <- final_v2_spliceai_bin %>% 
	left_join(gnomad_v2_constraint_small)

final_v2_spliceai_constraint_summ <- final_v2_spliceai_constraint %>% 
	group_by(LOEUF_bin_10) %>% 
	summarise(
		spliceai_bin_ns_n = sum(spliceai_bin_ns, na.rm=T), 
		spliceai_bin_weak_n = sum(spliceai_bin_weak, na.rm=T), 
		spliceai_bin_moderate_n = sum(spliceai_bin_moderate, na.rm=T), 
		spliceai_bin_strong_n = sum(spliceai_bin_strong, na.rm=T)
	) %>% 
	mutate(
		spliceai_bin_n = spliceai_bin_ns_n + spliceai_bin_weak_n + spliceai_bin_moderate_n + spliceai_bin_strong_n
	) %>% 
	rowwise() %>% 
	mutate(
		spliceai_bin_ns_p = binconf(spliceai_bin_ns_n, spliceai_bin_n)[[1]], 
		spliceai_bin_ns_l = binconf(spliceai_bin_ns_n, spliceai_bin_n)[[2]], 
		spliceai_bin_ns_u = binconf(spliceai_bin_ns_n, spliceai_bin_n)[[3]], 

		spliceai_bin_weak_p = binconf(spliceai_bin_weak_n, spliceai_bin_n)[[1]], 
		spliceai_bin_weak_l = binconf(spliceai_bin_weak_n, spliceai_bin_n)[[2]], 
		spliceai_bin_weak_u = binconf(spliceai_bin_weak_n, spliceai_bin_n)[[3]], 

		spliceai_bin_moderate_p = binconf(spliceai_bin_moderate_n, spliceai_bin_n)[[1]], 
		spliceai_bin_moderate_l = binconf(spliceai_bin_moderate_n, spliceai_bin_n)[[2]], 
		spliceai_bin_moderate_u = binconf(spliceai_bin_moderate_n, spliceai_bin_n)[[3]], 

		spliceai_bin_strong_p = binconf(spliceai_bin_strong_n, spliceai_bin_n)[[1]], 
		spliceai_bin_strong_l = binconf(spliceai_bin_strong_n, spliceai_bin_n)[[2]], 
		spliceai_bin_strong_u = binconf(spliceai_bin_strong_n, spliceai_bin_n)[[3]]
	) %>% 
	ungroup()

ggbarplot(final_v2_spliceai_constraint_summ %>% filter(!is.na(LOEUF_bin_10)), 
	x = "LOEUF_bin_10", y = "spliceai_bin_ns_p", fill="LOEUF_bin_10") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_10, 
		ymin=spliceai_bin_ns_l, 
		ymax=spliceai_bin_ns_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF decile") + ylab("Proportion SpliceAI < 0.01") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/final_v2_spliceai_constrained_genes_spliceai_ns.pdf", scale=0.5)

ggbarplot(final_v2_spliceai_constraint_summ %>% filter(!is.na(LOEUF_bin_10)), 
	x = "LOEUF_bin_10", y = "spliceai_bin_weak_p", fill="LOEUF_bin_10") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_10, 
		ymin=spliceai_bin_weak_l, 
		ymax=spliceai_bin_weak_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF decile") + ylab("Proportion SpliceAI 0.01-0.2") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/final_v2_spliceai_constrained_genes_spliceai_weak.pdf", scale=0.5)

ggbarplot(final_v2_spliceai_constraint_summ %>% filter(!is.na(LOEUF_bin_10)), 
	x = "LOEUF_bin_10", y = "spliceai_bin_moderate_p", fill="LOEUF_bin_10") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_10, 
		ymin=spliceai_bin_moderate_l, 
		ymax=spliceai_bin_moderate_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF decile") + ylab("Proportion SpliceAI 0.2-0.5") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/final_v2_spliceai_constrained_genes_spliceai_moderate.pdf", scale=0.5)

ggbarplot(final_v2_spliceai_constraint_summ %>% filter(!is.na(LOEUF_bin_10)), 
	x = "LOEUF_bin_10", y = "spliceai_bin_strong_p", fill="LOEUF_bin_10") + 
	scale_fill_viridis_d() + 
	geom_errorbar(aes(LOEUF_bin_10, 
		ymin=spliceai_bin_strong_l, 
		ymax=spliceai_bin_strong_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + xlab("LOEUF decile") + ylab("Proportion SpliceAI >= 0.5") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) + 
	theme(aspect.ratio=1.5) + 
	# theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
ggsave("../../results/additional_analyses_constrained_genes/final_v2_spliceai_constrained_genes_spliceai_strong.pdf", scale=0.5)
