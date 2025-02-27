#!/bin/R

library(tidyverse)
library(data.table)
library(ggthemes)
library(ggpubr)

# load files
adaptive_enrichments_list <- list()
adaptive_enrichments_list[["eQTL"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/adaptively_introgressed/adaptively_introgressed_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_eQTLs_NULL-enrichment_score_table.txt"))
adaptive_enrichments_list[["sQTL"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/adaptively_introgressed/adaptively_introgressed_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_sQTLs_NULL-enrichment_score_table.txt"))

archaic_enrichments_list <- list()
archaic_enrichments_list[["eQTL"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_eQTLs_NULL-enrichment_score_table.txt"))
archaic_enrichments_list[["sQTL"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_LDtagN_bin_control_standard_sQTLs_NULL-enrichment_score_table.txt"))

# join files
adaptive_enrichments_final <- bind_rows(adaptive_enrichments_list, .id="region")
adaptive_enrichments_final$region <- factor(adaptive_enrichments_final$region, levels=c("eQTL", "sQTL"))
write_tsv(adaptive_enrichments_final, gzfile("../../results/enrichment_GTEx_QTLs_v2/direct_compare/adaptively_introgressed-region-enrichment_scores_table_pivot_LDtagN_pointerror_direct.txt.gz"))

archaic_enrichments_final <- bind_rows(archaic_enrichments_list, .id="region")
archaic_enrichments_final$region <- factor(archaic_enrichments_final$region, levels=c("eQTL", "sQTL"))
write_tsv(archaic_enrichments_final, gzfile("../../results/enrichment_GTEx_QTLs_v2/direct_compare/archaic_introgressed-region-enrichment_scores_table_pivot_LDtagN_pointerror_direct.txt.gz"))

temp_enrichments_final <- list()
temp_enrichments_final[["Introgressed"]] <- archaic_enrichments_final
temp_enrichments_final[["Adaptively Introgressed"]] <- adaptive_enrichments_final
enrichments_final <- bind_rows(temp_enrichments_final, .id="subset")

enrichments_final_summ <- enrichments_final %>% 
	group_by(subset, region) %>% 
	dplyr::summarise(
		enrichment_mean = mean(value), 
		enrichment_sd = sd(value), 
		enrichment_LCI = mean(value)-2*sd(value), 
		enrichment_UCI = mean(value)+2*sd(value), 
		enrichment_min = min(value), 
		enrichment_max = max(value), 
		enrichment_list = list(value), 
		enrichment_p_up = (sum(value<0))/(n()),
		enrichment_p_down = (sum(value>0))/(n())
	) %>% 
	ungroup() %>% 
	rowwise() %>% 
	mutate(enrichment_p = 2*min(enrichment_p_up, enrichment_p_down)) %>% 
	mutate(enrichment_p_text = ifelse(enrichment_p<0.001, "<0.001", 
		paste(signif(enrichment_p, 2), sep=""))) %>% 
	ungroup()

enrichments_final_summ_diff <- enrichments_final %>% 
	pivot_wider(names_from="subset", values_from="value") %>% 
	mutate(value_diff = `Adaptively Introgressed`-`Introgressed`) %>% 
	group_by(region) %>% 
	dplyr::summarise(
		enrichment_mean = mean(value_diff), 
		enrichment_sd = sd(value_diff), 
		enrichment_LCI = mean(value_diff)-2*sd(value_diff), 
		enrichment_UCI = mean(value_diff)+2*sd(value_diff), 
		enrichment_min = min(value_diff), 
		enrichment_max = max(value_diff), 
		enrichment_list = list(value_diff), 
		enrichment_p_up = (sum(value_diff<0))/(n()),
		enrichment_p_down = (sum(value_diff>0))/(n())
	) %>% 
	ungroup() %>% 
	rowwise() %>% 
	mutate(enrichment_p = 2*min(enrichment_p_up, enrichment_p_down)) %>% 
	mutate(enrichment_p_text = ifelse(enrichment_p<0.001, "<0.001", 
		paste(signif(enrichment_p, 2), sep=""))) %>% 
	ungroup()

# visualize
my_comparisons <- list(c("Adaptively Introgressed", "Introgressed"))
stat_test <- enrichments_final_summ_diff %>% filter(region=="eQTL") %>% 
	mutate(
		group1 = "Introgressed", group2 = "Adaptively Introgressed",
		y.position=max(filter(enrichments_final_summ, region=="eQTL")$enrichment_max)+0.03)
ggboxplot(enrichments_final %>% filter(region=="eQTL"), 
	x="subset", y="value", color="#cb96a1") + 
	# stat_compare_means(comparisons=my_comparisons, size=3) + 
	stat_pvalue_manual(stat_test, label="enrichment_p_text", size=3) + 
	geom_hline(yintercept=0, linetype="dashed", color="grey70") + 
	geom_text(data=enrichments_final_summ %>% filter(region=="eQTL"), 
		mapping=aes(x=subset, y=enrichment_min-0.03, label=enrichment_p_text), size=3) + 
	theme_base() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
	theme(legend.position = "none") + theme(aspect.ratio=1.5) + theme(plot.background=element_blank()) + 
	xlab("") + ylab("log2 odds-ratio") + theme(plot.margin = margin(0,0,0,1, "in")) + ylim(c(-0.5, 1.5))
ggsave("../../results/enrichment_GTEx_QTLs_v2/direct_compare/direct_compare_AF_LDtagN_adaptive_vs_introgressed_eQTL.pdf", device=cairo_pdf, scale=0.65)

my_comparisons <- list(c("Adaptively Introgressed", "Introgressed"))
stat_test <- enrichments_final_summ_diff %>% filter(region=="sQTL") %>% 
	mutate(
		group1 = "Introgressed", group2 = "Adaptively Introgressed",
		y.position=max(filter(enrichments_final_summ, region=="sQTL")$enrichment_max)+0.03)
ggboxplot(enrichments_final %>% filter(region=="sQTL"), 
	x="subset", y="value", color="#669e85") + 
	# stat_compare_means(comparisons=my_comparisons, size=3) + 
	stat_pvalue_manual(stat_test, label="enrichment_p_text", size=3) + 
	geom_hline(yintercept=0, linetype="dashed", color="grey70") + 
	geom_text(data=enrichments_final_summ %>% filter(region=="sQTL"), 
		mapping=aes(x=subset, y=enrichment_min-0.03, label=enrichment_p_text), size=3) + 
	theme_base() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
	theme(legend.position = "none") + theme(aspect.ratio=1.5) + theme(plot.background=element_blank()) + 
	xlab("") + ylab("log2 odds-ratio") + theme(plot.margin = margin(0,0,0,1, "in")) + ylim(c(-0.5, 1.5))
ggsave("../../results/enrichment_GTEx_QTLs_v2/direct_compare/direct_compare_AF_LDtagN_adaptive_vs_introgressed_sQTL.pdf", device=cairo_pdf, scale=0.65)
