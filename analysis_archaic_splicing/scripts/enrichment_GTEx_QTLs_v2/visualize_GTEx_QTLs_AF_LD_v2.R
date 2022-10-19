#!/bin/R

library(tidyverse)
library(data.table)
library(ggthemes)
library(MetBrewer)

# AI vs NAI comparisons
final_filtered_hub <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub.txt.gz"))
ALL_1KGP_not_AI_hub <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub.txt.gz"))

# retrieve AF and LD
final_filtered_browning_hub <- final_filtered_hub %>% 
	filter(hub_in_browning_2018_CEU | hub_in_browning_2018_IBS | hub_in_browning_2018_GBR | hub_in_browning_2018_FIN | hub_in_browning_2018_TSI | hub_in_ldvernot_akey_2016_EUR) %>% 
	dplyr::select(hub_variant_ID, EUR_AF, EUR_AC, EUR_AC_bin, EUR_hapR2tag, EUR_hapR2tag_bin) %>% 
	mutate(Description = "Introgressed")
final_filtered_ldvernot_hub<- final_filtered_hub %>% 
	filter(hub_in_browning_2018_CEU | hub_in_browning_2018_IBS | hub_in_browning_2018_GBR | hub_in_browning_2018_FIN | hub_in_browning_2018_TSI | hub_in_ldvernot_akey_2016_EUR) %>% 
	filter(hub_in_gittelman_2016 | hub_in_racimo_2017_AI) %>% 
	dplyr::select(hub_variant_ID, EUR_AF, EUR_AC, EUR_AC_bin, EUR_hapR2tag, EUR_hapR2tag_bin) %>% 
	mutate(Description = "Adaptively introgressed")
ALL_1KGP_not_AI_short_hub <- ALL_1KGP_not_AI_hub %>% 
	dplyr::select(hub_variant_ID, EUR_AF, EUR_AC, EUR_AC_bin, EUR_hapR2tag, EUR_hapR2tag_bin) %>% 
	mutate(Description = "Non-introgressed controls")
collated_AF_LD_AI_NAI_hub <- bind_rows(final_filtered_browning_hub, final_filtered_ldvernot_hub, ALL_1KGP_not_AI_short_hub)

# visualize
collated_AF_LD_AI_NAI_SFS <- collated_AF_LD_AI_NAI_hub %>% 
	group_by(Description, EUR_AF) %>% 
	summarise(EUR_AF_count = n()) %>% 
	ungroup() %>% 
	group_by(Description) %>% 
	mutate(EUR_AF_prop = EUR_AF_count/sum(EUR_AF_count, na.rm=T)) %>% 
	ungroup()

ggplot(collated_AF_LD_AI_NAI_SFS, aes(EUR_AF, EUR_AF_prop, color=Description)) + 
	geom_line() + 
	theme_base() + theme(aspect.ratio=0.5) + 
	xlab("EUR AF") + ylab("Proportion of SNPs") + 
	theme(plot.background=element_blank()) + scale_color_met_d("Juarez")
ggsave("../../results/enrichment_GTEx_QTLs_v2/figures/collated_AF_LD_AI_NAI_hub_sfs_EUR_AF_v2.pdf")

ggplot(collated_AF_LD_AI_NAI_hub, aes(EUR_AF, color=Description)) + 
	stat_ecdf() + 
	theme_base() + theme(aspect.ratio=0.5) + 
	xlab("EUR AF") + ylab("Empirical CDF") + 
	theme(plot.background=element_blank()) + scale_color_met_d("Juarez")
ggsave("../../results/enrichment_GTEx_QTLs_v2/figures/collated_AF_LD_AI_NAI_hub_ecdf_EUR_AF_v2.pdf")

# visualize
collated_AF_LD_AI_NAI_SFS <- collated_AF_LD_AI_NAI_hub %>% 
	group_by(Description, EUR_hapR2tag) %>% 
	summarise(EUR_hapR2tag_count = n()) %>% 
	ungroup() %>% 
	group_by(Description) %>% 
	mutate(EUR_hapR2tag_prop = EUR_hapR2tag_count/sum(EUR_hapR2tag_count, na.rm=T)) %>% 
	ungroup()

ggplot(collated_AF_LD_AI_NAI_SFS, aes(EUR_hapR2tag, EUR_hapR2tag_prop, color=Description)) + 
	geom_line() + 
	theme_base() + theme(aspect.ratio=0.5) + 
	xlab("EUR N_(R2>0.2, <1Mb)") + ylab("Proportion of SNPs") + 
	theme(plot.background=element_blank()) + scale_color_met_d("Juarez") + scale_x_log10()
ggsave("../../results/enrichment_GTEx_QTLs_v2/figures/collated_AF_LD_AI_NAI_hub_sfs_EUR_hapR2tag_v2.pdf")

ggplot(collated_AF_LD_AI_NAI_hub, aes(EUR_hapR2tag, color=Description)) + 
	stat_ecdf() + 
	theme_base() + theme(aspect.ratio=0.5) + 
	xlab("EUR N_(R2>0.2, <1Mb)") + ylab("Empirical CDF") + 
	theme(plot.background=element_blank()) + scale_color_met_d("Juarez") + scale_x_log10()
ggsave("../../results/enrichment_GTEx_QTLs_v2/figures/collated_AF_LD_AI_NAI_hub_ecdf_EUR_hapR2tag_v2.pdf")
