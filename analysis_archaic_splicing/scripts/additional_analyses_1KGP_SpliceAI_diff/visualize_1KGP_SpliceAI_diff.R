#!/bin/R

# Compare SpliceAI scores by 

# Load packages
library(tidyverse)
library(data.table)
library(ggpubr)
library(Hmisc)

# Load tables
ALL_1KGP_phase3_common_hub_spliceai <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_common_hub_spliceai_gs_raw_dedup.txt.gz")) %>% 
	filter(!is.na(spliceai_max))

# Summarise
ALL_1KGP_phase3_common_hub_spliceai_small <- ALL_1KGP_phase3_common_hub_spliceai %>% 
	mutate(EUR = ((EUR_AF>=0.01) & (EUR_AF<=0.99))) %>% 
	mutate(EAS = ((EAS_AF>=0.01) & (EAS_AF<=0.99))) %>% 
	mutate(AFR = ((AFR_AF>=0.01) & (AFR_AF<=0.99))) %>% 
	mutate(SAS = ((SAS_AF>=0.01) & (SAS_AF<=0.99))) %>% 
	mutate(AMR = ((AMR_AF>=0.01) & (AMR_AF<=0.99))) %>% 
	mutate(spliceai_bin_ns = (spliceai_max < 0.01)) %>% 
	mutate(spliceai_bin_weak = (spliceai_max >= 0.01)&(spliceai_max < 0.2)) %>% 
	mutate(spliceai_bin_moderate = (spliceai_max >= 0.2)&(spliceai_max < 0.5)) %>% 
	mutate(spliceai_bin_strong = (spliceai_max >= 0.5)) %>% 
	dplyr::select(EUR, EAS, AFR, SAS, AMR, starts_with("spliceai_bin"))

ALL_1KGP_phase3_common_hub_spliceai_long <- ALL_1KGP_phase3_common_hub_spliceai_small %>% 
	pivot_longer(cols=c("EUR", "EAS", "AFR", "SAS", "AMR"), 
		names_to="spliceai_pop", values_to="spliceai_in_pop") %>% 
	filter(spliceai_in_pop)

ALL_1KGP_phase3_common_hub_spliceai_summ <- ALL_1KGP_phase3_common_hub_spliceai_long %>% 
	group_by(spliceai_pop) %>% 
	summarise(
		spliceai_bin_ns_n = sum(spliceai_bin_ns, na.rm=T), 
		spliceai_bin_weak_n = sum(spliceai_bin_weak, na.rm=T), 
		spliceai_bin_moderate_n = sum(spliceai_bin_moderate, na.rm=T), 
		spliceai_bin_strong_n = sum(spliceai_bin_strong, na.rm=T)
	) %>% 
	ungroup() %>% 
	rowwise() %>% 
		mutate(spliceai_bin_n = spliceai_bin_ns_n+spliceai_bin_weak_n+spliceai_bin_moderate_n+spliceai_bin_strong_n) %>% 
		mutate(spliceai_bin_ns_p = binconf(spliceai_bin_ns_n, spliceai_bin_n)[1]) %>% 
		mutate(spliceai_bin_ns_l = binconf(spliceai_bin_ns_n, spliceai_bin_n)[2]) %>% 
		mutate(spliceai_bin_ns_u = binconf(spliceai_bin_ns_n, spliceai_bin_n)[3]) %>% 
		mutate(spliceai_bin_weak_p = binconf(spliceai_bin_weak_n, spliceai_bin_n)[1]) %>% 
		mutate(spliceai_bin_weak_l = binconf(spliceai_bin_weak_n, spliceai_bin_n)[2]) %>% 
		mutate(spliceai_bin_weak_u = binconf(spliceai_bin_weak_n, spliceai_bin_n)[3]) %>% 
		mutate(spliceai_bin_moderate_p = binconf(spliceai_bin_moderate_n, spliceai_bin_n)[1]) %>% 
		mutate(spliceai_bin_moderate_l = binconf(spliceai_bin_moderate_n, spliceai_bin_n)[2]) %>% 
		mutate(spliceai_bin_moderate_u = binconf(spliceai_bin_moderate_n, spliceai_bin_n)[3]) %>% 
		mutate(spliceai_bin_strong_p = binconf(spliceai_bin_strong_n, spliceai_bin_n)[1]) %>% 
		mutate(spliceai_bin_strong_l = binconf(spliceai_bin_strong_n, spliceai_bin_n)[2]) %>% 
		mutate(spliceai_bin_strong_u = binconf(spliceai_bin_strong_n, spliceai_bin_n)[3]) %>% 
	ungroup()

# Visualize
p <- ggbarplot(ALL_1KGP_phase3_common_hub_spliceai_summ %>% arrange(spliceai_pop), 
	x = "spliceai_pop", y = "spliceai_bin_weak_p", fill="#FFC0CB") + 
	geom_errorbar(aes(spliceai_pop, 
		ymin=spliceai_bin_weak_l, 
		ymax=spliceai_bin_weak_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + 
	xlab("1KGP super-population") + ylab("Proportion SpliceAI 0.01-0.2") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)) + 
	theme(aspect.ratio=1.5) + 
	theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
p
ggsave("../../results/additional_analyses_1KGP_SpliceAI_diff/ALL_1KGP_phase3_common_hub_spliceai_summ_weak.pdf", scale=0.55)

p <- ggbarplot(ALL_1KGP_phase3_common_hub_spliceai_summ %>% arrange(spliceai_pop), 
	x = "spliceai_pop", y = "spliceai_bin_moderate_p", fill="#FFC0CB") + 
	geom_errorbar(aes(spliceai_pop, 
		ymin=spliceai_bin_moderate_l, 
		ymax=spliceai_bin_moderate_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + 
	xlab("1KGP super-population") + ylab("Proportion SpliceAI 0.2-0.5") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)) + 
	theme(aspect.ratio=1.5) + 
	theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
p
ggsave("../../results/additional_analyses_1KGP_SpliceAI_diff/ALL_1KGP_phase3_common_hub_spliceai_summ_moderate.pdf", scale=0.55)

p <- ggbarplot(ALL_1KGP_phase3_common_hub_spliceai_summ %>% arrange(spliceai_pop), 
	x = "spliceai_pop", y = "spliceai_bin_strong_p", fill="#FFC0CB") + 
	geom_errorbar(aes(spliceai_pop, 
		ymin=spliceai_bin_strong_l, 
		ymax=spliceai_bin_strong_u), stat="identity", width=0.2) + 
	theme(legend.position = "none") + 
	xlab("1KGP super-population") + ylab("Proportion SpliceAI >= 0.5") + 
	theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5)) + 
	theme(aspect.ratio=1.5) + 
	theme(axis.title.x=element_blank()) + 
	theme(rect = element_rect(fill = "transparent"))
p
ggsave("../../results/additional_analyses_1KGP_SpliceAI_diff/ALL_1KGP_phase3_common_hub_spliceai_summ_strong.pdf", scale=0.55)
