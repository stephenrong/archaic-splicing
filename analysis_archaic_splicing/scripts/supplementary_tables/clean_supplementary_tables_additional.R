#!/bin/R

library(tidyverse)
library(data.table)

# load tables
supplement_mapsy_results_merged <- read_tsv("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_MaPSy_scores_main_experiment_SpliceAI_DeltaEI_VEP.txt.gz")
supplement_variant_sets_final <- read_tsv("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_Human_variant_sets_and_SpliceAI.txt.gz")

supplement_mapsy_results_merged_finemap_ukb <- read_tsv("../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_ukb.txt.gz")
supplement_variant_sets_final_finemap_ukb <- read_tsv("../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_ukb.txt.gz")
supplement_mapsy_results_merged_finemap_ukb <- supplement_mapsy_results_merged_finemap_ukb %>% 
	filter(hub_in_final_study_modern|hub_in_final_study_archaic|hub_in_final_study_nean|hub_in_final_study_deni|hub_in_final_study_introgressed|hub_in_final_study_adaptive)
supplement_variant_sets_final_finemap_ukb <- supplement_variant_sets_final_finemap_ukb %>% 
	filter(hub_in_final_study_modern|hub_in_final_study_archaic|hub_in_final_study_nean|hub_in_final_study_deni|hub_in_final_study_introgressed|hub_in_final_study_adaptive)
supplement_mapsy_results_merged_finemap_ukb <- supplement_mapsy_results_merged_finemap_ukb[,c(112,209:236)]
supplement_variant_sets_final_finemap_ukb <- supplement_variant_sets_final_finemap_ukb[,c(4,158:185)]
names(supplement_mapsy_results_merged_finemap_ukb) <- c("main_variant_ID", paste("FineMapped", names(supplement_mapsy_results_merged_finemap_ukb)[2:29], sep="_"))
names(supplement_variant_sets_final_finemap_ukb) <- c("main_variant_ID", paste("FineMapped", names(supplement_variant_sets_final_finemap_ukb)[2:29], sep="_"))
supplement_mapsy_results_merged_finemap_ukb <- supplement_mapsy_results_merged_finemap_ukb %>% 
	mutate(`FineMapped_cohort`="UKB") %>% 
	dplyr::select(main_variant_ID, `FineMapped_cohort`, `FineMapped_method`, `FineMapped_description`, `FineMapped_pip`, `FineMapped_cs_id`) %>% 
	dplyr::rename(`FineMapped_trait`=`FineMapped_description`)
supplement_variant_sets_final_finemap_ukb <- supplement_variant_sets_final_finemap_ukb %>% 
	mutate(`FineMapped_cohort`="UKB") %>% 
	dplyr::select(main_variant_ID, `FineMapped_cohort`, `FineMapped_method`, `FineMapped_description`, `FineMapped_pip`, `FineMapped_cs_id`) %>% 
	dplyr::rename(`FineMapped_trait`=`FineMapped_description`)

supplement_mapsy_results_merged_finemap_bbj_FINEMAP <- read_tsv("../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_bbj_FINEMAP.txt.gz") %>% 
	mutate(cohort = "BBJ", method="SUSIE")
supplement_variant_sets_final_finemap_bbj_FINEMAP <- read_tsv("../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_bbj_FINEMAP.txt.gz") %>% 
	mutate(cohort = "BBJ", method="SUSIE")
supplement_mapsy_results_merged_finemap_bbj_SUSIE <- read_tsv("../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_bbj_SUSIE.txt.gz") %>% 
	mutate(cohort = "BBJ", method="SUSIE") %>% dplyr::select(-starts_with("alpha"), -starts_with("lbf"))
supplement_variant_sets_final_finemap_bbj_SUSIE <- read_tsv("../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_bbj_SUSIE.txt.gz") %>% 
	mutate(cohort = "BBJ", method="SUSIE") %>% dplyr::select(-starts_with("alpha"), -starts_with("lbf"))

supplement_mapsy_results_merged_finemap_bbj_FINEMAP <- supplement_mapsy_results_merged_finemap_bbj_FINEMAP %>% 
	filter(hub_in_final_study_modern|hub_in_final_study_archaic|hub_in_final_study_nean|hub_in_final_study_deni|hub_in_final_study_introgressed|hub_in_final_study_adaptive)
supplement_variant_sets_final_finemap_bbj_FINEMAP <- supplement_variant_sets_final_finemap_bbj_FINEMAP %>% 
	filter(hub_in_final_study_modern|hub_in_final_study_archaic|hub_in_final_study_nean|hub_in_final_study_deni|hub_in_final_study_introgressed|hub_in_final_study_adaptive)
supplement_mapsy_results_merged_finemap_bbj_FINEMAP <- supplement_mapsy_results_merged_finemap_bbj_FINEMAP[,c(112,209:235)]
supplement_variant_sets_final_finemap_bbj_FINEMAP <- supplement_variant_sets_final_finemap_bbj_FINEMAP[,c(4,158:184)]
names(supplement_mapsy_results_merged_finemap_bbj_FINEMAP) <- c("main_variant_ID", paste("FineMapped", names(supplement_mapsy_results_merged_finemap_bbj_FINEMAP)[2:28], sep="_"))
names(supplement_variant_sets_final_finemap_bbj_FINEMAP) <- c("main_variant_ID", paste("FineMapped", names(supplement_variant_sets_final_finemap_bbj_FINEMAP)[2:28], sep="_"))

supplement_mapsy_results_merged_finemap_bbj_SUSIE <- supplement_mapsy_results_merged_finemap_bbj_SUSIE %>% 
	filter(hub_in_final_study_modern|hub_in_final_study_archaic|hub_in_final_study_nean|hub_in_final_study_deni|hub_in_final_study_introgressed|hub_in_final_study_adaptive)
supplement_variant_sets_final_finemap_bbj_SUSIE <- supplement_variant_sets_final_finemap_bbj_SUSIE %>% 
	filter(hub_in_final_study_modern|hub_in_final_study_archaic|hub_in_final_study_nean|hub_in_final_study_deni|hub_in_final_study_introgressed|hub_in_final_study_adaptive)
supplement_mapsy_results_merged_finemap_bbj_SUSIE <- supplement_mapsy_results_merged_finemap_bbj_SUSIE[,c(112,209:235)]
supplement_variant_sets_final_finemap_bbj_SUSIE <- supplement_variant_sets_final_finemap_bbj_SUSIE[,c(4,158:184)]
names(supplement_mapsy_results_merged_finemap_bbj_SUSIE) <- c("main_variant_ID", paste("FineMapped", names(supplement_mapsy_results_merged_finemap_bbj_SUSIE)[2:28], sep="_"))
names(supplement_variant_sets_final_finemap_bbj_SUSIE) <- c("main_variant_ID", paste("FineMapped", names(supplement_variant_sets_final_finemap_bbj_SUSIE)[2:28], sep="_"))

supplement_mapsy_results_merged_finemap_bbj_FINEMAP <- supplement_mapsy_results_merged_finemap_bbj_FINEMAP %>% 
	mutate(`FineMapped_cohort`="BBJ", `FineMapped_method`="FINEMAP") %>% 
	dplyr::select(main_variant_ID, `FineMapped_cohort`, `FineMapped_method`, `FineMapped_trait`, `FineMapped_pip`, `FineMapped_cs_id`)
supplement_variant_sets_final_finemap_bbj_FINEMAP <- supplement_variant_sets_final_finemap_bbj_FINEMAP %>% 
	mutate(`FineMapped_cohort`="BBJ", `FineMapped_method`="FINEMAP") %>% 
	dplyr::select(main_variant_ID, `FineMapped_cohort`, `FineMapped_method`, `FineMapped_trait`, `FineMapped_pip`, `FineMapped_cs_id`)

supplement_mapsy_results_merged_finemap_bbj_SUSIE <- supplement_mapsy_results_merged_finemap_bbj_SUSIE %>% 
	mutate(`FineMapped_cohort`="BBJ", `FineMapped_method`="SUSIE") %>% 
	dplyr::select(main_variant_ID, `FineMapped_cohort`, `FineMapped_method`, `FineMapped_trait`, `FineMapped_pip`, `FineMapped_cs_id`)
supplement_variant_sets_final_finemap_bbj_SUSIE <- supplement_variant_sets_final_finemap_bbj_SUSIE %>% 
	mutate(`FineMapped_cohort`="BBJ", `FineMapped_method`="SUSIE") %>% 
	dplyr::select(main_variant_ID, `FineMapped_cohort`, `FineMapped_method`, `FineMapped_trait`, `FineMapped_pip`, `FineMapped_cs_id`)

supplement_mapsy_results_merged_finemap_bbj <- bind_rows(supplement_mapsy_results_merged_finemap_bbj_FINEMAP, supplement_mapsy_results_merged_finemap_bbj_SUSIE)
supplement_variant_sets_final_finemap_bbj <- bind_rows(supplement_variant_sets_final_finemap_bbj_FINEMAP, supplement_variant_sets_final_finemap_bbj_SUSIE)

# top hits mapsy
main_in_final_study_modern_top <- supplement_mapsy_results_merged %>% 
	filter(mpralm.sigvar, main_in_final_study_modern) %>% arrange(abs(mpralm.ANCDER.logFC)) %>% tail(10) %>% 
	mutate(main_in_final_study_modern_top10 = TRUE)
main_in_final_study_archaic_top <- supplement_mapsy_results_merged %>% 
	filter(mpralm.sigvar, main_in_final_study_archaic) %>% arrange(abs(mpralm.ANCDER.logFC)) %>% tail(10) %>% 
	mutate(main_in_final_study_archaic_top10 = TRUE)
main_in_final_study_nean_top <- supplement_mapsy_results_merged %>% 
	filter(mpralm.sigvar, main_in_final_study_nean) %>% arrange(abs(mpralm.ANCDER.logFC)) %>% tail(10) %>% 
	mutate(main_in_final_study_nean_top10 = TRUE)
main_in_final_study_deni_top <- supplement_mapsy_results_merged %>% 
	filter(mpralm.sigvar, main_in_final_study_deni) %>% arrange(abs(mpralm.ANCDER.logFC)) %>% tail(10) %>% 
	mutate(main_in_final_study_deni_top10 = TRUE)
main_in_final_study_introgressed_top <- supplement_mapsy_results_merged %>% 
	filter(mpralm.sigvar, main_in_final_study_introgressed) %>% arrange(abs(mpralm.ANCDER.logFC)) %>% tail(10) %>% 
	mutate(main_in_final_study_introgressed_top10 = TRUE)
main_in_final_study_adaptive_top <- supplement_mapsy_results_merged %>% 
	filter(mpralm.sigvar, main_in_final_study_adaptive) %>% arrange(abs(mpralm.ANCDER.logFC)) %>% tail(Inf) %>% 
	mutate(main_in_final_study_adaptive_topall = TRUE)

supplement_mapsy_results_merged_top <- bind_rows(
	main_in_final_study_modern_top,
	main_in_final_study_archaic_top,
	main_in_final_study_nean_top,
	main_in_final_study_deni_top,
	main_in_final_study_introgressed_top,
	main_in_final_study_adaptive_top
)

supplement_mapsy_results_merged_dup <- supplement_mapsy_results_merged_top %>% 
	mutate(temp = duplicated(main_variant_ID)) %>% filter(temp)
supplement_mapsy_results_merged_top <- supplement_mapsy_results_merged_top %>% 
	mutate(main_in_final_study_top_multiple_sets = (main_variant_ID %in% supplement_mapsy_results_merged_dup$main_variant_ID))
write_tsv(supplement_mapsy_results_merged_top, "../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_MaPSy_scores_top_hits_per_variant_set.txt.gz")

# lineage-specific hits
supplement_mapsy_results_merged_lineage <- supplement_mapsy_results_merged %>% 
	filter(main_in_final_study_modern|main_in_final_study_archaic|main_in_final_study_nean|main_in_final_study_deni|main_in_final_study_introgressed|main_in_final_study_adaptive) %>% 
	filter(!main_in_final_study_introgressed, !main_in_final_study_adaptive) %>% 
	filter(is.na(`1KGP_AF`), is.na(gnomAD_AF)) %>% 
	filter(mpralm.sigvar, SpliceAI_max>=0.1) %>% 
	arrange(mpralm.ANCDER.logFC) %>% 
	mutate(`MaPSy_strong_+_SpliceAI>=0.2` = ifelse(mpralm.sigclass=="strong" & SpliceAI_max>=0.2, TRUE, FALSE)) %>% 
	arrange(`MaPSy_strong_+_SpliceAI>=0.2`)
write_tsv(supplement_mapsy_results_merged_lineage, "../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_MaPSy_scores_fixed_lineage_specific.txt.gz")

# fine map mapsy
supplement_mapsy_results_merged_finemap_ukb <- supplement_mapsy_results_merged %>% 
	right_join(bind_rows(supplement_mapsy_results_merged_finemap_ukb, supplement_mapsy_results_merged_finemap_bbj)) %>% arrange(desc(FineMapped_cohort), FineMapped_trait)
write_tsv(supplement_mapsy_results_merged_finemap_ukb, "../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_MaPSy_scores_overlap_fine_map.txt.gz")

# fine map spliceai
supplement_variant_sets_final_finemap_ukb <- supplement_variant_sets_final %>% 
	right_join(bind_rows(supplement_variant_sets_final_finemap_ukb, supplement_variant_sets_final_finemap_bbj)) %>% arrange(desc(FineMapped_cohort), FineMapped_trait)
write_tsv(supplement_variant_sets_final_finemap_ukb, "../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_SpliceAI_scores_overlap_fine_map.txt.gz")
