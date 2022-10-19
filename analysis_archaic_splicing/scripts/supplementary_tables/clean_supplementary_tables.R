#!/bin/R

library(tidyverse)
library(data.table)

# variant sets
supplement_variant_sets <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub.txt.gz"))
supplement_variant_sets_final <- supplement_variant_sets %>% 
	dplyr::select(hub_reference_genome,hub_variant_CHROM,hub_variant_POS,hub_variant_ID,hub_variant_REF,hub_variant_ALT,hub_variant_ANC,hub_variant_DER,hub_variant_CASE,hub_in_prufer_2014_archaic,hub_in_prufer_2014_modern,hub_in_castellano_2014_archaic,hub_in_castellano_2014_modern,hub_in_castellano_2014_nean,hub_in_vernot_akey_2016,hub_in_vernot_akey_2016_EAS,hub_in_vernot_akey_2016_EUR,hub_in_vernot_akey_2016_MEL,hub_in_vernot_akey_2016_SAS,hub_in_ldvernot_akey_2016,hub_in_ldvernot_akey_2016_EAS,hub_in_ldvernot_akey_2016_EUR,hub_in_ldvernot_akey_2016_MEL,hub_in_ldvernot_akey_2016_SAS,hub_in_gittelman_2016,hub_in_gittelman_2016_EAS,hub_in_gittelman_2016_EUR,hub_in_gittelman_2016_MEL,hub_in_gittelman_2016_SAS,hub_in_browning_2018,hub_in_browning_2018_BEB,hub_in_browning_2018_CDX,hub_in_browning_2018_CEU,hub_in_browning_2018_CHB,hub_in_browning_2018_CHS,hub_in_browning_2018_CLM,hub_in_browning_2018_FIN,hub_in_browning_2018_GBR,hub_in_browning_2018_GIH,hub_in_browning_2018_IBS,hub_in_browning_2018_ITU,hub_in_browning_2018_JPT,hub_in_browning_2018_KHV,hub_in_browning_2018_MXL,hub_in_browning_2018_PEL,hub_in_browning_2018_PJL,hub_in_browning_2018_PUR,hub_in_browning_2018_STU,hub_in_browning_2018_TSI,hub_in_browning_2018_Papuans,hub_in_racimo_2017_AI,hub_in_final_study_modern,hub_in_final_study_archaic,hub_in_final_study_nean,hub_in_final_study_deni,hub_in_final_study_introgressed,hub_in_final_study_adaptive,CHROM,POS,ID,REF,ALT,QUAL,FILTER,AC,AF,AN,NS,DP,EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF,EAS_AC,AMR_AC,AFR_AC,EUR_AC,SAS_AC,EUR_AC_bin,EUR_hapR2tag,EUR_hapR2tag_bin,gnomAD_ID,gnomAD_QUAL,gnomAD_FILTER,gnomAD_AC,gnomAD_AN,gnomAD_AF,gnomAD_AC_afr,gnomAD_AN_afr,gnomAD_AF_afr,archaic_mask_inter,altai_denisovan_AC,altai_denisovan_AN,altai_neanderthal_AC,altai_neanderthal_AN,chagyrskaya_neanderthal_AC,chagyrskaya_neanderthal_AN,vindija_neanderthal_AC,vindija_neanderthal_AN,archaic_AC,archaic_AN,neanderthal_AC,neanderthal_AN,denisovan_AC,denisovan_AN) %>% 
	dplyr::rename(main_reference_genome=hub_reference_genome,main_variant_CHROM=hub_variant_CHROM,main_variant_POS=hub_variant_POS,main_variant_ID=hub_variant_ID,main_variant_REF=hub_variant_REF,main_variant_ALT=hub_variant_ALT,main_variant_ANC=hub_variant_ANC,main_variant_DER=hub_variant_DER,main_variant_CASE=hub_variant_CASE,main_in_prufer_2014_archaic=hub_in_prufer_2014_archaic,main_in_prufer_2014_modern=hub_in_prufer_2014_modern,main_in_castellano_2014_archaic=hub_in_castellano_2014_archaic,main_in_castellano_2014_modern=hub_in_castellano_2014_modern,main_in_castellano_2014_nean=hub_in_castellano_2014_nean,main_in_vernot_akey_2016=hub_in_vernot_akey_2016,main_in_vernot_akey_2016_EAS=hub_in_vernot_akey_2016_EAS,main_in_vernot_akey_2016_EUR=hub_in_vernot_akey_2016_EUR,main_in_vernot_akey_2016_MEL=hub_in_vernot_akey_2016_MEL,main_in_vernot_akey_2016_SAS=hub_in_vernot_akey_2016_SAS,main_in_ldvernot_akey_2016=hub_in_ldvernot_akey_2016,main_in_ldvernot_akey_2016_EAS=hub_in_ldvernot_akey_2016_EAS,main_in_ldvernot_akey_2016_EUR=hub_in_ldvernot_akey_2016_EUR,main_in_ldvernot_akey_2016_MEL=hub_in_ldvernot_akey_2016_MEL,main_in_ldvernot_akey_2016_SAS=hub_in_ldvernot_akey_2016_SAS,main_in_gittelman_2016=hub_in_gittelman_2016,main_in_gittelman_2016_EAS=hub_in_gittelman_2016_EAS,main_in_gittelman_2016_EUR=hub_in_gittelman_2016_EUR,main_in_gittelman_2016_MEL=hub_in_gittelman_2016_MEL,main_in_gittelman_2016_SAS=hub_in_gittelman_2016_SAS,main_in_browning_2018=hub_in_browning_2018,main_in_browning_2018_BEB=hub_in_browning_2018_BEB,main_in_browning_2018_CDX=hub_in_browning_2018_CDX,main_in_browning_2018_CEU=hub_in_browning_2018_CEU,main_in_browning_2018_CHB=hub_in_browning_2018_CHB,main_in_browning_2018_CHS=hub_in_browning_2018_CHS,main_in_browning_2018_CLM=hub_in_browning_2018_CLM,main_in_browning_2018_FIN=hub_in_browning_2018_FIN,main_in_browning_2018_GBR=hub_in_browning_2018_GBR,main_in_browning_2018_GIH=hub_in_browning_2018_GIH,main_in_browning_2018_IBS=hub_in_browning_2018_IBS,main_in_browning_2018_ITU=hub_in_browning_2018_ITU,main_in_browning_2018_JPT=hub_in_browning_2018_JPT,main_in_browning_2018_KHV=hub_in_browning_2018_KHV,main_in_browning_2018_MXL=hub_in_browning_2018_MXL,main_in_browning_2018_PEL=hub_in_browning_2018_PEL,main_in_browning_2018_PJL=hub_in_browning_2018_PJL,main_in_browning_2018_PUR=hub_in_browning_2018_PUR,main_in_browning_2018_STU=hub_in_browning_2018_STU,main_in_browning_2018_TSI=hub_in_browning_2018_TSI,main_in_browning_2018_Papuans=hub_in_browning_2018_Papuans,main_in_racimo_2017_overlap=hub_in_racimo_2017_AI,main_in_final_study_modern=hub_in_final_study_modern,main_in_final_study_archaic=hub_in_final_study_archaic,main_in_final_study_nean=hub_in_final_study_nean,main_in_final_study_deni=hub_in_final_study_deni,main_in_final_study_introgressed=hub_in_final_study_introgressed,main_in_final_study_adaptive=hub_in_final_study_adaptive,'1KGP_CHROM'=CHROM,'1KGP_POS'=POS,'1KGP_ID'=ID,'1KGP_REF'=REF,'1KGP_ALT'=ALT,'1KGP_QUAL'=QUAL,'1KGP_FILTER'=FILTER,'1KGP_AC'=AC,'1KGP_AF'=AF,'1KGP_AN'=AN,'1KGP_NS'=NS,'1KGP_DP'=DP,'1KGP_EAS_AF'=EAS_AF,'1KGP_AMR_AF'=AMR_AF,'1KGP_AFR_AF'=AFR_AF,'1KGP_EUR_AF'=EUR_AF,'1KGP_SAS_AF'=SAS_AF,'1KGP_EAS_AC'=EAS_AC,'1KGP_AMR_AC'=AMR_AC,'1KGP_AFR_AC'=AFR_AC,'1KGP_EUR_AC'=EUR_AC,'1KGP_SAS_AC'=SAS_AC,'1KGP_EUR_AC_bin'=EUR_AC_bin,'1KGP_EUR_LD'=EUR_hapR2tag,'1KGP_EUR_LD_bin'=EUR_hapR2tag_bin,archaic_mask_intersection=archaic_mask_inter)
supplement_variant_spliceai <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
supplement_variant_spliceai_final <- supplement_variant_spliceai %>% 
	dplyr::select(hub_reference_genome,hub_variant_ID,ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL,spliceai_max) %>% 
	dplyr::rename(main_reference_genome=hub_reference_genome,main_variant_ID=hub_variant_ID,SpliceAI_ALLELE=ALLELE,SpliceAI_SYMBOL=SYMBOL,SpliceAI_DS_AG=DS_AG,SpliceAI_DS_AL=DS_AL,SpliceAI_DS_DG=DS_DG,SpliceAI_DS_DL=DS_DL,SpliceAI_DP_AG=DP_AG,SpliceAI_DP_AL=DP_AL,SpliceAI_DP_DG=DP_DG,SpliceAI_DP_DL=DP_DL,SpliceAI_max=spliceai_max)
supplement_variant_sets_final <- supplement_variant_sets_final %>% 
	left_join(supplement_variant_spliceai_final)
write_tsv(supplement_variant_sets_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_Human_variant_sets_and_SpliceAI.txt.gz"))

# 	tack on vep
supplement_variant_sets_vep <- as_tibble(fread("../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep.txt.gz"))
supplement_variant_sets_vep <- supplement_variant_sets_vep %>% 
	dplyr::select(hub_variant_ID,stop_gained,stop_lost,missense_variant,synonymous_variant,splice_donor_variant,splice_acceptor_variant,splice_region_variant,`5_prime_UTR_variant`,`3_prime_UTR_variant`,intron_variant,upstream_gene_variant,downstream_gene_variant,intergenic_variant) %>% 
	dplyr::rename(main_variant_ID=hub_variant_ID,VEP_stop_gained=stop_gained,VEP_stop_lost=stop_lost,VEP_missense_variant=missense_variant,VEP_synonymous_variant=synonymous_variant,VEP_splice_donor_variant=splice_donor_variant,VEP_splice_acceptor_variant=splice_acceptor_variant,VEP_splice_region_variant=splice_region_variant,`VEP_5_prime_UTR_variant`=`5_prime_UTR_variant`,`VEP_3_prime_UTR_variant`=`3_prime_UTR_variant`,VEP_intron_variant=intron_variant,VEP_upstream_gene_variant=upstream_gene_variant,VEP_downstream_gene_variant=downstream_gene_variant,VEP_intergenic_variant=intergenic_variant)
supplement_variant_sets_final <- supplement_variant_sets_final %>% 
	left_join(supplement_variant_sets_vep)
write_tsv(supplement_variant_sets_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_Human_variant_sets_and_SpliceAI.txt.gz"))

# mapsy results
supplement_mapsy_results <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
supplement_mapsy_results_final <- supplement_mapsy_results %>% 
	dplyr::select(source,construct_type,exon_seqnames,exon_start,exon_end,exon_width,exon_strand,exon_exon_id,exon_exon_name,exon_exon_rank,exon_exon_internal,exon_gene_id,exon_gene_type,exon_gene_name,exon_transcript_id,exon_transcript_type,exon_transcript_name,variant_seqnames,variant_start,variant_end,variant_width,variant_strand,variant_REF,variant_ALT,variant_QUAL,variant_FILTER,variant_name,variant_id,variant_type,ref_input_readc_a_wt,ref_input_readc_b_wt,ref_input_readc_c_wt,ref_input_readc_a_mt,ref_input_readc_b_mt,ref_input_readc_c_mt,ref_output_readc_a_wt,ref_output_readc_b_wt,ref_output_readc_c_wt,ref_output_readc_a_mt,ref_output_readc_b_mt,ref_output_readc_c_mt,filter_input_readcs,filter_output_readcs,mpralm.logFC,mpralm.CI.L,mpralm.CI.R,mpralm.t,mpralm.P.Value,mpralm.adj.Pval,mpralm.sigvar,mpralm.sigclass,mpralm.ANCDER.logFC,hub_variant_ID,hub_reference_genome,hub_variant_CHROM,hub_variant_POS,hub_variant_REF,hub_variant_ALT,hub_variant_ANC,hub_variant_DER,hub_variant_CASE) %>% 
	dplyr::rename(main_variant_ID=hub_variant_ID,main_reference_genome=hub_reference_genome,main_variant_CHROM=hub_variant_CHROM,main_variant_POS=hub_variant_POS,main_variant_REF=hub_variant_REF,main_variant_ALT=hub_variant_ALT,main_variant_ANC=hub_variant_ANC,main_variant_DER=hub_variant_DER,main_variant_CASE=hub_variant_CASE)
supplement_mapsy_results_final <- supplement_mapsy_results_final %>% 
	left_join(supplement_variant_sets_final)

# validation mapsy
supplement_mapsy_validation <- as_tibble(fread("../../results/validate_half_exons/validate_half_exons_order/Neanderthal_order_mapsy_table_short.txt.gz"))
supplement_mapsy_validation_final <- supplement_mapsy_validation
write_tsv(supplement_mapsy_validation, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_MaPSY_scores_half_exon_validation_experiment.txt.gz"))

# mapsy annotations
# 	hexamer scores
supplement_annotations_hexamer <- as_tibble(fread("../../results/distribution_splice_prediction/mapsy_variant_table_hexamer.txt.gz"))
supplement_annotations_hexamer_final <- supplement_annotations_hexamer %>% 
	dplyr::select(variant_id,exon_gene_name,ss_strand,ss_wt_elevenmer,ss_mt_elevenmer,ss_EI_mt,ss_EI_wt,ss_EI,ss_exonicA3SS_mt,ss_exonicA3SS_wt,ss_exonicA3SS,ss_exonicA5SS_mt,ss_exonicA5SS_wt,ss_exonicA5SS,ss_intronicA3SS_mt,ss_intronicA3SS_wt,ss_intronicA3SS,ss_intronicA5SS_mt,ss_intronicA5SS_wt,ss_intronicA5SS) %>% 
	dplyr::rename(Hexamer_strand=ss_strand,Hexamer_wt_elevenmer=ss_wt_elevenmer,Hexamer_mt_elevenmer=ss_mt_elevenmer,Hexamer_mt_EI=ss_EI_mt,Hexamer_wt_EI=ss_EI_wt,Hexamer_delta_EI=ss_EI,Hexamer_mt_exonicA3SS=ss_exonicA3SS_mt,Hexamer_wt_exonicA3SS=ss_exonicA3SS_wt,Hexamer_delta_exonicA3SS=ss_exonicA3SS,Hexamer_mt_exonicA5SS=ss_exonicA5SS_mt,Hexamer_wt_exonicA5SS=ss_exonicA5SS_wt,Hexamer_delta_exonicA5SS=ss_exonicA5SS,Hexamer_mt_intronicA3SS=ss_intronicA3SS_mt,Hexamer_wt_intronicA3SS=ss_intronicA3SS_wt,Hexamer_delta_intronicA3SS=ss_intronicA3SS,Hexamer_mt_intronicA5SS=ss_intronicA5SS_mt,Hexamer_wt_intronicA5SS=ss_intronicA5SS_wt,Hexamer_delta_intronicA5SS=ss_intronicA5SS)
supplement_annotations_hexamer_final <- supplement_mapsy_results_final %>% 
	left_join(supplement_annotations_hexamer_final)

# 	vep consequences
supplement_annotations_vep <- as_tibble(fread("../../results/distribution_splice_prediction/mapsy_variant_table_vep_long.txt.gz"))
supplement_annotations_vep_final <- supplement_annotations_vep %>% 
	dplyr::select(variant_id,exon_gene_name,VEP) %>% 
	dplyr::rename(VEP_consequence=VEP)
supplement_annotations_vep_final <- supplement_annotations_vep_final %>% 
	left_join(supplement_annotations_vep_final)

# 	spliceai gs raw
supplement_annotations_spliceai <- as_tibble(fread("../../results/distribution_splice_prediction/mapsy_variant_table_spliceai_gs_raw.txt.gz"))
supplement_annotations_spliceai_final <- supplement_annotations_spliceai %>% 
	dplyr::select(variant_id,exon_gene_name,ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL,spliceai_max) %>% 
	dplyr::rename(SpliceAI_ALLELE=ALLELE,SpliceAI_SYMBOL=SYMBOL,SpliceAI_DS_AG=DS_AG,SpliceAI_DS_AL=DS_AL,SpliceAI_DS_DG=DS_DG,SpliceAI_DS_DL=DS_DL,SpliceAI_DP_AG=DP_AG,SpliceAI_DP_AL=DP_AL,SpliceAI_DP_DG=DP_DG,SpliceAI_DP_DL=DP_DL,SpliceAI_max=spliceai_max)
supplement_annotations_spliceai_final <- supplement_mapsy_results_final %>% 
	left_join(supplement_annotations_spliceai_final)

# 	merge and save
supplement_mapsy_results_merged <- supplement_mapsy_results_final %>% 
	full_join(supplement_annotations_spliceai_final) %>% 
	full_join(supplement_annotations_hexamer_final) # %>% 
	# full_join(supplement_annotations_vep_final)

# 

# 
# 
# 

write_tsv(supplement_mapsy_results_merged, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_MaPSy_scores_main_experiment_SpliceAI_DeltaEI_VEP.txt.gz"))

# gtex sqtls and eqtls
supplement_sqtls_collated <- as_tibble(fread("../../results/preprocess_GTEx_QTLs/collate_GTEx_sQTLs_hg19.txt.gz"))
names(supplement_sqtls_collated) <- gsub("hub", "main", names(supplement_sqtls_collated))
supplement_eqtls_collated <- as_tibble(fread("../../results/preprocess_GTEx_QTLs/collate_GTEx_eQTLs_hg19.txt.gz"))
names(supplement_eqtls_collated) <- gsub("hub", "main", names(supplement_eqtls_collated))
write_tsv(supplement_sqtls_collated, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_GTEx_sQTL_aggregated_variants.txt.gz"))
write_tsv(supplement_eqtls_collated, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_GTEx_eQTL_aggregated_variants.txt.gz"))

# # 1kgp splice ai scores
# supplement_spliceai_weak <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_weak.txt.gz"))
# supplement_spliceai_moderate <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_moderate.txt.gz"))
# supplement_spliceai_strong <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_SpliceAI_bins_strong.txt.gz"))

# supplement_spliceai_weak_final <- supplement_spliceai_weak %>% 
# 	dplyr::select(hub_reference_genome,hub_variant_CHROM,hub_variant_POS,hub_variant_ID,hub_variant_REF,hub_variant_ALT,hub_variant_ANC,hub_variant_DER,hub_variant_CASE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,AC,AF,AN,NS,DP,EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF,EAS_AC,AMR_AC,AFR_AC,EUR_AC,SAS_AC,EUR_AC_bin,EUR_hapR2tag,EUR_hapR2tag_bin,ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL,spliceai_max) %>% 
# 	dplyr::rename(main_reference_genome=hub_reference_genome,main_variant_CHROM=hub_variant_CHROM,main_variant_POS=hub_variant_POS,main_variant_ID=hub_variant_ID,main_variant_REF=hub_variant_REF,main_variant_ALT=hub_variant_ALT,main_variant_ANC=hub_variant_ANC,main_variant_DER=hub_variant_DER,main_variant_CASE=hub_variant_CASE,`1KGP_CHROM`=CHROM,`1KGP_POS`=POS,`1KGP_ID`=ID,`1KGP_REF`=REF,`1KGP_ALT`=ALT,`1KGP_QUAL`=QUAL,`1KGP_FILTER`=FILTER,`1KGP_AC`=AC,`1KGP_AF`=AF,`1KGP_AN`=AN,`1KGP_NS`=NS,`1KGP_DP`=DP,`1KGP_EAS_AF`=EAS_AF,`1KGP_AMR_AF`=AMR_AF,`1KGP_AFR_AF`=AFR_AF,`1KGP_EUR_AF`=EUR_AF,`1KGP_SAS_AF`=SAS_AF,`1KGP_EAS_AC`=EAS_AC,`1KGP_AMR_AC`=AMR_AC,`1KGP_AFR_AC`=AFR_AC,`1KGP_EUR_AC`=EUR_AC,`1KGP_SAS_AC`=SAS_AC,`1KGP_EUR_AC_bin`=EUR_AC_bin,`1KGP_EUR_LD`=EUR_hapR2tag,`1KGP_EUR_LD_bin`=EUR_hapR2tag_bin,SpliceAI_ALLELE=ALLELE,SpliceAI_SYMBOL=SYMBOL,SpliceAI_DS_AG=DS_AG,SpliceAI_DS_AL=DS_AL,SpliceAI_DS_DG=DS_DG,SpliceAI_DS_DL=DS_DL,SpliceAI_DP_AG=DP_AG,SpliceAI_DP_AL=DP_AL,SpliceAI_DP_DG=DP_DG,SpliceAI_DP_DL=DP_DL,SpliceAI_max=spliceai_max) %>% 
# 	mutate(SpliceAI_weak=TRUE)
# supplement_spliceai_moderate_final <- supplement_spliceai_moderate %>% 
# 	dplyr::select(hub_reference_genome,hub_variant_CHROM,hub_variant_POS,hub_variant_ID,hub_variant_REF,hub_variant_ALT,hub_variant_ANC,hub_variant_DER,hub_variant_CASE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,AC,AF,AN,NS,DP,EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF,EAS_AC,AMR_AC,AFR_AC,EUR_AC,SAS_AC,EUR_AC_bin,EUR_hapR2tag,EUR_hapR2tag_bin,ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL,spliceai_max) %>% 
# 	dplyr::rename(main_reference_genome=hub_reference_genome,main_variant_CHROM=hub_variant_CHROM,main_variant_POS=hub_variant_POS,main_variant_ID=hub_variant_ID,main_variant_REF=hub_variant_REF,main_variant_ALT=hub_variant_ALT,main_variant_ANC=hub_variant_ANC,main_variant_DER=hub_variant_DER,main_variant_CASE=hub_variant_CASE,`1KGP_CHROM`=CHROM,`1KGP_POS`=POS,`1KGP_ID`=ID,`1KGP_REF`=REF,`1KGP_ALT`=ALT,`1KGP_QUAL`=QUAL,`1KGP_FILTER`=FILTER,`1KGP_AC`=AC,`1KGP_AF`=AF,`1KGP_AN`=AN,`1KGP_NS`=NS,`1KGP_DP`=DP,`1KGP_EAS_AF`=EAS_AF,`1KGP_AMR_AF`=AMR_AF,`1KGP_AFR_AF`=AFR_AF,`1KGP_EUR_AF`=EUR_AF,`1KGP_SAS_AF`=SAS_AF,`1KGP_EAS_AC`=EAS_AC,`1KGP_AMR_AC`=AMR_AC,`1KGP_AFR_AC`=AFR_AC,`1KGP_EUR_AC`=EUR_AC,`1KGP_SAS_AC`=SAS_AC,`1KGP_EUR_AC_bin`=EUR_AC_bin,`1KGP_EUR_LD`=EUR_hapR2tag,`1KGP_EUR_LD_bin`=EUR_hapR2tag_bin,SpliceAI_ALLELE=ALLELE,SpliceAI_SYMBOL=SYMBOL,SpliceAI_DS_AG=DS_AG,SpliceAI_DS_AL=DS_AL,SpliceAI_DS_DG=DS_DG,SpliceAI_DS_DL=DS_DL,SpliceAI_DP_AG=DP_AG,SpliceAI_DP_AL=DP_AL,SpliceAI_DP_DG=DP_DG,SpliceAI_DP_DL=DP_DL,SpliceAI_max=spliceai_max) %>% 
# 	mutate(SpliceAI_moderate=TRUE)
# supplement_spliceai_strong_final <- supplement_spliceai_strong %>% 
# 	dplyr::select(hub_reference_genome,hub_variant_CHROM,hub_variant_POS,hub_variant_ID,hub_variant_REF,hub_variant_ALT,hub_variant_ANC,hub_variant_DER,hub_variant_CASE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,AC,AF,AN,NS,DP,EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF,EAS_AC,AMR_AC,AFR_AC,EUR_AC,SAS_AC,EUR_AC_bin,EUR_hapR2tag,EUR_hapR2tag_bin,ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL,spliceai_max) %>% 
# 	dplyr::rename(main_reference_genome=hub_reference_genome,main_variant_CHROM=hub_variant_CHROM,main_variant_POS=hub_variant_POS,main_variant_ID=hub_variant_ID,main_variant_REF=hub_variant_REF,main_variant_ALT=hub_variant_ALT,main_variant_ANC=hub_variant_ANC,main_variant_DER=hub_variant_DER,main_variant_CASE=hub_variant_CASE,`1KGP_CHROM`=CHROM,`1KGP_POS`=POS,`1KGP_ID`=ID,`1KGP_REF`=REF,`1KGP_ALT`=ALT,`1KGP_QUAL`=QUAL,`1KGP_FILTER`=FILTER,`1KGP_AC`=AC,`1KGP_AF`=AF,`1KGP_AN`=AN,`1KGP_NS`=NS,`1KGP_DP`=DP,`1KGP_EAS_AF`=EAS_AF,`1KGP_AMR_AF`=AMR_AF,`1KGP_AFR_AF`=AFR_AF,`1KGP_EUR_AF`=EUR_AF,`1KGP_SAS_AF`=SAS_AF,`1KGP_EAS_AC`=EAS_AC,`1KGP_AMR_AC`=AMR_AC,`1KGP_AFR_AC`=AFR_AC,`1KGP_EUR_AC`=EUR_AC,`1KGP_SAS_AC`=SAS_AC,`1KGP_EUR_AC_bin`=EUR_AC_bin,`1KGP_EUR_LD`=EUR_hapR2tag,`1KGP_EUR_LD_bin`=EUR_hapR2tag_bin,SpliceAI_ALLELE=ALLELE,SpliceAI_SYMBOL=SYMBOL,SpliceAI_DS_AG=DS_AG,SpliceAI_DS_AL=DS_AL,SpliceAI_DS_DG=DS_DG,SpliceAI_DS_DL=DS_DL,SpliceAI_DP_AG=DP_AG,SpliceAI_DP_AL=DP_AL,SpliceAI_DP_DG=DP_DG,SpliceAI_DP_DL=DP_DL,SpliceAI_max=spliceai_max) %>% 
# 	mutate(SpliceAI_strong=TRUE)

# supplement_spliceai_final <- supplement_spliceai_weak_final %>% 
# 	full_join(supplement_spliceai_moderate_final) %>% 
# 	full_join(supplement_spliceai_strong_final) %>% 
# 	mutate(SpliceAI_weak = ifelse(is.na(SpliceAI_weak), FALSE, SpliceAI_weak)) %>% 
# 	mutate(SpliceAI_moderate = ifelse(is.na(SpliceAI_moderate), FALSE, SpliceAI_moderate)) %>% 
# 	mutate(SpliceAI_strong = ifelse(is.na(SpliceAI_strong), FALSE, SpliceAI_strong))
# write_tsv(supplement_spliceai_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_1KGP_EUR_MAF0.01_SpliceAI_categories.txt.gz"))

# # 1kgp eur maf 0.01
# supplement_1KGP <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub.txt.gz"))
# supplement_1KGP_final <- supplement_1KGP %>% 
# 	dplyr::select(hub_reference_genome,hub_variant_CHROM,hub_variant_POS,hub_variant_ID,hub_variant_REF,hub_variant_ALT,hub_variant_ANC,hub_variant_DER,hub_variant_CASE,CHROM,POS,ID,REF,ALT,QUAL,FILTER,AC,AF,AN,NS,DP,EAS_AF,AMR_AF,AFR_AF,EUR_AF,SAS_AF,EAS_AC,AMR_AC,AFR_AC,EUR_AC,SAS_AC,EUR_AC_bin,EUR_hapR2tag,EUR_hapR2tag_bin) %>% 
# 	dplyr::rename(main_reference_genome=hub_reference_genome,main_variant_CHROM=hub_variant_CHROM,main_variant_POS=hub_variant_POS,main_variant_ID=hub_variant_ID,main_variant_REF=hub_variant_REF,main_variant_ALT=hub_variant_ALT,main_variant_ANC=hub_variant_ANC,main_variant_DER=hub_variant_DER,main_variant_CASE=hub_variant_CASE,`1KGP_CHROM`=CHROM,`1KGP_POS`=POS,`1KGP_ID`=ID,`1KGP_REF`=REF,`1KGP_ALT`=ALT,`1KGP_QUAL`=QUAL,`1KGP_FILTER`=FILTER,`1KGP_AC`=AC,`1KGP_AF`=AF,`1KGP_AN`=AN,`1KGP_NS`=NS,`1KGP_DP`=DP,`1KGP_EAS_AF`=EAS_AF,`1KGP_AMR_AF`=AMR_AF,`1KGP_AFR_AF`=AFR_AF,`1KGP_EUR_AF`=EUR_AF,`1KGP_SAS_AF`=SAS_AF,`1KGP_EAS_AC`=EAS_AC,`1KGP_AMR_AC`=AMR_AC,`1KGP_AFR_AC`=AFR_AC,`1KGP_EUR_AC`=EUR_AC,`1KGP_SAS_AC`=SAS_AC,`1KGP_EUR_AC_bin`=EUR_AC_bin,`1KGP_EUR_LD`=EUR_hapR2tag,`1KGP_EUR_LD_bin`=EUR_hapR2tag_bin)
# write_tsv(supplement_1KGP_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_1KGP_EUR_MAF0.01_annotated.txt.gz"))

# # 	add on column of spliceai scores
# supplement_1KGP_spliceai <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_hub_spliceai_gs_raw_dedup.gz"))
# supplement_1KGP_spliceai <- supplement_1KGP_spliceai %>% 
# 	dplyr::select(hub_variant_ID,ALLELE,SYMBOL,DS_AG,DS_AL,DS_DG,DS_DL,DP_AG,DP_AL,DP_DG,DP_DL,spliceai_max) %>% 
# 	dplyr::rename(main_variant_ID=hub_variant_ID,SpliceAI_ALLELE=ALLELE,SpliceAI_SYMBOL=SYMBOL,SpliceAI_DS_AG=DS_AG,SpliceAI_DS_AL=DS_AL,SpliceAI_DS_DG=DS_DG,SpliceAI_DS_DL=DS_DL,SpliceAI_DP_AG=DP_AG,SpliceAI_DP_AL=DP_AL,SpliceAI_DP_DG=DP_DG,SpliceAI_DP_DL=DP_DL,SpliceAI_max=spliceai_max)
# supplement_1KGP_final <- supplement_1KGP_final %>% 
# 	left_join(supplement_1KGP_spliceai)
# write_tsv(supplement_1KGP_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_1KGP_EUR_MAF0.01_annotated.txt.gz"))

# # 	add column for not-introgressed
# supplement_1KGP_notAI <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub.txt.gz"))
# supplement_1KGP_notAI <- supplement_1KGP_notAI %>% 
# 	dplyr::select(hub_variant_ID) %>% 
# 	dplyr::rename(main_variant_ID=hub_variant_ID)
# supplement_1KGP_final <- supplement_1KGP_final %>% 
# 	mutate(Non_introgressed = ifelse(main_variant_ID %in% supplement_1KGP_notAI$main_variant_ID, TRUE, FALSE))
# write_tsv(supplement_1KGP_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_1KGP_EUR_MAF0.01_annotated.txt.gz"))
# # supplement_1KGP_final <- as_tibble(fread("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_1KGP_EUR_MAF0.01_annotated.txt.gz"))

# # 	tack on vep
# supplement_1KGP_vep <- as_tibble(fread("../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.txt.gz"))
# supplement_1KGP_vep <- supplement_1KGP_vep %>% 
# 	dplyr::select(hub_variant_ID,stop_gained,stop_lost,missense_variant,synonymous_variant,splice_donor_variant,splice_acceptor_variant,splice_region_variant,`5_prime_UTR_variant`,`3_prime_UTR_variant`,intron_variant,upstream_gene_variant,downstream_gene_variant,intergenic_variant) %>% 
# 	dplyr::rename(main_variant_ID=hub_variant_ID,VEP_stop_gained=stop_gained,VEP_stop_lost=stop_lost,VEP_missense_variant=missense_variant,VEP_synonymous_variant=synonymous_variant,VEP_splice_donor_variant=splice_donor_variant,VEP_splice_acceptor_variant=splice_acceptor_variant,VEP_splice_region_variant=splice_region_variant,`VEP_5_prime_UTR_variant`=`5_prime_UTR_variant`,`VEP_3_prime_UTR_variant`=`3_prime_UTR_variant`,VEP_intron_variant=intron_variant,VEP_upstream_gene_variant=upstream_gene_variant,VEP_downstream_gene_variant=downstream_gene_variant,VEP_intergenic_variant=intergenic_variant)
# supplement_1KGP_final <- supplement_1KGP_final %>% 
# 	left_join(supplement_1KGP_vep)
# write_tsv(supplement_1KGP_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_1KGP_EUR_MAF0.01_annotated.txt.gz"))

# # 	add on columns of spliceai categories
# supplement_1KGP_final <- supplement_1KGP_final %>% 
# 	left_join(supplement_spliceai_final %>% 
# 		dplyr::select(main_variant_ID,SpliceAI_weak,SpliceAI_moderate,SpliceAI_strong)) %>% 
# 	mutate(SpliceAI_weak = ifelse(is.na(SpliceAI_weak), FALSE, SpliceAI_weak)) %>% 
# 	mutate(SpliceAI_moderate = ifelse(is.na(SpliceAI_moderate), FALSE, SpliceAI_moderate)) %>% 
# 	mutate(SpliceAI_strong = ifelse(is.na(SpliceAI_strong), FALSE, SpliceAI_strong))
# write_tsv(supplement_1KGP_final, gzfile("../../results/supplementary_tables/Rong_et_al_Supplementary_Data_SX_1KGP_EUR_MAF0.01_annotated.txt.gz"))
