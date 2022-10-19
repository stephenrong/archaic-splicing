#!/bin/R

library(tidyverse)
library(data.table)
library(wrapr)

# list of tissues
list_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

# process each sQTLs
print("sQTLs")
collate_GTEx_sQTLs_hg19 <- NULL
for (tissue in list_tissues) {
	print(tissue)
	collate_GTEx_sQTLs_hg19_temp <- as_tibble(fread(paste("../../results/preprocess_GTEx_QTLs/sQTLs_tb_final_hg19-", tissue, ".txt.gz", sep=""))) %>% dplyr::select(starts_with("hub")) %>% 
		group_by(hub_reference_genome, hub_variant_CHROM, hub_variant_POS, hub_variant_ID, hub_variant_REF, hub_variant_ALT, hub_variant_ANC, hub_variant_DER, hub_variant_CASE) %>% 
		summarise(QTL_count = n()) %>% ungroup()
	# add to list
	collate_GTEx_sQTLs_hg19[[paste(tissue, "sQTL", sep="_")]] <- collate_GTEx_sQTLs_hg19_temp
}

collate_GTEx_sQTLs_hg19_final <- collate_GTEx_sQTLs_hg19 %>% bind_rows(.id="QTL_tissue") %>% pivot_wider(names_from=QTL_tissue, values_from=QTL_count)

# save collated sQTLs
write_tsv(collate_GTEx_sQTLs_hg19_final, gzfile("../../results/preprocess_GTEx_QTLs/collate_GTEx_sQTLs_hg19.txt.gz"))

# process each eQTLs
print("eQTLs")
collate_GTEx_eQTLs_hg19 <- NULL
for (tissue in list_tissues) {
	print(tissue)
	collate_GTEx_eQTLs_hg19_temp <- as_tibble(fread(paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg19-", tissue, ".txt.gz", sep=""))) %>% dplyr::select(starts_with("hub")) %>% 
		group_by(hub_reference_genome, hub_variant_CHROM, hub_variant_POS, hub_variant_ID, hub_variant_REF, hub_variant_ALT, hub_variant_ANC, hub_variant_DER, hub_variant_CASE) %>% 
		summarise(QTL_count = n()) %>% ungroup()
	# add to list
	collate_GTEx_eQTLs_hg19[[paste(tissue, "eQTL", sep="_")]] <- collate_GTEx_eQTLs_hg19_temp
}

collate_GTEx_eQTLs_hg19_final <- collate_GTEx_eQTLs_hg19 %>% bind_rows(.id="QTL_tissue") %>% pivot_wider(names_from=QTL_tissue, values_from=QTL_count)

# save collated eQTLs
write_tsv(collate_GTEx_eQTLs_hg19_final, gzfile("../../results/preprocess_GTEx_QTLs/collate_GTEx_eQTLs_hg19.txt.gz"))
