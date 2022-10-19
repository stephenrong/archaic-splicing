#!/bin/R
library(tidyverse)
library(data.table)
library(plyranges)
source("../get_helper.R")
library(EnsDb.Hsapiens.v86)
# hg38tohg19 <- import.chain("../../data/annotate_lift_over/hg38ToHg19.over.chain")

# list of tissues
list_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

# GTEx v8 eQTLs
# eQTLs_tb_final <- NULL
for (tissue in list_tissues) {
	print(tissue)
	file_name <- paste("../../data/annotate_GTEx_eQTLs/GTEx_Analysis_v8_eQTL/", tissue, ".v8.signif_variant_gene_pairs.txt.gz", sep="")
	eQTLs_tb <- as_tibble(fread(file_name))
	names(eQTLs_tb) <- paste("QTLs_", names(eQTLs_tb), sep="")

	# get hub columns, qtl columns, and remap gene ids
	eQTLs_tb_tissue <- eQTLs_tb %>% 
		separate(QTLs_variant_id, into=c("hub_variant_CHROM", "hub_variant_POS", "hub_variant_REF", "hub_variant_ALT", "hub_reference_genome"), sep="_", remove=F) %>% 
		mutate(hub_reference_genome = "hg38") %>% 
		mutate(hub_variant_CHROM = gsub("chr", "", hub_variant_CHROM)) %>% 
		mutate(hub_variant_POS = as.numeric(hub_variant_POS)) %>% 
		mutate_hub_variant_ID() %>% 
		mutate(QTLs_CHROM = hub_variant_CHROM) %>% 
		mutate(QTLs_POS = hub_variant_POS) %>% 
		mutate(QTLs_REF = hub_variant_REF) %>% 
		mutate(QTLs_ALT = hub_variant_ALT) %>% 
		mutate(QTLs_type = "eQTL") %>% 
		mutate(QTLs_tissue = tissue) %>% 
		mutate(QTLs_gene_id = gsub("\\..*", "", QTLs_gene_id))

	# map gene ids to gene names
	ensembl_map_eQTL <- ensembldb::select(EnsDb.Hsapiens.v86,
		keys=unique(eQTLs_tb_tissue$QTLs_gene_id), 
		keytype="GENEID", columns=c("SYMBOL", "GENEID")
	) %>% as_tibble() %>% 
		dplyr::rename(QTLs_gene_id = GENEID, QTLs_gene_name = SYMBOL)
	table(table(ensembl_map_eQTL$QTLs_gene_id))
	table(table(ensembl_map_eQTL$QTLs_gene_name))
	# don't bother with repeats
	ensembl_map_eQTL <- ensembl_map_eQTL %>% 
		dplyr::filter(!duplicated(QTLs_gene_name), !duplicated(QTLs_gene_id))

	# get col
	eQTLs_tb_tissue <- eQTLs_tb_tissue %>% 
		left_join(ensembl_map_eQTL) %>% 
		mutate(QTLs_gene_name = ifelse(is.na(QTLs_gene_name)|
			(QTLs_gene_name=="NA"), "", QTLs_gene_name))

	# get ancestral
	eQTLs_tb_tissue <- eQTLs_tb_tissue %>% 
		mutate_hub_variant_ANC() %>% 
		mutate_hub_variant_DER() %>% 
		mutate_hub_variant_ID() %>% 
		split_ancestral()

	# sort cols
	eQTLs_tb_tissue <- eQTLs_tb_tissue %>% 
		dplyr::select(hub_reference_genome, hub_variant_CHROM, hub_variant_POS, hub_variant_ID, 
			hub_variant_REF, hub_variant_ALT, hub_variant_ANC, hub_variant_DER, hub_variant_CASE, 
			starts_with("QTLs_")) %>% 
		standard_sort()

	# save tissue QTL file
	write_tsv(eQTLs_tb_tissue, gzfile(paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg38-", tissue, ".txt.gz", sep="")))
	
	# save tissue VCF file
	hub_to_vcf(eQTLs_tb_tissue, paste("eQTLs_tb_final_hg38-", tissue, sep=""), 
		paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg38-", tissue, ".vcf", sep=""), 
		BSgenome.Hsapiens.UCSC.hg38, partial=FALSE, compress=TRUE)

	# liftover VCF file
	hg38ToHg39 <- "../../data/annotate_lift_over/hg38ToHg19.over.chain"
	hg19fasta <- "../../data/fasta/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz"
	in_vcf_gz <- paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg38-", tissue, ".vcf.gz", sep="")
	out_vcf <- paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg19-", tissue, ".vcf", sep="")
	system(paste("CrossMap.py vcf", hg38ToHg39, in_vcf_gz, hg19fasta, out_vcf, sep=" "))
	system(paste("bgzip -f", out_vcf, sep=" "))

	# load and clean up
	out_vcf_gz <- paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg19-", tissue, ".vcf.gz", sep="")
	eQTLs_tb_tissue_lift <- read.vcfR(out_vcf_gz)@fix %>% as_tibble()
	INFO_cols <- strsplit(gsub("=.*$", "", gsub("=.*?;", ",", eQTLs_tb_tissue_lift$INFO[1])), ",")[[1]]
	eQTLs_tb_tissue_lift <- eQTLs_tb_tissue_lift %>% 
		separate(col=INFO, into=INFO_cols, sep=";")
	for (col in INFO_cols) {
		let(c(VAR=col), eQTLs_tb_tissue_lift <- eQTLs_tb_tissue_lift %>% 
			mutate(VAR = gsub(".*=", "", VAR)))
	}
	eQTLs_tb_tissue_final <- eQTLs_tb_tissue_lift %>% 
		dplyr::select(-c(QUAL, FILTER)) %>% 
		dplyr::rename(
			hub_variant_CHROM = CHROM,
			hub_variant_POS = POS, 
			hub_variant_ID = ID, 
			hub_variant_REF = REF, 
			hub_variant_ALT = ALT) %>% 
		mutate(hub_reference_genome = "hg19") %>% 
		mutate(hub_variant_POS = as.numeric(hub_variant_POS)) %>% 
		mutate_hub_variant_ID()  # overwrite

	# get ancestral
	eQTLs_tb_tissue_final <- eQTLs_tb_tissue_final %>% 
		mutate_hub_variant_ANC() %>% 
		mutate_hub_variant_DER() %>% 
		mutate_hub_variant_ID() %>% 
		split_ancestral()

	# sort cols
	eQTLs_tb_tissue_final <- eQTLs_tb_tissue_final %>% 
		dplyr::select(hub_reference_genome, hub_variant_CHROM, hub_variant_POS, hub_variant_ID, 
			hub_variant_REF, hub_variant_ALT, hub_variant_ANC, hub_variant_DER, hub_variant_CASE, 
			starts_with("QTLs_")) %>% 
		standard_sort()

	# save tissue QTL file
	write_tsv(eQTLs_tb_tissue_final, gzfile(paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg19-", tissue, ".txt.gz", sep="")))
	
	# save tissue VCF file
	hub_to_vcf(eQTLs_tb_tissue_final, paste("eQTLs_tb_final_hg19-", tissue, sep=""), 
		paste("../../results/preprocess_GTEx_QTLs/eQTLs_tb_final_hg19-", tissue, ".vcf", sep=""), 
		BSgenome.Hsapiens.UCSC.hg38, partial=FALSE, compress=TRUE)
}
