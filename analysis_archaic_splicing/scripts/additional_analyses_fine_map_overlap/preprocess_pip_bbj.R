#!/bin/R

library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)

trait_tb <- as_tibble(fread("../../data/finemap_overlap/bbj/manifest_Kanai_2021_fine_mapped_bbj.txt"))
prefix_file <- "../../data/finemap_overlap/bbj/pheweb.jp/"
pheno_list <- trait_tb$Phenotype

bbj_fine_mapping_FINEMAP_list <- list()
parse_into_bed <- function(pheno) {
	file_name <- paste(prefix_file, filter(trait_tb, Phenotype==pheno)$`FINEMAP file`, sep="")
	temp_list <- as_tibble(fread(file_name)) %>% 
		mutate(chromosome = as.character(chromosome)) %>% 
		filter(pip >= 0.05)
	return(temp_list)
}

for (pheno in pheno_list) {
	print(pheno)
	bbj_fine_mapping_FINEMAP_list[[pheno]] <- parse_into_bed(pheno)
}

bbj_fine_mapping_FINEMAP_final <- bind_rows(bbj_fine_mapping_FINEMAP_list, .id="trait")
print(bbj_fine_mapping_FINEMAP_list)
print(bbj_fine_mapping_FINEMAP_final)
write_tsv(bbj_fine_mapping_FINEMAP_final, gzfile("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_FINEMAP_final.txt.gz"))
bbj_fine_mapping_FINEMAP_final <- NULL

bbj_fine_mapping_SUSIE_list <- list()
parse_into_bed <- function(pheno) {
	file_name <- paste(prefix_file, filter(trait_tb, Phenotype==pheno)$`SUSIE file`, sep="")
	temp_list <- as_tibble(fread(file_name)) %>% 
		mutate(chromosome = as.character(chromosome)) %>% 
		filter(pip >= 0.05)
	return(temp_list)
}

for (pheno in pheno_list) {
	print(pheno)
	bbj_fine_mapping_SUSIE_list[[pheno]] <- parse_into_bed(pheno)
}

bbj_fine_mapping_SUSIE_final <- bind_rows(bbj_fine_mapping_SUSIE_list, .id="trait")
print(bbj_fine_mapping_SUSIE_list)
print(bbj_fine_mapping_SUSIE_final)
write_tsv(bbj_fine_mapping_SUSIE_final, gzfile("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_SUSIE_final.txt.gz"))
bbj_fine_mapping_SUSIE_final <- NULL
