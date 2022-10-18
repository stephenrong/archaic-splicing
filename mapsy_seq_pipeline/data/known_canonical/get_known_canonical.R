#!/bin/R

library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)
source("../../scripts/mapsy_helper_v1.R")

# get gencode basic
gencode_v32_basic <- get_gff("../gff/gencode.v32.basic.annotation.gff3.gz")
gencode_v32_lift37_basic <- get_gff("../gff/gencode.v32lift37.basic.annotation.gff3.gz")

# get canonical hg38
# 	get canonical transcript by appris
gencode_v32_basic_canonical_transcripts <- gencode_v32_basic %>% 
	as_tibble() %>% filter(type == "transcript") %>% 
	mutate(appris_ranking = case_when(
		grepl("appris_principal_1", tag) ~ 1,
		grepl("appris_principal_2", tag) ~ 2, 
		grepl("appris_principal_3", tag) ~ 3, 
		grepl("appris_principal_4", tag) ~ 4, 
		grepl("appris_principal_5", tag) ~ 5, 
		grepl("appris_alternative_1", tag) ~ 6, 
		grepl("appris_alternative_2", tag) ~ 7, 
		TRUE ~ Inf
	)) %>% 
	filter(is.finite(appris_ranking))
# 	for genes with multiple transcripts, choose longest transcript
gencode_v32_basic_canonical_transcripts <- gencode_v32_basic_canonical_transcripts %>% 
	group_by(gene_name) %>% 
	arrange(appris_ranking, desc(width)) %>% 
	mutate(canonical_transcript = (!duplicated(gene_name))) %>% filter(canonical_transcript) %>% 
	ungroup()
# 	get canonical exons
gencode_v32_basic_canonical_exons <- gencode_v32_basic %>% 
	as_tibble() %>% filter(type == "exon") %>% 
	filter(transcript_id %in% gencode_v32_basic_canonical_transcripts$transcript_id) %>% 
	dplyr::select(seqnames, start, end, width, strand, gene_name, transcript_id, transcript_name, exon_number, exon_id) %>% 
	dplyr::rename(canonical_transcript_id=transcript_id, canonical_transcript_name=transcript_name, canonical_exon_number=exon_number, canonical_exon_id=exon_id) %>% 
	mutate(canonical_exon = TRUE)

# save
write_tsv(gencode_v32_basic_canonical_transcripts, gzfile("gencode_v32_basic_canonical_transcripts.txt.gz"))
write_tsv(gencode_v32_basic_canonical_exons, gzfile("gencode_v32_basic_canonical_exons.txt.gz"))

# get canonical hg3y
# 	get canonical transcript by appris
gencode_v32_lift37_basic_canonical_transcripts <- gencode_v32_lift37_basic %>% 
	as_tibble() %>% filter(type == "transcript") %>% 
	mutate(appris_ranking = case_when(
		grepl("appris_principal_1", tag) ~ 1,
		grepl("appris_principal_2", tag) ~ 2, 
		grepl("appris_principal_3", tag) ~ 3, 
		grepl("appris_principal_4", tag) ~ 4, 
		grepl("appris_principal_5", tag) ~ 5, 
		grepl("appris_alternative_1", tag) ~ 6, 
		grepl("appris_alternative_2", tag) ~ 7, 
		TRUE ~ Inf
	)) %>% 
	filter(is.finite(appris_ranking))
# 	for genes with multiple transcripts, choose longest transcript
gencode_v32_lift37_basic_canonical_transcripts <- gencode_v32_lift37_basic_canonical_transcripts %>% 
	group_by(gene_name) %>% 
	arrange(appris_ranking, desc(width)) %>% 
	mutate(canonical_transcript = (!duplicated(gene_name))) %>% filter(canonical_transcript) %>% 
	ungroup()
# 	get canonical exons
gencode_v32_lift37_basic_canonical_exons <- gencode_v32_lift37_basic %>% 
	as_tibble() %>% filter(type == "exon") %>% 
	filter(transcript_id %in% gencode_v32_lift37_basic_canonical_transcripts$transcript_id) %>% 
	dplyr::select(seqnames, start, end, width, strand, gene_name, transcript_id, transcript_name, exon_number, exon_id) %>% 
	dplyr::rename(canonical_transcript_id=transcript_id, canonical_transcript_name=transcript_name, canonical_exon_number=exon_number, canonical_exon_id=exon_id) %>% 
	mutate(canonical_exon = TRUE)

# save
write_tsv(gencode_v32_lift37_basic_canonical_transcripts, gzfile("gencode_v32_lift37_basic_canonical_transcripts.txt.gz"))
write_tsv(gencode_v32_lift37_basic_canonical_exons, gzfile("gencode_v32_lift37_basic_canonical_exons.txt.gz"))
