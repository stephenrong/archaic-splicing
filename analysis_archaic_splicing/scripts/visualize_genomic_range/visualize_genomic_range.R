#!/bin/R

library(tidyverse)
library(data.table)
library(plyranges)
library(ggbio)
library(viridis)
library(wrapr)
library(Gviz)
library(erma)
library(BSgenome.Hsapiens.UCSC.hg19)

reverse_complement <- function(x) {
	return(reverse(chartr(old="ATCGatcg", new="TAGCtagc", x)))
}

visualize_genomic_range <- function(vis_variant_id, win_half_size_gtex=NA, win_half_size_splice=NA, win_half_size_seq=NA, out_folder) {  # , col_REF="#FBBD7A", col_ALT="#377FA8") {
	# parse variant ID
	vis_variant_id_split <- strsplit(vis_variant_id, ":|_|/")[[1]]
	vis_CHROM <- gsub("chr", "", vis_variant_id_split[1])
	vis_POS <- as.numeric(vis_variant_id_split[2])
	vis_REF <- vis_variant_id_split[3]
	vis_ALT <- vis_variant_id_split[4]

	# load mapsy
	mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
	final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
	vis_SYMBOL1 <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_gene_name  # get gene name
	vis_SYMBOL2 <- final_v2_spliceai %>% filter(hub_variant_ID == vis_variant_id) %>% .$SYMBOL  # get gene name
	vis_SYMBOL <- unique(c(vis_SYMBOL1, vis_SYMBOL2))
	# vis_SYMBOL <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_gene_name  # get gene name
	vis_STRAND <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_strand  # get gene name
	mapsy_variant_table <- mapsy_variant_table %>%  # restrict to target gene
		filter(hub_variant_CHROM == vis_CHROM)
	temp <- mapsy_variant_table %>% 
		filter(hub_variant_ID == vis_variant_id) %>%  # .[,169:174]
		dplyr::select(hub_in_final_study_introgressed, hub_in_final_study_adaptive, hub_in_final_study_modern, hub_in_final_study_archaic, hub_in_final_study_nean, hub_in_final_study_deni)
	temp_list <- c()
	for (i in names(temp)) {
		if (temp[[i]][[1]]) {
			temp_list <- c(temp_list, i)
		}
	}
	# if (("hub_in_final_study_adaptive" %in% temp_list) | ("hub_in_final_study_introgressed" %in% temp_list)) {
	mapsy_variant_table <- mapsy_variant_table %>% 
		filter(exon_gene_name %in% c(vis_SYMBOL))
	# }
	temp_filter <- paste(temp_list[[1]], collapse="|")
	print(temp_filter)
	mapsy_variant_table <- mapsy_variant_table %>% 
		filter(!!rlang::parse_expr(temp_filter))

	# create mapsy grange
	mapsy_variant_table_gr <- 
		GRanges(seqnames=mapsy_variant_table$hub_variant_CHROM, IRanges(start=mapsy_variant_table$hub_variant_POS, end=mapsy_variant_table$hub_variant_POS), # strand=mapsy_variant_table$exon_strand, 
			hub_variant_ID=mapsy_variant_table$hub_variant_ID, hub_variant_CHROM=mapsy_variant_table$hub_variant_CHROM, hub_variant_POS=mapsy_variant_table$hub_variant_POS, hub_variant_REF=mapsy_variant_table$hub_variant_REF, hub_variant_ALT=mapsy_variant_table$hub_variant_ALT,
			mpralm.ANCDER.logFC=mapsy_variant_table$mpralm.ANCDER.logFC)

	# get gene models
	gene_anno_gr <- as_tibble(fread("../../../final-mapsy-geisinger/data/known_canonical/gencode_v32_lift37_basic_canonical_exons.txt.gz")) %>% 
		mutate(seqnames = gsub("chr", "", seqnames)) %>% 
		mutate(SYMBOL=gene_name, transcript=gene_name) %>% 
		filter(seqnames == vis_CHROM) %>% 
		GRanges()

	# useful for defining windows
	gene_min <- min(start(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_max <- max(end(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_width <- gene_max - gene_min

	# load spliceai
	# # 	based on new SpliceAI scores
	# final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(hub_variant_CHROM == vis_CHROM) %>% 
		mutate(spliceai_max = spliceai_max)
	# if (("hub_in_final_study_adaptive" %in% temp_list) | ("hub_in_final_study_introgressed" %in% temp_list)) {
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(SYMBOL %in% c(vis_SYMBOL))
	# }
	final_v2_spliceai_filter <- final_v2_spliceai
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!is.na(spliceai_max))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!!rlang::parse_expr(temp_filter))

	final_v2_spliceai_gr <- 
		GRanges(seqnames=final_v2_spliceai$hub_variant_CHROM, IRanges(start=final_v2_spliceai$hub_variant_POS, end=final_v2_spliceai$hub_variant_POS), # strand=final_v2_spliceai$exon_strand, 
			hub_variant_ID=final_v2_spliceai$hub_variant_ID, hub_variant_CHROM=final_v2_spliceai$hub_variant_CHROM, hub_variant_POS=final_v2_spliceai$hub_variant_POS, hub_variant_REF=final_v2_spliceai$hub_variant_REF, hub_variant_ALT=final_v2_spliceai$hub_variant_ALT,
			spliceai_max=final_v2_spliceai$spliceai_max)

	# load LD pairs
	ld_pairs <- as_tibble(fread(paste("../../results/linkage_hapR2_SNPs/ALL.chr", vis_CHROM, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz", sep="")))
	ld_pairs <- ld_pairs %>% 
		filter(CHR == vis_CHROM, (POS1 == vis_POS)|(POS2 == vis_POS)) 
	ld_pairs_1 <- ld_pairs %>% 
		mutate(hub_variant_CHROM = as.character(CHR), hub_variant_POS = POS1, hub_variant_POS2 = POS2)
	ld_pairs_2 <- ld_pairs %>% 
		mutate(hub_variant_CHROM = as.character(CHR), hub_variant_POS = POS2, hub_variant_POS2 = POS1)
	ld_pairs_temp <- bind_rows(ld_pairs_1, ld_pairs_2) %>% 
		dplyr::select(-CHR, -POS1, -POS2)

	# load sQTLs
	list_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
	collate_GTEx_sQTLs_hg19 <- NULL
	for (tissue in list_tissues) {
		print(tissue)
		collate_GTEx_sQTLs_hg19_temp <- as_tibble(fread(paste("../../results/preprocess_GTEx_QTLs/sQTLs_tb_final_hg19-", tissue, ".txt.gz", sep=""))) %>% 
			filter(QTLs_gene_name %in% c(vis_SYMBOL))
		collate_GTEx_sQTLs_hg19[[paste(tissue, "sQTL", sep="_")]] <- collate_GTEx_sQTLs_hg19_temp
	}
	collate_GTEx_sQTLs_hg19_long <- collate_GTEx_sQTLs_hg19 %>% bind_rows() %>% 
		mutate(QTLs_tissue = gsub("_", " ", QTLs_tissue)) # %>% 
		# filter(hub_variant_ID %in% final_v2_spliceai_filter$hub_variant_ID)

	collate_GTEx_sQTLs_hg19_long_gr <- 
		GRanges(seqnames=collate_GTEx_sQTLs_hg19_long$hub_variant_CHROM, IRanges(start=collate_GTEx_sQTLs_hg19_long$hub_variant_POS, end=collate_GTEx_sQTLs_hg19_long$hub_variant_POS), # strand=collate_GTEx_sQTLs_hg19_long$exon_strand, 
			hub_variant_ID=collate_GTEx_sQTLs_hg19_long$hub_variant_ID, hub_variant_CHROM=collate_GTEx_sQTLs_hg19_long$hub_variant_CHROM, hub_variant_POS=collate_GTEx_sQTLs_hg19_long$hub_variant_POS, hub_variant_REF=collate_GTEx_sQTLs_hg19_long$hub_variant_REF, hub_variant_ALT=collate_GTEx_sQTLs_hg19_long$hub_variant_ALT,
			tissue=collate_GTEx_sQTLs_hg19_long$QTLs_tissue, pval_nominal=collate_GTEx_sQTLs_hg19_long$QTLs_pval_nominal, gene_name=collate_GTEx_sQTLs_hg19_long$QTLs_gene_name) %>% 
		as_tibble() %>% unique() %>% GRanges()

	# restrict by LD
	collate_GTEx_sQTLs_hg19_long_gr_ld <- collate_GTEx_sQTLs_hg19_long_gr %>% 
		as_tibble() %>% left_join(as_tibble(ld_pairs_temp)) %>% GRanges()

	# pivot wide tissue sQTLs
	collate_GTEx_sQTLs_hg19_long_gr_ld <- collate_GTEx_sQTLs_hg19_long_gr_ld %>% 
		mutate(pval_nominal = -log10(pval_nominal)) %>% 
		as_tibble() %>% pivot_wider(names_from=tissue, values_from=pval_nominal) %>% GRanges()

	# get shared tracks
	itrack <- IdeogramTrack(genome="hg19", chromosome=vis_CHROM, cex=1)
	gtrack <- GenomeAxisTrack(cex=1, fontsize=8)
	fcol <- c(A="darkgrey", C="darkgrey", G="darkgrey", T="darkgrey")
	# fcol[ifelse(vis_STRAND=="+", vis_REF, reverse_complement(vis_REF))] <- col_REF
	# fcol[ifelse(vis_STRAND=="+", vis_ALT, reverse_complement(vis_ALT))] <- col_ALT
	strack <- SequenceTrack(Hsapiens, chromosome=vis_CHROM, cex=0.7, fontcolor=fcol)
	grtrack <- GeneRegionTrack(gene_anno_gr, transcriptAnnotation="transcript", name="Gene", fontsize=12, fontsize.group=12)

	# plot sequence
	htrackseq <- HighlightTrack(trackList=list(strack, strack), start=c(vis_POS), width=0, chromosome=vis_CHROM, alpha=0.2)
	win_from <- vis_POS-win_half_size_seq
	win_to <- vis_POS+win_half_size_seq
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_seq.pdf", sep=""), 
		height=0.5, width=2.5)
	plotTracks(list( 
			htrackseq
		), 
		sizes=c(10,10),
		from=win_from, 
		to=win_to+1, 
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot ideogram
	ext_win_splice <- 1.5*max(abs(vis_POS-gene_min), abs(gene_max-vis_POS))
	win_from <- ifelse(!is.na(win_half_size_splice), vis_POS-ext_win_splice, vis_POS-ext_win_splice)
	win_to <- ifelse(!is.na(win_half_size_splice), vis_POS+ext_win_splice, vis_POS+ext_win_splice)

	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_ideo.pdf", sep=""), 
		height=0.5, width=4)
	plotTracks(list(
			itrack
		), 
		sizes=c(10),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# get data tracks
	# 	mapsy
	ylim_mapsy <- c(min(0, min(mapsy_variant_table_gr$mpralm.ANCDER.logFC)-0.2), max(0, max(mapsy_variant_table_gr$mpralm.ANCDER.logFC)+0.2))
	mapsytrack <- DataTrack(mapsy_variant_table_gr[,"mpralm.ANCDER.logFC"], type="h", name="MaPSy\nlog2 FC", col="#92c5de", fontsize=12, ylim=ylim_mapsy, cex=1)
	mapsytrack2 <- DataTrack(mapsy_variant_table_gr[,"mpralm.ANCDER.logFC"], name="MaPSy\nlog2 FC", col="#92c5de", fontsize=12, ylim=ylim_mapsy, cex=0.8, baseline=0, col.baseline="#92c5de", lty.baseline=2)
	mapsytrack3 <- DataTrack(filter(mapsy_variant_table_gr, start==vis_POS)[,"mpralm.ANCDER.logFC"], name="MaPSy log2 FC", col="#92c5de", fill="#FFFFFF", pch=23, cex=1.2, fontsize=12, ylim=ylim_mapsy)
	mapsyoverlay <- OverlayTrack(trackList=list(mapsytrack, mapsytrack2, mapsytrack3))

	# 	spliceai
	ylim_spliceai <- c(0, max(filter(final_v2_spliceai_gr, (as.numeric(hub_variant_POS) >= as.numeric(win_from)) & (as.numeric(hub_variant_POS) <= as.numeric(win_to)))$spliceai_max)+0.1)
	spliceaitrack <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", type="h", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=1)
	spliceaitrack2 <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=0.8, baseline=0, col.baseline="#3a6bb2", lty.baseline=2)
	if (length(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"]) > 0) {
		spliceaitrack3 <- DataTrack(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"], name="SpliceAI max score", col="#3a6bb2", fill="#FFFFFF",  pch=23, cex=1.2, fontsize=12, ylim=ylim_spliceai)
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2, spliceaitrack3))
	} else {
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2))
	}

	# plot exon splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_exon.pdf", sep=""), 
		height=0.75, width=4)
	plotTracks(list(
			htracksplice
		), 
		sizes=c(10),
		from=vis_POS-200,
		to=vis_POS+200,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot gene splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack, mapsyoverlay, spliceaioverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_splice.pdf", sep=""), 
		height=3, width=5)
	plotTracks(list(
			gtrack, 
			htracksplice
		), 
		sizes=c(10,10,20,20),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot gtex tracks if avail
	if (length(collate_GTEx_sQTLs_hg19_long_gr_ld)>=1) {
		tissues_avail <- intersect(gsub("_", " ", list_tissues), names(mcols(collate_GTEx_sQTLs_hg19_long_gr_ld)))
		ylim_sqtlld <- c(0, max(unlist(mcols(collate_GTEx_sQTLs_hg19_long_gr_ld[,tissues_avail])), na.rm=T)+1)
		sqtlldtrack <- DataTrack(collate_GTEx_sQTLs_hg19_long_gr_ld[,tissues_avail], name="GTEx sQTL\n-log10(p-value)", groups=tissues_avail, cex=0.6, fontsize=12, ylim=ylim_sqtlld, cex.legend=0.5, fontsize.legend=ifelse(length(tissues_avail)>8, ifelse(length(tissues_avail)>16, 8, 10), 12))
		sqtl_info <- filter(collate_GTEx_sQTLs_hg19_long_gr_ld, start==vis_POS)[,tissues_avail]
		if (length(sqtl_info) > 0) {
			sqtlldtrack2 <- DataTrack(sqtl_info, name="GTEx sQTL\n-log10(p-value)", groups=tissues_avail, fill="#FFFFFF", pch=23, cex=0.9, fontsize=12, ylim=ylim_sqtlld, cex.legend=0.5, fontsize.legend=ifelse(length(tissues_avail)>8, ifelse(length(tissues_avail)>16, 8, 10), 12))
			sqtlldoverlay <- OverlayTrack(trackList=list(sqtlldtrack, sqtlldtrack2))

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_gtex.pdf", sep=""), 
				height=3, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, spliceaioverlay, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_gtex.pdf", sep=""), 
				height=3.5, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,15,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, mapsyoverlay, spliceaioverlay, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_all.pdf", sep=""), 
				height=4, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,15,15,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()
		}
	}
}


visualize_genomic_range_manual_spliceai <- function(vis_variant_id, win_half_size_gtex=NA, win_half_size_splice=NA, win_half_size_seq=NA, out_folder, manual_spliceai=NA) {  # , col_REF="#FBBD7A", col_ALT="#377FA8") {
	# parse variant ID
	vis_variant_id_split <- strsplit(vis_variant_id, ":|_|/")[[1]]
	vis_CHROM <- gsub("chr", "", vis_variant_id_split[1])
	vis_POS <- as.numeric(vis_variant_id_split[2])
	vis_REF <- vis_variant_id_split[3]
	vis_ALT <- vis_variant_id_split[4]

	# load mapsy
	mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
	final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
	vis_SYMBOL1 <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_gene_name  # get gene name
	vis_SYMBOL2 <- final_v2_spliceai %>% filter(hub_variant_ID == vis_variant_id) %>% .$SYMBOL  # get gene name
	vis_SYMBOL <- unique(c(vis_SYMBOL1, vis_SYMBOL2))
	# vis_SYMBOL <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_gene_name  # get gene name
	vis_STRAND <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_strand  # get gene name
	mapsy_variant_table <- mapsy_variant_table %>%  # restrict to target gene
		filter(hub_variant_CHROM == vis_CHROM)
	temp <- mapsy_variant_table %>% 
		filter(hub_variant_ID == vis_variant_id) %>%  # .[,169:174]
		dplyr::select(hub_in_final_study_introgressed, hub_in_final_study_adaptive, hub_in_final_study_modern, hub_in_final_study_archaic, hub_in_final_study_nean, hub_in_final_study_deni)
	temp_list <- c()
	for (i in names(temp)) {
		if (temp[[i]][[1]]) {
			temp_list <- c(temp_list, i)
		}
	}
	# if (("hub_in_final_study_adaptive" %in% temp_list) | ("hub_in_final_study_introgressed" %in% temp_list)) {
	mapsy_variant_table <- mapsy_variant_table %>% 
		filter(exon_gene_name %in% c(vis_SYMBOL))
	# }
	temp_filter <- paste(temp_list[[1]], collapse="|")
	print(temp_filter)
	mapsy_variant_table <- mapsy_variant_table %>% 
		filter(!!rlang::parse_expr(temp_filter))

	# create mapsy grange
	mapsy_variant_table_gr <- 
		GRanges(seqnames=mapsy_variant_table$hub_variant_CHROM, IRanges(start=mapsy_variant_table$hub_variant_POS, end=mapsy_variant_table$hub_variant_POS), # strand=mapsy_variant_table$exon_strand, 
			hub_variant_ID=mapsy_variant_table$hub_variant_ID, hub_variant_CHROM=mapsy_variant_table$hub_variant_CHROM, hub_variant_POS=mapsy_variant_table$hub_variant_POS, hub_variant_REF=mapsy_variant_table$hub_variant_REF, hub_variant_ALT=mapsy_variant_table$hub_variant_ALT,
			mpralm.ANCDER.logFC=mapsy_variant_table$mpralm.ANCDER.logFC)


	# get gene models
	gene_anno_gr <- as_tibble(fread("../../../final-mapsy-geisinger/data/known_canonical/gencode_v32_lift37_basic_canonical_exons.txt.gz")) %>% 
		mutate(seqnames = gsub("chr", "", seqnames)) %>% 
		mutate(SYMBOL=gene_name, transcript=gene_name) %>% 
		filter(seqnames == vis_CHROM) %>% 
		GRanges()

	# useful for defining windows
	gene_min <- min(start(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_max <- max(end(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_width <- gene_max - gene_min

	# load spliceai
	# # 	based on new SpliceAI scores
	# final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(hub_variant_CHROM == vis_CHROM) %>% 
		mutate(spliceai_max = spliceai_max)
	# if (("hub_in_final_study_adaptive" %in% temp_list) | ("hub_in_final_study_introgressed" %in% temp_list)) {
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(SYMBOL %in% c(vis_SYMBOL))
	# }
	final_v2_spliceai_filter <- final_v2_spliceai
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!is.na(spliceai_max))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!!rlang::parse_expr(temp_filter))

	final_v2_spliceai_gr <- 
		GRanges(seqnames=final_v2_spliceai$hub_variant_CHROM, IRanges(start=final_v2_spliceai$hub_variant_POS, end=final_v2_spliceai$hub_variant_POS), # strand=final_v2_spliceai$exon_strand, 
			hub_variant_ID=final_v2_spliceai$hub_variant_ID, hub_variant_CHROM=final_v2_spliceai$hub_variant_CHROM, hub_variant_POS=final_v2_spliceai$hub_variant_POS, hub_variant_REF=final_v2_spliceai$hub_variant_REF, hub_variant_ALT=final_v2_spliceai$hub_variant_ALT,
			spliceai_max=final_v2_spliceai$spliceai_max)

	# manual repair spliceai value
	if (!is.na(manual_spliceai)) {
		if (length(final_v2_spliceai_gr %>% filter(hub_variant_ID==vis_variant_id))==0) {
			if (length(final_v2_spliceai_gr) > 0) {
				temp <- data.frame(final_v2_spliceai_gr)
				entry <- c(seqnames=vis_CHROM, start=vis_POS, end=vis_POS, width=1, strand="*", 
					hub_variant_ID=vis_variant_id, hub_variant_CHROM=vis_CHROM, hub_variant_POS=vis_POS, 
					hub_variant_REF=vis_REF, hub_variant_ALT=vis_ALT, spliceai_max=manual_spliceai)
				temp <- rbind(temp, entry) %>% as_tibble() %>% arrange(seqnames, start, end)
				final_v2_spliceai_gr <- GRanges(temp)	
			} else {
				entry <- c(seqnames=vis_CHROM, start=vis_POS, end=vis_POS, width=1, strand="*", 
					hub_variant_ID=vis_variant_id, hub_variant_CHROM=vis_CHROM, hub_variant_POS=vis_POS, 
					hub_variant_REF=vis_REF, hub_variant_ALT=vis_ALT, spliceai_max=manual_spliceai)
				temp <- rbind(entry) %>% as_tibble() %>% arrange(seqnames, start, end)
				final_v2_spliceai_gr <- GRanges(temp)
			}
		}
	}

	# load LD pairs
	ld_pairs <- as_tibble(fread(paste("../../results/linkage_hapR2_SNPs/ALL.chr", vis_CHROM, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz", sep="")))
	ld_pairs <- ld_pairs %>% 
		filter(CHR == vis_CHROM, (POS1 == vis_POS)|(POS2 == vis_POS)) 
	ld_pairs_1 <- ld_pairs %>% 
		mutate(hub_variant_CHROM = as.character(CHR), hub_variant_POS = POS1, hub_variant_POS2 = POS2)
	ld_pairs_2 <- ld_pairs %>% 
		mutate(hub_variant_CHROM = as.character(CHR), hub_variant_POS = POS2, hub_variant_POS2 = POS1)
	ld_pairs_temp <- bind_rows(ld_pairs_1, ld_pairs_2) %>% 
		dplyr::select(-CHR, -POS1, -POS2)

	# load sQTLs
	list_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
	collate_GTEx_sQTLs_hg19 <- NULL
	for (tissue in list_tissues) {
		print(tissue)
		collate_GTEx_sQTLs_hg19_temp <- as_tibble(fread(paste("../../results/preprocess_GTEx_QTLs/sQTLs_tb_final_hg19-", tissue, ".txt.gz", sep=""))) %>% 
			filter(QTLs_gene_name %in% c(vis_SYMBOL))
		collate_GTEx_sQTLs_hg19[[paste(tissue, "sQTL", sep="_")]] <- collate_GTEx_sQTLs_hg19_temp
	}
	collate_GTEx_sQTLs_hg19_long <- collate_GTEx_sQTLs_hg19 %>% bind_rows() %>% 
		mutate(QTLs_tissue = gsub("_", " ", QTLs_tissue)) # %>% 
		# filter(hub_variant_ID %in% final_v2_spliceai_filter$hub_variant_ID)

	collate_GTEx_sQTLs_hg19_long_gr <- 
		GRanges(seqnames=collate_GTEx_sQTLs_hg19_long$hub_variant_CHROM, IRanges(start=collate_GTEx_sQTLs_hg19_long$hub_variant_POS, end=collate_GTEx_sQTLs_hg19_long$hub_variant_POS), # strand=collate_GTEx_sQTLs_hg19_long$exon_strand, 
			hub_variant_ID=collate_GTEx_sQTLs_hg19_long$hub_variant_ID, hub_variant_CHROM=collate_GTEx_sQTLs_hg19_long$hub_variant_CHROM, hub_variant_POS=collate_GTEx_sQTLs_hg19_long$hub_variant_POS, hub_variant_REF=collate_GTEx_sQTLs_hg19_long$hub_variant_REF, hub_variant_ALT=collate_GTEx_sQTLs_hg19_long$hub_variant_ALT,
			tissue=collate_GTEx_sQTLs_hg19_long$QTLs_tissue, pval_nominal=collate_GTEx_sQTLs_hg19_long$QTLs_pval_nominal, gene_name=collate_GTEx_sQTLs_hg19_long$QTLs_gene_name) %>% 
		as_tibble() %>% unique() %>% GRanges()

	# restrict by LD
	collate_GTEx_sQTLs_hg19_long_gr_ld <- collate_GTEx_sQTLs_hg19_long_gr %>% 
		as_tibble() %>% left_join(as_tibble(ld_pairs_temp)) %>% GRanges()

	# pivot wide tissue sQTLs
	collate_GTEx_sQTLs_hg19_long_gr_ld <- collate_GTEx_sQTLs_hg19_long_gr_ld %>% 
		mutate(pval_nominal = -log10(pval_nominal)) %>% 
		as_tibble() %>% pivot_wider(names_from=tissue, values_from=pval_nominal) %>% GRanges()

	# get shared tracks
	itrack <- IdeogramTrack(genome="hg19", chromosome=vis_CHROM, cex=1)
	gtrack <- GenomeAxisTrack(cex=1, fontsize=8)
	fcol <- c(A="darkgrey", C="darkgrey", G="darkgrey", T="darkgrey")
	# fcol[ifelse(vis_STRAND=="+", vis_REF, reverse_complement(vis_REF))] <- col_REF
	# fcol[ifelse(vis_STRAND=="+", vis_ALT, reverse_complement(vis_ALT))] <- col_ALT
	strack <- SequenceTrack(Hsapiens, chromosome=vis_CHROM, cex=0.7, fontcolor=fcol)
	grtrack <- GeneRegionTrack(gene_anno_gr, transcriptAnnotation="transcript", name="Gene", fontsize=12, fontsize.group=12)

	# plot sequence
	htrackseq <- HighlightTrack(trackList=list(strack, strack), start=c(vis_POS), width=0, chromosome=vis_CHROM, alpha=0.2)
	win_from <- vis_POS-win_half_size_seq
	win_to <- vis_POS+win_half_size_seq
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_seq.pdf", sep=""), 
		height=0.5, width=2.5)
	plotTracks(list( 
			htrackseq
		), 
		sizes=c(10,10),
		from=win_from, 
		to=win_to+1, 
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot ideogram
	ext_win_splice <- 1.5*max(abs(vis_POS-gene_min), abs(gene_max-vis_POS))
	win_from <- ifelse(!is.na(win_half_size_splice), vis_POS-ext_win_splice, vis_POS-ext_win_splice)
	win_to <- ifelse(!is.na(win_half_size_splice), vis_POS+ext_win_splice, vis_POS+ext_win_splice)

	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_ideo.pdf", sep=""), 
		height=0.5, width=4)
	plotTracks(list(
			itrack
		), 
		sizes=c(10),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# get data tracks
	# 	mapsy
	ylim_mapsy <- c(min(0, min(mapsy_variant_table_gr$mpralm.ANCDER.logFC)-0.2), max(0, max(mapsy_variant_table_gr$mpralm.ANCDER.logFC)+0.2))
	mapsytrack <- DataTrack(mapsy_variant_table_gr[,"mpralm.ANCDER.logFC"], type="h", name="MaPSy\nlog2 FC", col="#92c5de", fontsize=12, ylim=ylim_mapsy, cex=1)
	mapsytrack2 <- DataTrack(mapsy_variant_table_gr[,"mpralm.ANCDER.logFC"], name="MaPSy\nlog2 FC", col="#92c5de", fontsize=12, ylim=ylim_mapsy, cex=0.8, baseline=0, col.baseline="#92c5de", lty.baseline=2)
	mapsytrack3 <- DataTrack(filter(mapsy_variant_table_gr, start==vis_POS)[,"mpralm.ANCDER.logFC"], name="MaPSy log2 FC", col="#92c5de", fill="#FFFFFF", pch=23, cex=1.2, fontsize=12, ylim=ylim_mapsy)
	mapsyoverlay <- OverlayTrack(trackList=list(mapsytrack, mapsytrack2, mapsytrack3))

	# 	spliceai
	ylim_spliceai <- c(0, 1)
	spliceaitrack <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", type="h", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=1)
	spliceaitrack2 <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=0.8, baseline=0, col.baseline="#3a6bb2", lty.baseline=2)
	if (length(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"]) > 0) {
		spliceaitrack3 <- DataTrack(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"], name="SpliceAI max score", col="#3a6bb2", fill="#FFFFFF",  pch=23, cex=1.2, fontsize=12, ylim=ylim_spliceai)
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2, spliceaitrack3))
	} else {
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2))
	}

	# plot exon splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_exon.pdf", sep=""), 
		height=0.75, width=4)
	plotTracks(list(
			htracksplice
		), 
		sizes=c(10),
		from=vis_POS-200,
		to=vis_POS+200,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot gene splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack, mapsyoverlay, spliceaioverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_splice.pdf", sep=""), 
		height=3, width=5)
	plotTracks(list(
			gtrack, 
			htracksplice
		), 
		sizes=c(10,10,20,20),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot gtex tracks if avail
	if (length(collate_GTEx_sQTLs_hg19_long_gr_ld)>=1) {
		tissues_avail <- intersect(gsub("_", " ", list_tissues), names(mcols(collate_GTEx_sQTLs_hg19_long_gr_ld)))
		ylim_sqtlld <- c(0, max(unlist(mcols(collate_GTEx_sQTLs_hg19_long_gr_ld[,tissues_avail])), na.rm=T)+1)
		sqtlldtrack <- DataTrack(collate_GTEx_sQTLs_hg19_long_gr_ld[,tissues_avail], name="GTEx sQTL\n-log10(p-value)", groups=tissues_avail, cex=0.6, fontsize=12, ylim=ylim_sqtlld, cex.legend=0.5, fontsize.legend=ifelse(length(tissues_avail)>8, ifelse(length(tissues_avail)>16, 8, 10), 12))
		sqtl_info <- filter(collate_GTEx_sQTLs_hg19_long_gr_ld, start==vis_POS)[,tissues_avail]
		if (length(sqtl_info) > 0) {
			sqtlldtrack2 <- DataTrack(sqtl_info, name="GTEx sQTL\n-log10(p-value)", groups=tissues_avail, fill="#FFFFFF", pch=23, cex=0.9, fontsize=12, ylim=ylim_sqtlld, cex.legend=0.5, fontsize.legend=ifelse(length(tissues_avail)>8, ifelse(length(tissues_avail)>16, 8, 10), 12))
			sqtlldoverlay <- OverlayTrack(trackList=list(sqtlldtrack, sqtlldtrack2))

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_gtex.pdf", sep=""), 
				height=3, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, spliceaioverlay, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_gtex.pdf", sep=""), 
				height=3.5, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,15,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, mapsyoverlay, spliceaioverlay, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_all.pdf", sep=""), 
				height=4, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,15,15,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()
		}
	}
}

visualize_genomic_range_mapsyless <- function(vis_variant_id, win_half_size_gtex=NA, win_half_size_splice=NA, win_half_size_seq=NA, out_folder) {  # , col_REF="#FBBD7A", col_ALT="#377FA8") {
	# parse variant ID
	vis_variant_id_split <- strsplit(vis_variant_id, ":|_|/")[[1]]
	vis_CHROM <- gsub("chr", "", vis_variant_id_split[1])
	vis_POS <- as.numeric(vis_variant_id_split[2])
	vis_REF <- vis_variant_id_split[3]
	vis_ALT <- vis_variant_id_split[4]

	# get gene models
	gene_anno_gr <- as_tibble(fread("../../../final-mapsy-geisinger/data/known_canonical/gencode_v32_lift37_basic_canonical_exons.txt.gz")) %>% 
		mutate(seqnames = gsub("chr", "", seqnames)) %>% 
		mutate(SYMBOL=gene_name, transcript=gene_name) %>% 
		filter(seqnames == vis_CHROM) %>% 
		GRanges()

	# load spliceai
	# # 	based on new SpliceAI scores
	final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(hub_variant_CHROM == vis_CHROM) %>% 
		mutate(spliceai_max = spliceai_max)
	vis_SYMBOL <- final_v2_spliceai %>% filter(hub_variant_ID == vis_variant_id) %>% .$SYMBOL
	vis_STRAND <- filter(gene_anno_gr, gene_name==vis_SYMBOL) %>% as_tibble() %>% .$strand %>% as.character() %>% .[1]
	# print(vis_STRAND)

	# useful for defining windows
	gene_min <- min(start(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_max <- max(end(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_width <- gene_max - gene_min

	temp <- final_v2_spliceai %>% 
		filter(hub_variant_ID == vis_variant_id) %>%  # .[,169:174]
		dplyr::select(hub_in_final_study_introgressed, hub_in_final_study_adaptive, hub_in_final_study_modern, hub_in_final_study_archaic, hub_in_final_study_nean, hub_in_final_study_deni)
	temp_list <- c()
	for (i in names(temp)) {
		if (temp[[i]][[1]]) {
			temp_list <- c(temp_list, i)
		}
	}
	temp_filter <- paste(temp_list[[1]], collapse="|")
	print(temp_filter)
	# if (("hub_in_final_study_adaptive" %in% temp_list) | ("hub_in_final_study_introgressed" %in% temp_list)) {
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(SYMBOL %in% c(vis_SYMBOL))
	# }
	final_v2_spliceai_filter <- final_v2_spliceai
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!is.na(spliceai_max))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!!rlang::parse_expr(temp_filter))

	final_v2_spliceai_gr <- 
		GRanges(seqnames=final_v2_spliceai$hub_variant_CHROM, IRanges(start=final_v2_spliceai$hub_variant_POS, end=final_v2_spliceai$hub_variant_POS), # strand=final_v2_spliceai$exon_strand, 
			hub_variant_ID=final_v2_spliceai$hub_variant_ID, hub_variant_CHROM=final_v2_spliceai$hub_variant_CHROM, hub_variant_POS=final_v2_spliceai$hub_variant_POS, hub_variant_REF=final_v2_spliceai$hub_variant_REF, hub_variant_ALT=final_v2_spliceai$hub_variant_ALT,
			spliceai_max=final_v2_spliceai$spliceai_max)

	# load LD pairs
	ld_pairs <- as_tibble(fread(paste("../../results/linkage_hapR2_SNPs/ALL.chr", vis_CHROM, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz", sep="")))
	ld_pairs <- ld_pairs %>% 
		filter(CHR == vis_CHROM, (POS1 == vis_POS)|(POS2 == vis_POS)) 
	ld_pairs_1 <- ld_pairs %>% 
		mutate(hub_variant_CHROM = as.character(CHR), hub_variant_POS = POS1, hub_variant_POS2 = POS2)
	ld_pairs_2 <- ld_pairs %>% 
		mutate(hub_variant_CHROM = as.character(CHR), hub_variant_POS = POS2, hub_variant_POS2 = POS1)
	ld_pairs_temp <- bind_rows(ld_pairs_1, ld_pairs_2) %>% 
		dplyr::select(-CHR, -POS1, -POS2)

	# load sQTLs
	list_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
	collate_GTEx_sQTLs_hg19 <- NULL
	for (tissue in list_tissues) {
		print(tissue)
		collate_GTEx_sQTLs_hg19_temp <- as_tibble(fread(paste("../../results/preprocess_GTEx_QTLs/sQTLs_tb_final_hg19-", tissue, ".txt.gz", sep=""))) %>% 
			filter(QTLs_gene_name %in% c(vis_SYMBOL))
		collate_GTEx_sQTLs_hg19[[paste(tissue, "sQTL", sep="_")]] <- collate_GTEx_sQTLs_hg19_temp
	}
	collate_GTEx_sQTLs_hg19_long <- collate_GTEx_sQTLs_hg19 %>% bind_rows() %>% 
		mutate(QTLs_tissue = gsub("_", " ", QTLs_tissue)) # %>% 
		# filter(hub_variant_ID %in% final_v2_spliceai_filter$hub_variant_ID)

	collate_GTEx_sQTLs_hg19_long_gr <- 
		GRanges(seqnames=collate_GTEx_sQTLs_hg19_long$hub_variant_CHROM, IRanges(start=collate_GTEx_sQTLs_hg19_long$hub_variant_POS, end=collate_GTEx_sQTLs_hg19_long$hub_variant_POS), # strand=collate_GTEx_sQTLs_hg19_long$exon_strand, 
			hub_variant_ID=collate_GTEx_sQTLs_hg19_long$hub_variant_ID, hub_variant_CHROM=collate_GTEx_sQTLs_hg19_long$hub_variant_CHROM, hub_variant_POS=collate_GTEx_sQTLs_hg19_long$hub_variant_POS, hub_variant_REF=collate_GTEx_sQTLs_hg19_long$hub_variant_REF, hub_variant_ALT=collate_GTEx_sQTLs_hg19_long$hub_variant_ALT,
			tissue=collate_GTEx_sQTLs_hg19_long$QTLs_tissue, pval_nominal=collate_GTEx_sQTLs_hg19_long$QTLs_pval_nominal, gene_name=collate_GTEx_sQTLs_hg19_long$QTLs_gene_name) %>% 
		as_tibble() %>% unique() %>% GRanges()

	# restrict by LD
	collate_GTEx_sQTLs_hg19_long_gr_ld <- collate_GTEx_sQTLs_hg19_long_gr %>% 
		as_tibble() %>% left_join(as_tibble(ld_pairs_temp)) %>% GRanges()

	# pivot wide tissue sQTLs
	collate_GTEx_sQTLs_hg19_long_gr_ld <- collate_GTEx_sQTLs_hg19_long_gr_ld %>% 
		mutate(pval_nominal = -log10(pval_nominal)) %>% 
		as_tibble() %>% pivot_wider(names_from=tissue, values_from=pval_nominal) %>% GRanges()

	# get shared tracks
	itrack <- IdeogramTrack(genome="hg19", chromosome=vis_CHROM, cex=1)
	gtrack <- GenomeAxisTrack(cex=1, fontsize=8)
	fcol <- c(A="darkgrey", C="darkgrey", G="darkgrey", T="darkgrey")
	# fcol[ifelse(vis_STRAND=="+", vis_REF, reverse_complement(vis_REF))] <- col_REF
	# fcol[ifelse(vis_STRAND=="+", vis_ALT, reverse_complement(vis_ALT))] <- col_ALT
	strack <- SequenceTrack(Hsapiens, chromosome=vis_CHROM, cex=0.7, fontcolor=fcol)
	grtrack <- GeneRegionTrack(gene_anno_gr, transcriptAnnotation="transcript", name="Gene", fontsize=12, fontsize.group=12)

	# plot sequence
	htrackseq <- HighlightTrack(trackList=list(strack, strack), start=c(vis_POS), width=0, chromosome=vis_CHROM, alpha=0.2)
	win_from <- vis_POS-win_half_size_seq
	win_to <- vis_POS+win_half_size_seq
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_seq.pdf", sep=""), 
		height=0.5, width=2.5)
	plotTracks(list( 
			htrackseq
		), 
		sizes=c(10,10),
		from=win_from, 
		to=win_to+1, 
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot ideogram
	ext_win_splice <- 1.5*max(abs(vis_POS-gene_min), abs(gene_max-vis_POS))
	win_from <- ifelse(!is.na(win_half_size_splice), vis_POS-ext_win_splice, vis_POS-ext_win_splice)
	win_to <- ifelse(!is.na(win_half_size_splice), vis_POS+ext_win_splice, vis_POS+ext_win_splice)

	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_ideo.pdf", sep=""), 
		height=0.5, width=4)
	plotTracks(list(
			itrack
		), 
		sizes=c(10),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# get data tracks
	# 	spliceai
	ylim_spliceai <- c(0, max(filter(final_v2_spliceai_gr, (as.numeric(hub_variant_POS) >= as.numeric(win_from)) & (as.numeric(hub_variant_POS) <= as.numeric(win_to)))$spliceai_max)+0.1)
	spliceaitrack <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", type="h", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=1)
	spliceaitrack2 <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=0.8, baseline=0, col.baseline="#3a6bb2", lty.baseline=2)
	if (length(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"]) > 0) {
		spliceaitrack3 <- DataTrack(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"], name="SpliceAI max score", col="#3a6bb2", fill="#FFFFFF",  pch=23, cex=1.2, fontsize=12, ylim=ylim_spliceai)
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2, spliceaitrack3))
	} else {
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2))
	}

	# plot exon splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_exon.pdf", sep=""), 
		height=0.75, width=4)
	plotTracks(list(
			htracksplice
		), 
		sizes=c(10),
		from=vis_POS-200,
		to=vis_POS+200,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot gene splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack, spliceaioverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_splice.pdf", sep=""), 
		height=3, width=5)
	plotTracks(list(
			gtrack, 
			htracksplice
		), 
		sizes=c(10,10,20),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot gtex tracks if avail
	if (length(collate_GTEx_sQTLs_hg19_long_gr_ld)>=1) {
		tissues_avail <- intersect(gsub("_", " ", list_tissues), names(mcols(collate_GTEx_sQTLs_hg19_long_gr_ld)))
		ylim_sqtlld <- c(0, max(unlist(mcols(collate_GTEx_sQTLs_hg19_long_gr_ld[,tissues_avail])), na.rm=T)+1)
		sqtlldtrack <- DataTrack(collate_GTEx_sQTLs_hg19_long_gr_ld[,tissues_avail], name="GTEx sQTL\n-log10(p-value)", groups=tissues_avail, cex=0.6, fontsize=12, ylim=ylim_sqtlld, cex.legend=0.5, fontsize.legend=ifelse(length(tissues_avail)>8, ifelse(length(tissues_avail)>16, 8, 10), 12))
		sqtl_info <- filter(collate_GTEx_sQTLs_hg19_long_gr_ld, start==vis_POS)[,tissues_avail]
		if (length(sqtl_info) > 0) {
			sqtlldtrack2 <- DataTrack(sqtl_info, name="GTEx sQTL\n-log10(p-value)", groups=tissues_avail, fill="#FFFFFF", pch=23, cex=0.9, fontsize=12, ylim=ylim_sqtlld, cex.legend=0.5, fontsize.legend=ifelse(length(tissues_avail)>8, ifelse(length(tissues_avail)>16, 8, 10), 12))
			sqtlldoverlay <- OverlayTrack(trackList=list(sqtlldtrack, sqtlldtrack2))

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_gtex.pdf", sep=""), 
				height=2.5, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, spliceaioverlay, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_gtex.pdf", sep=""), 
				height=3, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,15,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()

			# 	highlight
			htrackgtex <- HighlightTrack(trackList=list(grtrack, spliceaioverlay, sqtlldoverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")

			# plot tracks
			win_from <- ifelse(!is.na(win_half_size_gtex), vis_POS-win_half_size_gtex, min(start(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			win_to <- ifelse(!is.na(win_half_size_gtex), vis_POS+win_half_size_gtex, max(end(collate_GTEx_sQTLs_hg19_long_gr_ld)))
			pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_all.pdf", sep=""), 
				height=3.5, width=5)
			plotTracks(list(
					gtrack, 
					htrackgtex
				), 
				sizes=c(10,10,15,30),
				from=win_from,
				to=win_to+1,
				reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
				complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
			)
			dev.off()
		}
	}
}

visualize_genomic_range_no_LD <- function(vis_variant_id, win_half_size_gtex=NA, win_half_size_splice=NA, win_half_size_seq=NA, out_folder) {  # , col_REF="#FBBD7A", col_ALT="#377FA8") {
	# parse variant ID
	vis_variant_id_split <- strsplit(vis_variant_id, ":|_|/")[[1]]
	vis_CHROM <- gsub("chr", "", vis_variant_id_split[1])
	vis_POS <- as.numeric(vis_variant_id_split[2])
	vis_REF <- vis_variant_id_split[3]
	vis_ALT <- vis_variant_id_split[4]

	# load mapsy
	mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
	final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
	vis_SYMBOL1 <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_gene_name  # get gene name
	vis_SYMBOL2 <- final_v2_spliceai %>% filter(hub_variant_ID == vis_variant_id) %>% .$SYMBOL  # get gene name
	vis_SYMBOL <- unique(c(vis_SYMBOL1, vis_SYMBOL2))
	# vis_SYMBOL <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_gene_name  # get gene name
	vis_STRAND <- mapsy_variant_table %>% filter(hub_variant_ID == vis_variant_id) %>% .$exon_strand  # get gene name
	mapsy_variant_table <- mapsy_variant_table %>%  # restrict to target gene
		filter(hub_variant_CHROM == vis_CHROM)
	temp <- mapsy_variant_table %>% 
		filter(hub_variant_ID == vis_variant_id) %>%  # .[,169:174]
		dplyr::select(hub_in_final_study_introgressed, hub_in_final_study_adaptive, hub_in_final_study_modern, hub_in_final_study_archaic, hub_in_final_study_nean, hub_in_final_study_deni)
	temp_list <- c()
	for (i in names(temp)) {
		if (temp[[i]][[1]]) {
			temp_list <- c(temp_list, i)
		}
	}
	# if (("hub_in_final_study_adaptive" %in% temp_list) | ("hub_in_final_study_introgressed" %in% temp_list)) {
	mapsy_variant_table <- mapsy_variant_table %>% 
		filter(exon_gene_name %in% c(vis_SYMBOL))
	# }
	temp_filter <- paste(temp_list[[1]], collapse="|")
	print(temp_filter)
	mapsy_variant_table <- mapsy_variant_table %>% 
		filter(!!rlang::parse_expr(temp_filter))

	# create mapsy grange
	mapsy_variant_table_gr <- 
		GRanges(seqnames=mapsy_variant_table$hub_variant_CHROM, IRanges(start=mapsy_variant_table$hub_variant_POS, end=mapsy_variant_table$hub_variant_POS), # strand=mapsy_variant_table$exon_strand, 
			hub_variant_ID=mapsy_variant_table$hub_variant_ID, hub_variant_CHROM=mapsy_variant_table$hub_variant_CHROM, hub_variant_POS=mapsy_variant_table$hub_variant_POS, hub_variant_REF=mapsy_variant_table$hub_variant_REF, hub_variant_ALT=mapsy_variant_table$hub_variant_ALT,
			mpralm.ANCDER.logFC=mapsy_variant_table$mpralm.ANCDER.logFC)


	# get gene models
	gene_anno_gr <- as_tibble(fread("../../../final-mapsy-geisinger/data/known_canonical/gencode_v32_lift37_basic_canonical_exons.txt.gz")) %>% 
		mutate(seqnames = gsub("chr", "", seqnames)) %>% 
		mutate(SYMBOL=gene_name, transcript=gene_name) %>% 
		filter(seqnames == vis_CHROM) %>% 
		GRanges()

	# useful for defining windows
	gene_min <- min(start(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_max <- max(end(filter(gene_anno_gr, gene_name==vis_SYMBOL)))
	gene_width <- gene_max - gene_min

	# load spliceai
	# # 	based on new SpliceAI scores
	# final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(hub_variant_CHROM == vis_CHROM) %>% 
		mutate(spliceai_max = spliceai_max)
	# if (("hub_in_final_study_adaptive" %in% temp_list) | ("hub_in_final_study_introgressed" %in% temp_list)) {
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(SYMBOL %in% c(vis_SYMBOL))
	# }
	final_v2_spliceai_filter <- final_v2_spliceai
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!is.na(spliceai_max))
	final_v2_spliceai <- final_v2_spliceai %>% 
		filter(!!rlang::parse_expr(temp_filter))

	final_v2_spliceai_gr <- 
		GRanges(seqnames=final_v2_spliceai$hub_variant_CHROM, IRanges(start=final_v2_spliceai$hub_variant_POS, end=final_v2_spliceai$hub_variant_POS), # strand=final_v2_spliceai$exon_strand, 
			hub_variant_ID=final_v2_spliceai$hub_variant_ID, hub_variant_CHROM=final_v2_spliceai$hub_variant_CHROM, hub_variant_POS=final_v2_spliceai$hub_variant_POS, hub_variant_REF=final_v2_spliceai$hub_variant_REF, hub_variant_ALT=final_v2_spliceai$hub_variant_ALT,
			spliceai_max=final_v2_spliceai$spliceai_max)

	# get shared tracks
	itrack <- IdeogramTrack(genome="hg19", chromosome=vis_CHROM, cex=1)
	gtrack <- GenomeAxisTrack(cex=1, fontsize=8)
	fcol <- c(A="darkgrey", C="darkgrey", G="darkgrey", T="darkgrey")
	# fcol[ifelse(vis_STRAND=="+", vis_REF, reverse_complement(vis_REF))] <- col_REF
	# fcol[ifelse(vis_STRAND=="+", vis_ALT, reverse_complement(vis_ALT))] <- col_ALT
	strack <- SequenceTrack(Hsapiens, chromosome=vis_CHROM, cex=0.7, fontcolor=fcol)
	grtrack <- GeneRegionTrack(gene_anno_gr, transcriptAnnotation="transcript", name="Gene", fontsize=12, fontsize.group=12)

	# plot sequence
	htrackseq <- HighlightTrack(trackList=list(strack, strack), start=c(vis_POS), width=0, chromosome=vis_CHROM, alpha=0.2)
	win_from <- vis_POS-win_half_size_seq
	win_to <- vis_POS+win_half_size_seq
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_seq.pdf", sep=""), 
		height=0.5, width=2.5)
	plotTracks(list( 
			htrackseq
		), 
		sizes=c(10,10),
		from=win_from, 
		to=win_to+1, 
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot ideogram
	ext_win_splice <- 1.5*max(abs(vis_POS-gene_min), abs(gene_max-vis_POS))
	win_from <- ifelse(!is.na(win_half_size_splice), vis_POS-ext_win_splice, vis_POS-ext_win_splice)
	win_to <- ifelse(!is.na(win_half_size_splice), vis_POS+ext_win_splice, vis_POS+ext_win_splice)

	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_ideo.pdf", sep=""), 
		height=0.5, width=4)
	plotTracks(list(
			itrack
		), 
		sizes=c(10),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# get data tracks
	# 	mapsy
	ylim_mapsy <- c(min(0, min(mapsy_variant_table_gr$mpralm.ANCDER.logFC)-0.2), max(0, max(mapsy_variant_table_gr$mpralm.ANCDER.logFC)+0.2))
	mapsytrack <- DataTrack(mapsy_variant_table_gr[,"mpralm.ANCDER.logFC"], type="h", name="MaPSy\nlog2 FC", col="#92c5de", fontsize=12, ylim=ylim_mapsy, cex=1)
	mapsytrack2 <- DataTrack(mapsy_variant_table_gr[,"mpralm.ANCDER.logFC"], name="MaPSy\nlog2 FC", col="#92c5de", fontsize=12, ylim=ylim_mapsy, cex=0.8, baseline=0, col.baseline="#92c5de", lty.baseline=2)
	mapsytrack3 <- DataTrack(filter(mapsy_variant_table_gr, start==vis_POS)[,"mpralm.ANCDER.logFC"], name="MaPSy log2 FC", col="#92c5de", fill="#FFFFFF", pch=23, cex=1.2, fontsize=12, ylim=ylim_mapsy)
	mapsyoverlay <- OverlayTrack(trackList=list(mapsytrack, mapsytrack2, mapsytrack3))

	# 	spliceai
	ylim_spliceai <- c(0, max(filter(final_v2_spliceai_gr, (as.numeric(hub_variant_POS) >= as.numeric(win_from)) & (as.numeric(hub_variant_POS) <= as.numeric(win_to)))$spliceai_max)+0.1)
	spliceaitrack <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", type="h", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=1)
	spliceaitrack2 <- DataTrack(final_v2_spliceai_gr[,"spliceai_max"], name="SpliceAI\nmax score", col="#3a6bb2", fontsize=12, ylim=ylim_spliceai, cex=0.8, baseline=0, col.baseline="#3a6bb2", lty.baseline=2)
	if (length(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"]) > 0) {
		spliceaitrack3 <- DataTrack(filter(final_v2_spliceai_gr, start==vis_POS)[,"spliceai_max"], name="SpliceAI max score", col="#3a6bb2", fill="#FFFFFF",  pch=23, cex=1.2, fontsize=12, ylim=ylim_spliceai)
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2, spliceaitrack3))
	} else {
		spliceaioverlay <- OverlayTrack(trackList=list(spliceaitrack, spliceaitrack2))
	}

	# plot exon splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_exon.pdf", sep=""), 
		height=0.75, width=4)
	plotTracks(list(
			htracksplice
		), 
		sizes=c(10),
		from=vis_POS-200,
		to=vis_POS+200,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()

	# plot gene splicing tracks
	htracksplice <- HighlightTrack(trackList=list(grtrack, mapsyoverlay, spliceaioverlay), start=vis_POS, end=vis_POS, chromosome=vis_CHROM, alpha=0.1, col="#000000")
	pdf(paste(out_folder, vis_SYMBOL, "_", gsub("/", "_", vis_variant_id), "_splice.pdf", sep=""), 
		height=3, width=5)
	plotTracks(list(
			gtrack, 
			htracksplice
		), 
		sizes=c(10,10,20,20),
		from=win_from,
		to=win_to+1,
		reverseStrand=ifelse(vis_STRAND=="+", FALSE, TRUE),
		complement=ifelse(vis_STRAND=="+", FALSE, TRUE)
	)
	dev.off()
}
