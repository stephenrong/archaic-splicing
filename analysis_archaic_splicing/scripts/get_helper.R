#!/bin/R

# load libraries
library(vcfR)
library(wrapr)
library(tidyverse)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.HsapiensGRCh37ancestor.Ensembl75.custom)
library(BSgenome.HsapiensGRCh38ancestor.Ensembl99.custom)

hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg38 <- BSgenome.Hsapiens.UCSC.hg38
hg19anc <- BSgenome.HsapiensGRCh37ancestor.Ensembl75.custom
hg38anc <- BSgenome.HsapiensGRCh38ancestor.Ensembl99.custom

# helper functions
extract_hub <- function(t) {
	t[which(grepl("hub_", names(t)))]
}

# mutate functions
mutate_hub_variant_REF <- function(input) {
	reference <- c("hg19" = hg19, "hg38" = hg38)[[input$hub_reference_genome[[1]]]]
	input %>% mutate(hub_variant_REF = toupper(as.character(getSeq(reference, 
		paste("chr", hub_variant_CHROM, sep=""), 
		hub_variant_POS, hub_variant_POS))))
}

mutate_hub_variant_ANC <- function(input) {
	ancestral <- c("hg19" = hg19anc, "hg38" = hg38anc)[[input$hub_reference_genome[[1]]]]
	input %>% mutate(hub_variant_ANC = toupper(as.character(getSeq(ancestral, 
		paste("chr", hub_variant_CHROM, sep=""), 
		hub_variant_POS, hub_variant_POS))))
}

mutate_hub_variant_DER <- function(input) {
	input %>% mutate(hub_variant_DER = ifelse(
		(hub_variant_REF == hub_variant_ANC) & 
		(hub_variant_ALT != hub_variant_ANC), 
		hub_variant_ALT, ifelse(
			(hub_variant_REF != hub_variant_ANC) & 
			(hub_variant_ALT == hub_variant_ANC),
			hub_variant_REF, "N")))
}

mutate_hub_variant_ID <- function(input) {
	input %>% mutate(hub_variant_ID = paste(
		hub_variant_CHROM, "_", 
		format(hub_variant_POS, scientific=FALSE, trim=TRUE), "_", 
		hub_variant_REF, "/", hub_variant_ALT, sep=""))
}

select_hub_variant <- function(input) {
	input %>% dplyr::select((names(.)[which(grepl("hub_", names(.)))]))
}

# expand out ancestral  # actually don't
split_ancestral <- function(input) {
	# consider thre cases
	temp1 <- input %>% 
		dplyr::filter(hub_variant_ANC == hub_variant_REF) %>% 
		mutate(hub_variant_CASE = 1)
	temp2 <- input %>% 
		dplyr::filter(!(hub_variant_ANC == hub_variant_REF) & 
			(hub_variant_ANC == hub_variant_ALT)) %>% 
		mutate(hub_variant_CASE = 2)
	temp3 <- input %>% 
		dplyr::filter(!(hub_variant_ANC == hub_variant_REF) & 
			!(hub_variant_ANC == hub_variant_ALT)) %>% 
		mutate(hub_variant_CASE = 3)
	# # split off if anc is different from alt and ref
	# temp4 <- temp3 %>% 
	# 	mutate(
	# 		hub_variant_ALT = hub_variant_ANC, 
	# 		hub_variant_CASE = 4)
	# bind and correct hub variant id
	bind_rows(temp1, temp2, temp3) %>%  # , temp4) %>% 
		mutate_hub_variant_ID()
}

drop_rows_all_na <- function(input) {
	# https://stackoverflow.com/questions/41609912/
	# remove-rows-where-all-variables-are-na-using-dplyr
	input %>% filter_all(any_vars(!is.na(.)))
}

standard_sort <- function(input) {
	input %>% arrange(hub_variant_CHROM, hub_variant_POS, 
		hub_variant_REF, hub_variant_ALT)
}

clean_up_hub <- function(input) {
	input %>% drop_rows_all_na() %>% 
		standard_sort() %>% unique()
}

tab <- function(input) {
	# print("hub_variant_REF")
	print(table(input["hub_variant_REF"]))
	# print("hub_variant_ALT")
	print(table(input["hub_variant_ALT"]))
	# print("hub_variant_ANC")
	print(table(input["hub_variant_ANC"]))
	# print("hub_variant_DER")
	print(table(input["hub_variant_DER"]))
	# print("hub_variant_REF == hub_variant_ALT")
	print(table(input["hub_variant_REF"] == input["hub_variant_ALT"]))
	# print("hub_variant_REF == hub_variant_ANC")
	print(table(input["hub_variant_REF"] == input["hub_variant_ANC"]))
	# print("hub_variant_ANC == hub_variant_DER")
	print(table(input["hub_variant_ANC"] == input["hub_variant_DER"]))
	# print("hub_variant_ALT == hub_variant_DER")
	print(table(input["hub_variant_ALT"] == input["hub_variant_DER"]))
}

# move hub to vcf function here
hub_to_vcf <- function(hub_tb, hub_tb_name, hub_tb_outfile, 
	hub_reference_genome, partial=TRUE, compress=TRUE) {
	# vcf_fixed
	vcf_fixed <- hub_tb %>% dplyr::select(
		hub_variant_CHROM, hub_variant_POS, 
		hub_variant_ID, hub_variant_REF, hub_variant_ALT)
	if (partial) {
		vcf_info <- hub_tb %>% 
			dplyr::select(hub_variant_ANC, hub_variant_DER, hub_variant_CASE)
	} else {
		vcf_info <- hub_tb %>% 
			dplyr::select(setdiff(names(hub_tb), 
				c(names(vcf_fixed), "hub_reference_genome")))
	}
	names(vcf_fixed) <- c("#CHROM", "POS", "ID", "REF", "ALT")
	vcf_fixed <- vcf_fixed %>% mutate(QUAL = ".", FILTER = ".", INFO = ".")

	# vcf_info
	vcf_info_temp <- vcf_info
	names(vcf_info_temp) <- gsub("hub_variant_", "", names(vcf_info_temp))
	for (i in names(vcf_info_temp)) {
		# let(c(VAR=i), vcf_info_temp <- vcf_info_temp %>% 
		# 	mutate(VAR = paste(i, VAR, sep="=")))
		vcf_info_temp[[i]] <- paste(i, vcf_info_temp[[i]], sep="=")
	}
	vcf_fixed$INFO <- vcf_info_temp %>% 
		unite(col=INFO, sep=";") %>% .$INFO

	vcf_fixed <- vcf_fixed %>% 
		arrange(`#CHROM`, POS, REF, ALT)

	# vcf_contig
	vcf_contig_temp <- as_tibble(rownames_to_column(as.data.frame(
		seqinfo(hub_reference_genome)), "seqname"))
	vcf_contig_temp <- vcf_contig_temp %>% 
		mutate(seqname = gsub("chr", "", seqname)) %>% 
		dplyr::filter(seqname %in% vcf_fixed$`#CHROM`)
	vcf_contig_temp <- vcf_contig_temp %>% 
		mutate(contig_meta = paste(
			"##contig=<", 
			"ID=", seqname, 
			",length=", seqlengths, 
			",assembly=", attributes(hub_reference_genome)$provider_version,
			",species=\"", attributes(hub_reference_genome)$organism, "\"",
			">", sep=""))
	vcf_contig <- vcf_contig_temp$contig_meta

	# vcf_meta
	# 	map r types to vcf format types
	map <- c("logical"="Integer", "integer"="Integer", "double"="Float","character"="String", 
		"complex"="String", "raw"="String", "list"="String", "factor"="String", 
		"ordered"="String", "Date"="String", "POSIXt"="String", "difftime"="String")
	vcf_meta <- c(
		"##fileformat=VCFv4.3",
		paste("##fileDate=", gsub("-", "", Sys.Date()), sep=""), 
		paste("##source=", hub_tb_name, sep=""),
		paste("##reference=", 
			attributes(hub_reference_genome)$release_name, sep=""),
		vcf_contig,
		as.vector(sapply(names(vcf_info_temp), function(x) {
			paste("##INFO=<ID=", x, ",Number=1,Type=", 
				map[class(vcf_info_temp[[x]])], 
				",Description=\"", x, "\">", sep="")
		}))
	)

	# output vcf file
	writeLines(vcf_meta, hub_tb_outfile)
	write_tsv(vcf_fixed, hub_tb_outfile, append=TRUE, col_names=TRUE)
	if (compress) {
		system(paste("bgzip -f ", hub_tb_outfile, sep=""))
		system(paste("tabix -p vcf ", hub_tb_outfile, ".gz", sep=""))
	}
}

hub_to_bed <- function(hub_tb, hub_tb_outfile) {
	hub_bed <- hub_tb %>% 
		mutate(
			bed_CHROM = hub_variant_CHROM, 
			bed_START = hub_variant_POS-1, 
			bed_END = hub_variant_POS
		) %>% 
		dplyr::select(bed_CHROM, bed_START, bed_END)
	write.table(hub_bed, file=hub_tb_outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}
