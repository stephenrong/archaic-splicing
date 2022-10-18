#!/usr/bin/env Rscript

# Helper functions for automating oligo sequence design.

# load packages
library(tidyverse)
library(plyranges)
library(BSgenome)
library(rtracklayer)
library(GenomicFeatures)
library(VariantAnnotation)
library(stringdist)
library(stringi)

# extend grange elements
extend_granges <- function(x, upstream=0, downstream=0) {
  # based on https://support.bioconductor.org/p/78652/
  if (any(strand(x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  return(trim(x))
}

# load and process fasta file
file_fasta_seqlevelsStyle <- function(loaded_fasta, style="UCSC") {
  seqlevelsStyle(names(loaded_fasta)) <- style
  return(loaded_fasta)
}

get_fasta <- function(file_fasta) {
  loaded_fasta <- rtracklayer::import(file_fasta, format="fasta")
  names(loaded_fasta) <- gsub(" .*", "", names(loaded_fasta))
  loaded_fasta <- file_fasta_seqlevelsStyle(loaded_fasta)
  return(loaded_fasta)  # DNAStringSet
}

# load and process gff3 file
loaded_gr_seqlevelsStyle <- function(loaded_gr, style="UCSC") {
  seqlevelsStyle(loaded_gr) <- style
  return(loaded_gr)
}

# standardize gff grange into exons grange
extract_exons <- function(loaded_gff, mode="gencode") {
  # filter protein coding genomic features
  if (mode == "gencode") {
    # handle gencode format gff files
    suppressWarnings(loaded_txdb <- makeTxDbFromGRanges(loaded_gff))
    loaded_exons_grl <- exonsBy(loaded_txdb, by="tx") %>% unlist()
    loaded_gff_exon <- loaded_gff %>% 
      as_tibble() %>% 
      filter(type == "exon") %>% 
      mutate(exon_name = gsub("exon:", "", ID)) %>% 
      dplyr::rename(exon_id_gencode = exon_id)
    loaded_exons_grl_joined <- loaded_exons_grl %>% 
      as_tibble() %>% left_join(loaded_gff_exon, 
        by=c("seqnames", "start", "end", "width", "strand", "exon_name")) %>% GRanges()
  } else if (mode == "ensembl") {
    # handle ensembl format gff files
    suppressWarnings(loaded_txdb <- makeTxDbFromGRanges(loaded_gff))
    loaded_exons_grl <- exonsBy(loaded_txdb, by="tx") %>% unlist()
    loaded_gff_exon <- loaded_gff %>% 
      as_tibble() %>% 
      filter(type == "exon") %>% 
      mutate(exon_name = Name) %>% 
      dplyr::rename(exon_id_ensembl = exon_id)
    loaded_exons_grl_joined <- loaded_exons_grl %>% 
      as_tibble() %>% left_join(loaded_gff_exon,
        by=c("seqnames", "start", "end", "width", "strand", "exon_name")) %>% GRanges()
  }
  return(loaded_exons_grl_joined)
}

get_gff <- function(file_gff, mode="gencode") {
  loaded_gff <- rtracklayer::import(file_gff, format="gff")
  # filter protein coding genomic features
  if (mode == "gencode") {
    loaded_gff <- loaded_gff %>% 
      filter(gene_type == "protein_coding") # %>% 
    # filter(seqnames != "chrM")  # remove mitochondrion
  } else if (mode == "ensembl") {
    loaded_gff <- loaded_gff %>%  # modify ids for convenience
      as_tibble() %>% mutate(ID = ifelse(is.na(ID), Parent, ID)) %>% 
      fill(biotype, gene_id, transcript_id) %>% 
      mutate(transcript_id = ifelse(type=="gene", NA, transcript_id))
    loaded_gff <- loaded_gff %>%   # fill protein coding only
      filter(biotype=="protein_coding") %>% GRanges() # %>% 
    # filter(seqnames != "chrM")  # remove mitochondrion
  }
  loaded_gff <- loaded_gr_seqlevelsStyle(loaded_gff)
  return(loaded_gff)
}

get_exons <- function(file_gff, mode="gencode") {
  loaded_gff <- get_gff(file_gff, mode="gencode")
  loaded_gff <- extract_exons(loaded_gff, mode="gencode")
  return(loaded_gff)
}

# helper functions for exon granges
#   extract short exons
extract_exons_by_width <- function(loaded_gff_exons, min_length=0, max_length=Inf) {
  if (min_length < 0) {
    stop("min_length < 0")
  }
  if (max_length < min_length) {
    stop("max_length < min_length")
  }
  loaded_gff_exons_by_width <- loaded_gff_exons %>% 
    filter(width >= min_length, width <= max_length)
  return(loaded_gff_exons_by_width)
}

#   extract three prime boundary
extract_exons_three_prime_boundary <- function(loaded_gff_exons, min_length=0, max_length=Inf) {
  loaded_gff_exons_by_width <- extract_exons_by_width(
    loaded_gff_exons, min_length, max_length)
  # exon next to 3' splice site, i.e. 5' end of exon
  loaded_gff_exons_three_prime_boundary <- extend_granges(
    loaded_gff_exons_by_width, 0, -(width(loaded_gff_exons_by_width)))
  return(loaded_gff_exons_three_prime_boundary)
}

#   extract five prime boundary
extract_exons_five_prime_boundary <- function(loaded_gff_exons, min_length=0, max_length=Inf) {
  loaded_gff_exons_by_width <- extract_exons_by_width(
    loaded_gff_exons, min_length, max_length)
   # exon next to 5' splice site, i.e. 3' end of exon
  loaded_gff_exons_five_prime_boundary <- extend_granges(
    loaded_gff_exons_by_width, -(width(loaded_gff_exons_by_width)), 0)
  return(loaded_gff_exons_five_prime_boundary)
}

# load and process vcf file
get_vcf <- function(file_vcf) {
  loaded_vcf <- readVcf(file_vcf)
  loaded_vcf <- rowRanges(loaded_vcf)
  loaded_vcf <- loaded_gr_seqlevelsStyle(loaded_vcf)
  # clean types and names
  loaded_vcf <- loaded_vcf %>% sort()
  loaded_vcf$REF <- as.character(loaded_vcf$REF)
  loaded_vcf$ALT <- data.frame(loaded_vcf$ALT)$value
  # add universal variant_id
  loaded_vcf$variant_name <- names(loaded_vcf)
  loaded_vcf <- loaded_vcf %>% 
    mutate(variant_id = paste(
      paste(seqnames, start, sep=":"), paste(REF, ALT, sep="/"), sep="_"))
  # label types SNP, IN, DEL
  loaded_vcf <- loaded_vcf %>% 
    mutate(variant_type = 
      ifelse((nchar(ALT) == nchar(REF)) & (nchar(REF) == 1), "SNP", 
        ifelse((nchar(ALT) == nchar(REF)) & (nchar(REF) > 1), "MNP", 
          ifelse(nchar(ALT) < nchar(REF), "DEL", 
            ifelse(nchar(ALT) > nchar(REF), "IN", NA)))))

  # filter non-ambiguous alleles
  loaded_vcf <- loaded_vcf %>% 
    as_tibble() %>% rowwise() %>% filter(str_detect(REF, "^[ATCGatcg]+$"), 
      str_detect(ALT, "^[ATCGatcg]+$")) %>% ungroup() %>% GRanges()
  return(loaded_vcf)
}

get_tbi <- function(file_tbi, loaded_exons) {
  extend_temp <- 230
  loaded_exons_temp <- loaded_exons %>% 
    extend_granges(extend_temp, extend_temp) %>% 
    mutate(strand = "*") %>% 
    reduce()
  loaded_vcf <- NULL
  header_tab <- headerTabix(gsub(
    ".tbi", "", file_tbi))$seqnames
  try({
    # default version
    loaded_exons_filter <- loaded_exons_temp %>% 
      filter(seqnames %in% header_tab)
    if (length(loaded_exons_filter)>0) {
      param <- ScanVcfParam(which=loaded_exons_filter)
      loaded_vcf <- readVcf(file_tbi, 
        genome=seqinfo(loaded_exons_filter), param=param)
    }
  }, silent=TRUE)
  if (is.null(loaded_vcf)) {
    # NCBI version
    try({
      loaded_exons_filter <- loaded_exons_temp %>% 
        loaded_gr_seqlevelsStyle("NCBI") %>% 
        filter(seqnames %in% header_tab)
      if (length(loaded_exons_filter)>0) {
        param <- ScanVcfParam(which=loaded_exons_filter)
        loaded_vcf <- readVcf(file_tbi, 
          genome=seqinfo(loaded_exons_filter), param=param)
      }
    }, silent=TRUE)
  }
  if (is.null(loaded_vcf)) {
    # UCSC version
    try({
      loaded_exons_filter <- loaded_exons_temp %>% 
        loaded_gr_seqlevelsStyle("UCSC") %>% 
        filter(seqnames %in% header_tab)
      if (length(loaded_exons_filter)>0) {
        param <- ScanVcfParam(which=loaded_exons_filter)
        loaded_vcf <- readVcf(file_tbi, 
          genome=seqinfo(loaded_exons_filter), param=param)
      }
    }, silent=TRUE)
  }
  if (is.null(loaded_vcf)) {
    # ENSEMBL version
    try({
      loaded_exons_filter <- loaded_exons_temp %>% 
        loaded_gr_seqlevelsStyle("Ensembl") %>% 
        filter(seqnames %in% header_tab)
      if (length(loaded_exons_filter)>0) {
        param <- ScanVcfParam(which=loaded_exons_filter)
        loaded_vcf <- readVcf(file_tbi, 
          genome=seqinfo(loaded_exons_filter), param=param)
      }
    }, silent=TRUE)
  }
  if (is.null(loaded_vcf)) {
    # dbSNP version
    try({
      loaded_exons_filter <- loaded_exons_temp %>% 
        loaded_gr_seqlevelsStyle("dbSNP") %>% 
        filter(seqnames %in% header_tab)
      if (length(loaded_exons_filter)>0) {
        param <- ScanVcfParam(which=loaded_exons_filter)
        loaded_vcf <- readVcf(file_tbi, 
          genome=seqinfo(loaded_exons_filter), param=param)
      }
    }, silent=TRUE)
  }
  # finalize
  if (!is.null(loaded_vcf)) {
    loaded_vcf <- rowRanges(loaded_vcf)
    loaded_vcf <- loaded_vcf %>% 
      loaded_gr_seqlevelsStyle() 
    # clean types and names
    loaded_vcf <- loaded_vcf %>% sort()
    loaded_vcf$REF <- as.character(loaded_vcf$REF)
    loaded_vcf$ALT <- data.frame(loaded_vcf$ALT)$value
    # add universal variant_id
    loaded_vcf$variant_name <- names(loaded_vcf)
    loaded_vcf <- loaded_vcf %>% 
      mutate(variant_id = paste(
        paste(seqnames, start, sep=":"), paste(REF, ALT, sep="/"), sep="_"))
    # label type SNP, IN, DEL
    loaded_vcf <- loaded_vcf %>% 
      mutate(variant_type = 
        ifelse((nchar(ALT) == nchar(REF)) & (nchar(REF) == 1), "SNP", 
          ifelse((nchar(ALT) == nchar(REF)) & (nchar(REF) > 1), "MNP", 
            ifelse(nchar(ALT) < nchar(REF), "DEL", 
              ifelse(nchar(ALT) > nchar(REF), "IN", NA)))))
  }

  # filter non-ambiguous alleles
  loaded_vcf <- loaded_vcf %>% 
    as_tibble() %>% rowwise() %>% filter(str_detect(REF, "^[ATCGatcg]+$"), 
      str_detect(ALT, "^[ATCGatcg]+$")) %>% ungroup() %>% GRanges()
  return(loaded_vcf)
}

# create the sequence corresponding to the variable variable
create_variable_sequence <- function(loaded_fasta, loaded_vcf, loaded_gff_exons, 
  loaded_gff_exons_var, loaded_gff_exons_buf, description=NA, design=NA, switch_case=F) {  
  # get overlaps
  overlaps <- findOverlaps(loaded_vcf, loaded_gff_exons_var, type="within")

  # get seqs_var
  seqs_var <- getSeq(loaded_fasta, loaded_gff_exons_var[subjectHits(overlaps)])
  seqs_buf <- getSeq(loaded_fasta, loaded_gff_exons_buf[subjectHits(overlaps)])

  # convert to variable variable table
  loaded_gff_exons_temp <- loaded_gff_exons %>% as_tibble() %>% 
    dplyr::select(seqnames, start, end, width, strand, exon_id, exon_name, exon_rank, exon_internal, 
      gene_id, gene_type, gene_name, transcript_id, transcript_type, transcript_name)
  exon_table <- loaded_gff_exons_temp[subjectHits(overlaps),]
  colnames(exon_table) <- paste("exon_", colnames(exon_table), sep="")
  
  # convert to variable variable table
  loaded_gff_exons_var_temp <- loaded_gff_exons_var %>% as_tibble() %>% 
    dplyr::select(seqnames, start, end, width, strand, exon_id, exon_name, exon_rank, exon_internal, 
      gene_id, gene_type, gene_name, transcript_id, transcript_type, transcript_name)
  variable_table <- loaded_gff_exons_var_temp[subjectHits(overlaps),]
  colnames(variable_table) <- paste("variable_", colnames(variable_table), sep="")

  # convert to buffer variable table
  loaded_gff_exons_buf_temp <- loaded_gff_exons_buf %>% as_tibble() %>% 
    dplyr::select(seqnames, start, end, width, strand, exon_id, exon_name, exon_rank, exon_internal, 
      gene_id, gene_type, gene_name, transcript_id, transcript_type, transcript_name)
  buffer_table <- loaded_gff_exons_buf_temp[subjectHits(overlaps),]
  colnames(buffer_table) <- paste("buffer_", colnames(buffer_table), sep="")
  
  # convert to overlap variant table
  variant_temp <- loaded_vcf[queryHits(overlaps)]
  variant_table <- as_tibble(variant_temp) %>% 
    dplyr::select(-paramRangeID)
  names(variant_table) <- paste("variant_", names(variant_table), sep="")
  names(variant_table) <- gsub("variant_variant", "variant", names(variant_table))

  # drop extraneous columns
  variable_table <- variable_table %>% 
    dplyr::select(-c(
    variable_exon_id, variable_exon_name, variable_exon_rank, 
      variable_exon_internal, variable_gene_id, variable_gene_type, variable_gene_name, 
      variable_transcript_id, variable_transcript_type, variable_transcript_name))
  buffer_table <- buffer_table %>% 
    dplyr::select(-c(
      buffer_exon_id, buffer_exon_name, buffer_exon_rank, 
      buffer_exon_internal, buffer_gene_id, buffer_gene_type, buffer_gene_name, 
      buffer_transcript_id, buffer_transcript_type, buffer_transcript_name))

  # join the above three tables
  variant_joined_table <- bind_cols(bind_cols(bind_cols(
    variant_table, exon_table), variable_table), buffer_table)
  
  # add metadata column
  variant_joined_table <- variant_joined_table %>% 
    mutate(description = description)
  variant_joined_table <- variant_joined_table %>% 
    mutate(design = design)

  # check if variant not allowed
  variant_joined_table <- variant_joined_table %>% 
    mutate(relative_variable_start = 
      ifelse(variable_strand=="+",
        variant_start - variable_start,
        ifelse(variable_strand=="-",
          variable_end - variant_end,
          NA))) %>% 
    mutate(relative_variable_end = 
      ifelse(variable_strand=="+",
        variant_end - variable_start,
        ifelse(variable_strand=="-",
          variable_end - variant_start,
          NA))) %>% 
    mutate(relative_buffer_start = 
      ifelse(buffer_strand=="+",
        variant_start - buffer_start,
        ifelse(buffer_strand=="-",
          buffer_end - variant_end,
          NA))) %>% 
    mutate(relative_buffer_end = 
      ifelse(buffer_strand=="+",
        variant_end - buffer_start,
        ifelse(buffer_strand=="-",
          buffer_end - variant_start,
          NA)))

  variant_joined_table <- variant_joined_table %>% 
    mutate(flag_var_in_variable = (
      variant_start < variable_start) | (variant_end > variable_end))

  # create wildtype allele
  variant_joined_table <- variant_joined_table %>% 
    mutate(variable_wt_allele = 
      ifelse(
        variable_strand=="+", variant_REF,
      ifelse(
        variable_strand=="-", stri_reverse(chartr("ATGC", "TACG", variant_REF)),
        NA)))
  
  # create mutant allele
  variant_joined_table <- variant_joined_table %>% 
    mutate(variable_mt_allele = 
      ifelse(
        variable_strand=="+", variant_ALT,
      ifelse(
        variable_strand=="-", stri_reverse(chartr("ATGC", "TACG", variant_ALT)),
        NA)))

  # create allele
  variant_joined_table <- variant_joined_table %>% 
    mutate(variable_alleles = paste(variable_wt_allele, variable_mt_allele, sep="/"))

  # create sequence
  variant_joined_table$variable_wt_seq <- 
    enframe(as.vector(seqs_var))$value
  
  # create wildtype sequence
  variant_joined_table <- variant_joined_table %>% 
    mutate(variable_wt_seq = ifelse(
      variable_strand=="+",
      paste(
        substr(
          variable_wt_seq, 1, variant_start - variable_start), chartr("ATGCatgc", "atgcATGC", variant_REF), 
        substr(
          variable_wt_seq, variant_width + variant_start - variable_start + 1, variable_width),
        sep=""),
      ifelse(
        variable_strand=="-",
        paste(
          substr(
            variable_wt_seq, 1, variable_end - variant_end), chartr(
              "ATGCatgc", "atgcATGC", stri_reverse(chartr("ATGC", "TACG", variant_REF))), 
          substr(
            variable_wt_seq, variant_width + variable_end - variant_end + 1, variable_width),
          sep=""),
        NA)))
  
  # create mutant sequence
  variant_joined_table <- variant_joined_table %>% 
    mutate(variable_mt_seq = ifelse(
      variable_strand=="+",
      paste(
        substr(
          variable_wt_seq, 1, variant_start - variable_start), chartr("ATGCatgc", "atgcATGC", variant_ALT), 
        substr(
          variable_wt_seq, variant_width + variant_start - variable_start + 1, variable_width),
        sep=""),
      ifelse(
        variable_strand=="-",
        paste(
          substr(
            variable_wt_seq, 1, variable_end - variant_end), chartr(
              "ATGCatgc", "atgcATGC", stri_reverse(chartr("ATGCatgc", "TACGtacg", variant_ALT))), 
          substr(
            variable_wt_seq, variant_width + variable_end - variant_end + 1, variable_width),
          sep=""),
        NA)))
  
  # check wildtype sequence
  variant_joined_table <- variant_joined_table %>%
    mutate(
      flag_ref_allele = (
        variant_REF !=
          ifelse(variable_strand == "+",
                 toupper(substr(variable_wt_seq, variant_start - variable_start + 1,
                                variant_start - variable_start + nchar(variant_REF))),
                 ifelse(variable_strand == "-",
                        stri_reverse(chartr("ATGCatgc", "TACGtacg", toupper(substr(
                          variable_wt_seq, variable_end - variant_end + 1,
                          variable_end - variant_end + nchar(variant_REF))))), NA))))
  # check ambiguous alleles
  # sum(variant_joined_table$flag_ref_allele, na.rm=T)
  
  # check mutation sequence
  variant_joined_table <- variant_joined_table %>%
    mutate(
      flag_alt_allele = (
        variant_ALT !=
          ifelse(variable_strand == "+",
                 toupper(substr(variable_mt_seq, variant_start - variable_start + 1,
                                variant_start - variable_start + nchar(variant_ALT))),
                 ifelse(variable_strand == "-",
                        stri_reverse(chartr("ATGCatgc", "TACGtacg", toupper(substr(
                          variable_mt_seq, variable_end - variant_end + 1,
                          variable_end - variant_end + nchar(variant_ALT))))), NA))))
  # check ambiguous alleles
  # sum(variant_joined_table$flag_alt_allele, na.rm=T)
  
  # mutant vs wildtype length
  variant_joined_table$variable_wt_width <- nchar(variant_joined_table$variable_wt_seq)
  variant_joined_table$variable_mt_width <- nchar(variant_joined_table$variable_mt_seq)
  variant_joined_table$flag_mt_wt_delta_width <- 
    (variant_joined_table$variable_mt_width - variant_joined_table$variable_wt_width) != 0
  variant_joined_table$flag_mt_wt_string_dist <- (stringdist(
    variant_joined_table$variable_mt_seq, variant_joined_table$variable_wt_seq)) != 1

  # add sequence for buffer variables
  variant_joined_table$buffer_seq_left <- 
    enframe(as.vector(seqs_buf))$value
  variant_joined_table <- variant_joined_table %>% 
    mutate(buffer_seq_left = 
      ifelse(variable_strand=="+",
        substr(buffer_seq_left, 1, variable_start-buffer_start),
      ifelse(variable_strand=="-",
        substr(buffer_seq_left, 1, buffer_end-variable_end),
      NA))) %>% 
    mutate(buffer_seq_left = tolower(buffer_seq_left))

  variant_joined_table$buffer_seq_right <- 
    enframe(as.vector(seqs_buf))$value
  variant_joined_table <- variant_joined_table %>% 
    mutate(buffer_seq_right = 
      ifelse(variable_strand=="+",
        substr(buffer_seq_right, (variable_end+1)-(buffer_start-1), nchar(buffer_seq_right)),
      ifelse(variable_strand=="-",
        substr(buffer_seq_right, (buffer_end+1)-(variable_start-1), nchar(buffer_seq_right)),
      NA))) %>% 
    mutate(buffer_seq_right = tolower(buffer_seq_right))

  variant_joined_table$buffer_width_left <- nchar(variant_joined_table$buffer_seq_left)
  variant_joined_table$buffer_width_right <- nchar(variant_joined_table$buffer_seq_right)

  # check concatenate is buffer_seq
  variant_joined_table$flag_buffer_seq <- !toupper(
    enframe(as.vector(seqs_buf))$value) == toupper(paste(
      variant_joined_table$buffer_seq_left, 
      variant_joined_table$variable_wt_seq, 
      variant_joined_table$buffer_seq_right, sep=""))

  # create final wt and mt sequence
  variant_joined_table$final_wt_seq <- 
    paste(variant_joined_table$buffer_seq_left, 
      variant_joined_table$variable_wt_seq, 
      variant_joined_table$buffer_seq_right, sep="")

  variant_joined_table$final_mt_seq <- 
    paste(variant_joined_table$buffer_seq_left, 
      variant_joined_table$variable_mt_seq, 
      variant_joined_table$buffer_seq_right, sep="")

  variant_joined_table$final_wt_width <- nchar(variant_joined_table$final_wt_seq)
  variant_joined_table$final_mt_width <- nchar(variant_joined_table$final_mt_seq)
  variant_joined_table$flag_final_width <- 
    (variant_joined_table$final_wt_width != variant_joined_table$final_mt_width)

  # switch case of wt and mt sequences
  if (switch_case) {
    variant_joined_table <- variant_joined_table %>%
      mutate(final_wt_seq = chartr("ATGCatgc", "atgcATGC", final_wt_seq)) %>%
      mutate(final_mt_seq = chartr("ATGCatgc", "atgcATGC", final_mt_seq))
  }

  # filter non-ambiguous sequences
  variant_joined_table <- variant_joined_table %>% 
    rowwise() %>% filter(str_detect(final_wt_seq, "^[ATCGatcg]+$"), 
      str_detect(final_mt_seq, "^[ATCGatcg]+$")) %>% ungroup()    
  # output final table
  return(variant_joined_table)
}

# create the sequence corresponding to nonvariable variables
create_nonvariable_sequence <- function(loaded_fasta, loaded_gff_exons, metadata_name, switch_case=F) {
  seqs_var <- getSeq(loaded_fasta, loaded_gff_exons)
  seqs_var <- enframe(as.vector(seqs_var))
  if (switch_case) {
    seqs_var <- seqs_var %>% 
      mutate(value = chartr("ATGCatgc", "atgcATGC", value))
  }
  names(seqs_var) <- metadata_name
  return(seqs_var)
}

# check any columns with flag in column name
check_flags_throw_warnings <- function(library_sequence_table) {
  for (column in which(grepl("flag", names(eval(parse(text=library_sequence_table)))))) {
    if (sum(is.numeric(eval(parse(text=library_sequence_table))[,column])) != 0) {
      warning(paste(library_sequence_table, "has problem with", 
        names(eval(parse(text=library_sequence_table))))[column])
    }
  }
}
