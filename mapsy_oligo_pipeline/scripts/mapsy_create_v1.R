#!/usr/bin/env Rscript

# Command line script which takes as input a fasta reference genome, a gff transcript annotation, a vcf variant file, a param file, and a list of gene names to filter, and generates oligos.

# load packages
print("loading preamble")
source("mapsy_helper_v1.R")
library(optparse)

# inputs
option_list = list(
  make_option(
    c("-p", "--param"), type="character", default=NULL, help="param file location", metavar="character"), 
  make_option(
    c("-f", "--fasta"), type="character", default=NULL, help="fasta file location", metavar="character"), 
  make_option(
    c("-g", "--gff"), type="character", default=NULL, help="gff3 file location", metavar="character"), 
  make_option(
    c("-v", "--vcf"), type="character", default=NULL, help="vcf file location", metavar="character"), 
  make_option(
    c("-e", "--genes"), type="character", default=NULL, help="gene list file location", metavar="character"), 
  make_option(
    c("-o", "--out_prefix"), type="character", default=NULL, help="output location prefix", metavar="character"),
  make_option(
    c("-d", "--desc_prefix"), type="character", default=NULL, help="description prefix", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

file_param <- opt$param
file_fasta <- opt$fasta
file_gff <- opt$gff
file_vcf <- opt$vcf
file_genes_list <- opt$genes
file_out_prefix <- opt$out_prefix
file_desc_prefix <- opt$desc_prefix

# load params
source(file_param)
if (!(oligo_length == intron_extend_three + exon_length_max + intron_extend_five)) {
  stop ('oligo_length != intron_extend_three + exon_legth_max + intron_extend_five')
}
if (!(exon_length_max == exon_variant_half + exon_variant_buffer + exon_variant_common)) {
  stop ('exon_length_max 1= exon_variant_half + exon_variant_buffer + exon_variant_common')
}

# load fasta
print("loading fasta file")
loaded_fasta <- get_fasta(file_fasta)

# process gff files
# to save time in future
# store processed versions
print("loading gff file")
file_gff_temp <- gsub(".gz", ".rds", file_gff)
if (!file.exists(file_gff_temp)) {
  loaded_exons <- get_exons(file_gff)
  saveRDS(loaded_exons, file_gff_temp)
} else {
  loaded_exons <- readRDS(file_gff_temp)
}

# filter by gene list if available
if (!is.null(file_genes_list)) {
  # load genes list
  loaded_genes <- as.character(read.table(file_genes_list)[,1])
  # check uniqueness of gene list
  log_loaded_genes <- enframe(loaded_genes) %>% 
    dplyr::select(-name) %>% 
    dplyr::rename(gene = value) %>% 
    arrange(gene) %>% 
    add_count(gene, name="gene_is_unique") %>% 
    mutate(gene_is_unique = (gene_is_unique == 1)) %>% 
    mutate(gene_is_in_gff = gene %in% loaded_exons$gene_name)
  # save uniqueness of gene list
  write_tsv(log_loaded_genes, 
    paste(file_out_prefix, "log_loaded_genes.txt", sep="_"))
  # now filter by gene list
  loaded_exons <- filter(loaded_exons, 
    gene_name %in% loaded_genes)
}  # else, use all exons in genome

# create a bed file that can be used to tabix vcf file
#   since many variants may be outside of testable regions
extend_temp <- 230
file_gff_bed <- gsub(".gz", ".bed", file_gff)
loaded_exons_temp <- extend_granges(
  loaded_exons, extend_temp, extend_temp)
seqlevelsStyle(loaded_exons_temp) <- "NCBI"
strand(loaded_exons_temp) <- "*"
loaded_exons_temp <- sort(trim(reduce(
  loaded_exons_temp)))
rtracklayer::export(loaded_exons_temp, 
  file_gff_bed, format="BED")

# add metadata on internal exons
#   on per transcript basis
loaded_exons <- loaded_exons %>% 
  as_tibble() %>%  # faster than plyranges
  group_by(transcript_name) %>% 
  mutate(exon_internal = !(exon_rank %in% 
    c(min(exon_rank), max(exon_rank)))) %>% 
  ungroup() %>% GRanges() %>% unique()

# load vcf file using exon ranges
if (grepl(".tbi$|.idx$", file_vcf)) {
  print("loading tbi file")
  loaded_vcf <- get_tbi(
    file_vcf, loaded_exons)
} else {
  print("loading vcf file")
  loaded_vcf <- get_vcf(file_vcf)
}

if (is.null(loaded_vcf)) {
  stop("no variants found in given ranges")
}

# filter SNPs only
loaded_vcf <- loaded_vcf %>% 
  filter(variant_type %in% c("SNP", "MNP"))

# design exon variant constructs
print("creating sequences")

# CREATE EXONIC OLIGO SEQUENCES
# 
# # 
# # # 
#   short variable variable
#     old short version
short_full_exon_variant_three <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max), 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max), 
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_full, exon_length_max), 
    oligo_length-intron_extend_five, intron_extend_five), 
  paste(file_desc_prefix, "exonic", sep="_"), "short_full_exon_variant_three")

#      new short version
short_full_exon_variant_five <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max), 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max), 
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_full, exon_length_max), 
    intron_extend_three, oligo_length-intron_extend_three), 
  paste(file_desc_prefix, "exonic", sep="_"), "short_full_exon_variant_five")

#   short half exon variable
temp_width_exon_min_half_exon_max <- width(extract_exons_by_width(
  loaded_exons, exon_length_min_half, exon_length_max))
#     old short three version
short_half_exon_variant_three <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_half, exon_length_max),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    0, temp_width_exon_min_half_exon_max - exon_variant_common - exon_variant_buffer),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    oligo_length-intron_extend_five, -exon_variant_common),
  paste(file_desc_prefix, "exonic", sep="_"), "short_half_exon_variant_three")
#     new short five version
short_half_exon_variant_five <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_half, exon_length_max),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    temp_width_exon_min_half_exon_max - exon_variant_common - exon_variant_buffer, 0),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    -exon_variant_common, oligo_length-intron_extend_three),
  paste(file_desc_prefix, "exonic", sep="_"), "short_half_exon_variant_five")

#   short half exon variable alt
#     old short five version
short_half_exon_variant_five_legacy <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_half, exon_length_max),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    temp_width_exon_min_half_exon_max - exon_variant_common - exon_variant_buffer, 0),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    -exon_variant_common, temp_width_exon_min_half_exon_max+intron_extend_five),
  paste(file_desc_prefix, "exonic", sep="_"), "short_half_exon_variant_five_legacy")

#   long three half exon variable
long_half_exon_variant_three <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_max+1, Inf),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    0, exon_variant_half),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    intron_extend_three, exon_variant_half+exon_variant_buffer),
  paste(file_desc_prefix, "exonic", sep="_"), "long_half_exon_variant_three")

#   long five half exon variable
long_half_exon_variant_five <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_max+1, Inf),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    exon_variant_half, 0),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    exon_variant_half+exon_variant_buffer, intron_extend_five),
  paste(file_desc_prefix, "exonic", sep="_"), "long_half_exon_variant_five")
# # # 
# # 
# 

# check flags, throw warnings
check_flags_throw_warnings("short_full_exon_variant_three")
check_flags_throw_warnings("short_full_exon_variant_five")
check_flags_throw_warnings("short_half_exon_variant_three")
check_flags_throw_warnings("short_half_exon_variant_five")
check_flags_throw_warnings("short_half_exon_variant_five_legacy")
check_flags_throw_warnings("long_half_exon_variant_three")
check_flags_throw_warnings("long_half_exon_variant_five")


# save output
write_tsv(short_full_exon_variant_three, 
  gzfile(paste(file_out_prefix, "short_full_exon_variant_three.txt.gz", sep="_")))
write_tsv(short_full_exon_variant_five, 
  gzfile(paste(file_out_prefix, "short_full_exon_variant_five.txt.gz", sep="_")))
write_tsv(short_half_exon_variant_three, 
  gzfile(paste(file_out_prefix, "short_half_exon_variant_three.txt.gz", sep="_")))
write_tsv(short_half_exon_variant_five, 
  gzfile(paste(file_out_prefix, "short_half_exon_variant_five.txt.gz", sep="_")))
write_tsv(short_half_exon_variant_five_legacy, 
  gzfile(paste(file_out_prefix, "short_half_exon_variant_five_legacy.txt.gz", sep="_")))
write_tsv(long_half_exon_variant_three, 
  gzfile(paste(file_out_prefix, "long_half_exon_variant_three.txt.gz", sep="_")))
write_tsv(long_half_exon_variant_five, 
  gzfile(paste(file_out_prefix, "long_half_exon_variant_five.txt.gz", sep="_")))

# CREATE INTRONIC OLIGO SEQUENCES
# 
# #
# # # 
# short variable variable
#     old short version
temp_width_intron_min_full_intron_max <- width(extract_exons_by_width(
  loaded_exons, exon_length_min_full, exon_length_max))

short_full_intron_variant_three <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max),
  c(extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_full, exon_length_max), 
    oligo_length-intron_extend_five-intron_variant_buffer, 
    0-temp_width_intron_min_full_intron_max)), 
  c(extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_full, exon_length_max), 
    oligo_length-intron_extend_five, intron_extend_five)), 
  paste(file_desc_prefix, "intronic", sep="_"), "short_full_intron_variant_three", switch_case=T)
# check 3' splice site sequence logo
# ggseqlogo(short_full_intron_variant_three$final_mt_seq)

#     new short version
short_full_intron_variant_five <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max),
  c(extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_full, exon_length_max), 
    0-temp_width_intron_min_full_intron_max, 
    oligo_length-intron_extend_three-intron_variant_buffer)), 
  c(extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_full, exon_length_max), 
    intron_extend_three, oligo_length-intron_extend_three)), 
  paste(file_desc_prefix, "intronic", sep="_"), "short_full_intron_variant_five", switch_case=T)
# check 5' splice site sequence logo
# ggseqlogo(short_full_intron_variant_five$final_mt_seq)

#   short half intron variable
temp_width_intron_min_half_intron_max <- width(extract_exons_by_width(
  loaded_exons, exon_length_min_half, exon_length_max))
short_half_intron_variant_three <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    oligo_length-intron_extend_five-intron_variant_buffer, 
    0-temp_width_intron_min_half_intron_max),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    oligo_length-intron_extend_five, -exon_variant_common),
  paste(file_desc_prefix, "intronic", sep="_"), "short_half_intron_variant_three", switch_case=T)
short_half_intron_variant_five <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_min_full, exon_length_max),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    0-temp_width_intron_min_half_intron_max, 
    oligo_length-intron_extend_three-intron_variant_buffer),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_min_half, exon_length_max), 
    -exon_variant_common, oligo_length-intron_extend_three),
  paste(file_desc_prefix, "intronic", sep="_"), "short_half_intron_variant_five", switch_case=T)

#   long three half intron variable
long_half_intron_variant_three <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_max+1, Inf),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    intron_extend_three-intron_variant_buffer, 0),
  extend_granges(extract_exons_three_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    intron_extend_three, exon_variant_half+exon_variant_buffer),
  paste(file_desc_prefix, "intronic", sep="_"), "long_half_intron_variant_three", switch_case=T)

#   long five half intron variable
long_half_intron_variant_five <- create_variable_sequence(loaded_fasta, loaded_vcf, 
  extract_exons_by_width(loaded_exons, exon_length_max+1, Inf),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    0, intron_extend_five-intron_variant_buffer),
  extend_granges(extract_exons_five_prime_boundary(
    loaded_exons, exon_length_max+1, Inf), 
    exon_variant_half+exon_variant_buffer, intron_extend_five),
  paste(file_desc_prefix, "intronic", sep="_"), "long_half_intron_variant_five", switch_case=T)
# # # 
# # 
# 

# check flags, throw warnings
check_flags_throw_warnings("short_full_intron_variant_three")
check_flags_throw_warnings("short_full_intron_variant_five")
check_flags_throw_warnings("short_half_intron_variant_three")
check_flags_throw_warnings("short_half_intron_variant_five")
check_flags_throw_warnings("long_half_intron_variant_three")
check_flags_throw_warnings("long_half_intron_variant_five")

# save output
write_tsv(short_full_intron_variant_three, 
  gzfile(paste(file_out_prefix, "short_full_intron_variant_three.txt.gz", sep="_")))
write_tsv(short_full_intron_variant_five, 
  gzfile(paste(file_out_prefix, "short_full_intron_variant_five.txt.gz", sep="_")))
write_tsv(short_half_intron_variant_three, 
  gzfile(paste(file_out_prefix, "short_half_intron_variant_three.txt.gz", sep="_")))
write_tsv(short_half_intron_variant_five, 
  gzfile(paste(file_out_prefix, "short_half_intron_variant_five.txt.gz", sep="_")))
write_tsv(long_half_intron_variant_three, 
  gzfile(paste(file_out_prefix, "long_half_intron_variant_three.txt.gz", sep="_")))
write_tsv(long_half_intron_variant_five, 
  gzfile(paste(file_out_prefix, "long_half_intron_variant_five.txt.gz", sep="_")))
