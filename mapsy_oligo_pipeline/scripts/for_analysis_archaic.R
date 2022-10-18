#!/usr/bin/env Rscript

# Create input and output reference genomes for alignment and variant to sequence association tables.

# load files
library(tidyverse)
library(data.table)

# switch case check
switchcase <- function(sequence) {
  return(chartr("ATGCatgc", "atgcATGC", sequence))
}

# join neanderthal tables
neanderthal_updated_table <- as_tibble(fread("../results/mapsy_orders/library_2/neanderthal_updated_library_exon_final.txt.gz"))
# neanderthal_updated_table <- neanderthal_updated_table %>% filter((construct_type != "short_full_exon_variant_five"))  # remove this

standardize_table <- function(table) {
  table %>% 
  mutate(Common_id_three = gsub("Common_", "Common_three_", Common_id_three)) %>% 
  mutate(Common_id_three = gsub("Common_three_three", "Common_three", Common_id_three)) %>% 
  mutate(Common_id_five = gsub("Common_", "Common_five_", Common_id_five)) %>% 
  mutate(Common_id_five = gsub("Common_five_five", "Common_five", Common_id_five)) %>% 
  .[,c(which(colnames(.)=="source"),which(colnames(.)!="source"))]
}

neanderthal_updated_table <- neanderthal_updated_table %>% 
  mutate(source = "Neanderthal_updated") %>% 
  standardize_table()

write_tsv(neanderthal_updated_table, gzfile("../results/for_analysis/neanderthal_updated_table.txt.gz"))

# fasta reference inputs
#   add ad81 exon intron, variable exon intron parts, extension, actin exon intron
# fasta reference outputs
#   add ad81 exon, variable exon, actin exon

# add minigene constructs
first_exon <- "GTTCGTCCTCACTCTCTTCCGCATCGCTGTCTGCGAGGGCCAGCTGTTGGG"
first_intron <- "GTGAGTACTCCCTCTCAAAAGCGGGCATGACTTCTGCCCTCGAGTTATTAACCCTCACTAAAGGCAGTAGTCAAGGGTTTCCTTGAAGCTTTCGTGCTCTTATGCGAACGTGA"
last_intron <- "AAAGGATTGCTGCTCTAAGTTTCATAAAGGCAGGGATTTTTGTCCATTTTAATTGATTTTTGTCAATTTCATTTAAAATTGTCATATAGTAAACACTCATTGTTTGTCCAAATGTCAAATGACTGGGTTTCCAGAACAGGTTGACACCTTTACTCCTGTTCTGTTACTGACCGGCCCTTCCCCTTCTGGGTCCTCCATTTCCTCATCTGGCACATGGACAGATGGCTGACCCATCCTGAAAGTCGCCTTGTTGTTCCTGCCCCAG"
last_exon <- "GTGCGGCAGCTGGTGCCTCGGAGGGACCAAGCTCTGACGGAGGAGCATGCCCGACAGCAGCACAATGAGAGGCTACGCAAGCAGTTTGGAGC"
variable_extension <- "TCGATTGTGCATAAC"

temp_list <- c("short_full_exon_variant_five", "short_half_exon_variant_five", 
  "short_half_exon_variant_five_legacy", "long_half_exon_variant_five")
neanderthal_updated_table_ref <- neanderthal_updated_table %>% 
  mutate(ref_var_intron_up_wt = ifelse(!(design %in% temp_list), 
    substr(Order_final_sequence_wt, 1, 25+55+(120-pmin(120, exon_width))),
    substr(Order_final_sequence_wt, 1, 25+55))) %>% 
  mutate(ref_var_exon_middle_wt = ifelse(!(design %in% temp_list), 
    substr(Order_final_sequence_wt, 25+55+(120-pmin(120, exon_width))+1, 25+55+120),
    substr(Order_final_sequence_wt, 25+55+1, 25+55+pmin(120, exon_width)))) %>% 
  mutate(ref_var_intron_down_wt = ifelse(!(design %in% temp_list), 
    substr(Order_final_sequence_wt, 25+55+120+1, nchar(Order_final_sequence_wt)), 
    substr(Order_final_sequence_wt, 25+55+pmin(120, exon_width)+1, nchar(Order_final_sequence_wt)))) %>% 
  mutate(ref_var_intron_up_mt = ifelse(!(design %in% temp_list), 
    substr(Order_final_sequence_mt, 1, 25+55+(120-pmin(120, exon_width))),
    substr(Order_final_sequence_mt, 1, 25+55))) %>% 
  mutate(ref_var_exon_middle_mt = ifelse(!(design %in% temp_list), 
    substr(Order_final_sequence_mt, 25+55+(120-pmin(120, exon_width))+1, 25+55+120),
    substr(Order_final_sequence_mt, 25+55+1, 25+55+pmin(120, exon_width)))) %>% 
  mutate(ref_var_intron_down_mt = ifelse(!(design %in% temp_list), 
    substr(Order_final_sequence_mt, 25+55+120+1, nchar(Order_final_sequence_mt)), 
    substr(Order_final_sequence_mt, 25+55+pmin(120, exon_width)+1, nchar(Order_final_sequence_mt))))

write_tsv(neanderthal_updated_table_ref %>% 
  dplyr::select(design, ref_var_intron_up_wt, ref_var_exon_middle_wt, ref_var_intron_down_wt, 
    ref_var_intron_up_mt, ref_var_exon_middle_mt, ref_var_intron_down_mt) %>% 
  group_by(design) %>% slice(c(1,1)), "../results/for_analysis/neanderthal_updated_table_ref_check.txt")

neanderthal_updated_table_ref <- neanderthal_updated_table_ref %>% 
  mutate(ref_input_debug_wt = paste(
    toupper(first_exon), tolower(first_intron), 
    ref_var_intron_up_wt, ref_var_exon_middle_wt, ref_var_intron_down_wt, 
    tolower(variable_extension), tolower(last_intron), toupper(last_exon),
    sep="_")) %>% 
  mutate(ref_input_debug_mt = paste(
    toupper(first_exon), tolower(first_intron), 
    ref_var_intron_up_mt, ref_var_exon_middle_mt, ref_var_intron_down_mt, 
    tolower(variable_extension), tolower(last_intron), toupper(last_exon),
    sep="_")) %>% 
  mutate(ref_output_debug_wt = paste(
    toupper(first_exon), switchcase(ref_var_exon_middle_wt), toupper(last_exon),
    sep="_")) %>% 
  mutate(ref_output_debug_mt = paste(
    toupper(first_exon), switchcase(ref_var_exon_middle_mt), toupper(last_exon),
    sep="_"))

neanderthal_updated_table_ref <- neanderthal_updated_table_ref %>% 
  mutate(ref_input_final_wt = gsub("_", "", ref_input_debug_wt)) %>% 
  mutate(ref_input_final_mt = gsub("_", "", ref_input_debug_mt)) %>% 
  mutate(ref_output_final_wt = gsub("_", "", ref_output_debug_wt)) %>% 
  mutate(ref_output_final_mt = gsub("_", "", ref_output_debug_mt))

neanderthal_updated_table_ref <- neanderthal_updated_table_ref %>% 
  mutate(ref_input_seq_wt = toupper(ref_input_final_wt)) %>% 
  mutate(ref_input_seq_mt = toupper(ref_input_final_mt)) %>% 
  mutate(ref_output_seq_wt = toupper(ref_output_final_wt)) %>% 
  mutate(ref_output_seq_mt = toupper(ref_output_final_mt))

# save fasta without dups
neanderthal_updated_input_seq <- unique(c(
  toupper(neanderthal_updated_table_ref$ref_input_final_wt), 
  toupper(neanderthal_updated_table_ref$ref_input_final_mt)))
neanderthal_updated_input_id <- paste("neanderthal_updated_input_", 
  as.character(c(1:length(neanderthal_updated_input_seq))), sep="")
neanderthal_updated_input_fasta <- c(rbind(
  paste(">", neanderthal_updated_input_id, sep=""), 
  neanderthal_updated_input_seq))
write(neanderthal_updated_input_fasta, gzfile("../results/for_analysis/neanderthal_updated_input.fasta.gz"))

neanderthal_updated_output_seq <- unique(c(
  toupper(neanderthal_updated_table_ref$ref_output_final_wt), 
  toupper(neanderthal_updated_table_ref$ref_output_final_mt)))
neanderthal_updated_output_id <- paste("neanderthal_updated_output_", 
  as.character(c(1:length(neanderthal_updated_output_seq))), sep="")
neanderthal_updated_output_fasta <- c(rbind(
  paste(">", neanderthal_updated_output_id, sep=""), 
  neanderthal_updated_output_seq))
neanderthal_updated_output_fasta <- c(
  neanderthal_updated_output_fasta, ">exon_skipped", 
  paste(first_exon, last_exon, sep=""))
write(neanderthal_updated_output_fasta, gzfile("../results/for_analysis/neanderthal_updated_output.fasta.gz"))

# finalize table
temp_input <- as_tibble(cbind(neanderthal_updated_input_id, neanderthal_updated_input_seq))
names(temp_input) <- c("id", "seq")

temp_output <- as_tibble(cbind(neanderthal_updated_output_id, neanderthal_updated_output_seq))
names(temp_output) <- c("id", "seq")

neanderthal_updated_table_ref <- neanderthal_updated_table_ref %>% 
  left_join(temp_input %>% dplyr::rename(ref_input_id_wt=id, ref_input_seq_wt=seq)) %>% 
  left_join(temp_input %>% dplyr::rename(ref_input_id_mt=id, ref_input_seq_mt=seq)) %>% 
  left_join(temp_output %>% dplyr::rename(ref_output_id_wt=id, ref_output_seq_wt=seq)) %>% 
  left_join(temp_output %>% dplyr::rename(ref_output_id_mt=id, ref_output_seq_mt=seq))

write_tsv(neanderthal_updated_table_ref, gzfile("../results/for_analysis/neanderthal_updated_table_ref.txt.gz"))
