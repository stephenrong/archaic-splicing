#!/bin/R

# Computes mpralm scores from MaPSy experiment.

library(mpra)
library(DescTools)

analyze_experiment <- function(experiment, mpralm_mode=TRUE, chisq_mode=TRUE) {
	# mpralm analysis
	print("preamble")
	experiment_dna <- as_tibble(experiment$input_counts_allelic + 1)  # pseudocount
	names(experiment_dna) <- as.vector(sapply(c(1:(ncol(experiment_dna)/2)), function(x) {
		c(paste("input_allele1_sample", x, sep=""), paste("input_allele2_sample", x, sep=""))}))
	experiment_rna <- as_tibble(experiment$output_spliced_counts_allelic + 1)  # pseudocount
	names(experiment_rna) <- as.vector(sapply(c(1:(ncol(experiment_rna)/2)), function(x) {
		c(paste("output_allele1_sample", x, sep=""), paste("output_allele2_sample", x, sep=""))}))
	experiment_eid <- as.vector(paste("variant_pair", c(1:nrow(experiment_dna)), sep="_"))
	
	experiment_mpralm_final <- as_tibble(cbind(
		experiment_eid, experiment_dna, experiment_rna)) %>% 
		mutate(experiment_eid = as.character(experiment_eid))

	# mpralm method
	print("mpralm method")
	if (mpralm_mode) {
		experiment_block <- as.vector(
			sapply(c(1:(ncol(experiment_dna)/2)), function(x) {c(x, x)}))
		experiment_mpraset <- MPRASet(
			DNA = experiment_dna, RNA = experiment_rna,
			eid = experiment_eid, eseq = NULL, barcode = NULL)
		experiment_design <- data.frame(
			intcpt = 1, allele1 = grepl("allele1", 
				colnames(experiment_mpraset)))

		# experiment_mpralm_fit <- mpralm_rna(
		experiment_mpralm_fit <- mpralm(
			object = experiment_mpraset,
			design = experiment_design, 
			normalize = TRUE,  aggregate = "mean", 
			block = experiment_block,
			model_type = "corr_groups", plot = TRUE)

		experiment_mpralm_toptab <- topTable(
			experiment_mpralm_fit, 
			coef = 2, number = Inf, confint = TRUE)
		names(experiment_mpralm_toptab) <- paste("mpralm.", 
			names(experiment_mpralm_toptab), sep="")

		experiment_mpralm_final <- experiment_mpralm_final %>% left_join(
			as_tibble(cbind(enframe(rownames(experiment_mpralm_toptab), 
				name=NULL, value="experiment_eid"), experiment_mpralm_toptab)))
	}

	print("fet method")
	# fisher's exact test
	if (chisq_mode) {
		names_experiment_mpralm_final <- names(experiment_mpralm_final)
		experiment_mpralm_final$input_allele1 <- experiment_mpralm_final %>% 
			select(contains("input_allele1")) %>% rowSums()
		experiment_mpralm_final$input_allele2 <- experiment_mpralm_final %>% 
			select(contains("input_allele2")) %>% rowSums()
		experiment_mpralm_final$output_allele1 <- experiment_mpralm_final %>% 
			select(contains("output_allele1")) %>% rowSums()
		experiment_mpralm_final$output_allele2 <- experiment_mpralm_final %>% 
			select(contains("output_allele2")) %>% rowSums() 
		experiment_mpralm_final_fisher <- experiment_mpralm_final %>% 
			select(output_allele1, output_allele2, input_allele1, input_allele2) %>% 
			transpose() %>% lapply(function(x) {x %>% unlist() %>% 
				matrix(nrow=2) %>% fisher.test()})
		experiment_mpralm_final$chisq.P.Value <- 
			lapply(experiment_mpralm_final_fisher, function(x) {
				x$p.value}) %>% unlist()
		experiment_mpralm_final$chisq.adj.P.Val <- 
			p.adjust(experiment_mpralm_final$chisq.P.Value, method="BH")
		experiment_mpralm_final$chisq.logFC <- 
			lapply(experiment_mpralm_final_fisher, function(x) {
				log2(x$estimate)}) %>% unlist()
		experiment_mpralm_final$chisq.CI.L <- 
			lapply(experiment_mpralm_final_fisher, function(x) {
				log2(x$conf.int[1])}) %>% unlist()
		experiment_mpralm_final$chisq.CI.R <- 
			lapply(experiment_mpralm_final_fisher, function(x) {
				log2(x$conf.int[2])}) %>% unlist()
		experiment_mpralm_final <- experiment_mpralm_final %>% 
			dplyr::select(names_experiment_mpralm_final,
				chisq.logFC, chisq.CI.L, chisq.CI.R, chisq.P.Value, chisq.adj.P.Val)
	}

	return(experiment_mpralm_final)
}
