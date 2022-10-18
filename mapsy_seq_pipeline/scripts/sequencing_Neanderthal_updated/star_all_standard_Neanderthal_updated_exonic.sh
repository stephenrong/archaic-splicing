#!/bin/sh

# Postprocess STAR alignment of updated sequencing data.

# preamble
Node=node01
Script=../sequencing_helper_files/star_all_standard_neanderthal.sh
InDirectory=/datasets2/lab_sequencing/Stephen_Chris_experiments/Neanderthal_mapsy_experiments/Neanderthal-HEK_293T-final_expanded_exonic/fastq
OutDirectory=../../data/Neanderthal-HEK_293T-final_expanded_exonic

# input
Index=../../data/custom_star_index/neanderthal_updated_input
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2A-input-a N2A-input-a
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2A-input-b N2A-input-b
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2A-input-c N2A-input-c

sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2B-input-a N2B-input-a
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2B-input-b N2B-input-b
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2B-input-c N2B-input-c

# output
Index=../../data/custom_star_index/neanderthal_updated_output
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2A-output-a N2A-output-a
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2A-output-b N2A-output-b
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2A-output-c N2A-output-c

sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2B-output-a N2B-output-a
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2B-output-b N2B-output-b
sbatch -w $Node $Script $Index $InDirectory $OutDirectory N2B-output-c N2B-output-c
