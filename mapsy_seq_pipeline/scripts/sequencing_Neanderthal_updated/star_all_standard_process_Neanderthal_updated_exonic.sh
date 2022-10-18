#!/bin/sh

# STAR alignment of pilot sequencing data.

# preamble
Node=node01
Script=../sequencing_helper_files/star_all_standard_process_neanderthal.sh
OutDirectory=../../data/Neanderthal-HEK_293T-final_expanded_exonic

# input
sbatch -w $Node $Script $OutDirectory N2A-input-a
sbatch -w $Node $Script $OutDirectory N2A-input-b
sbatch -w $Node $Script $OutDirectory N2A-input-c

sbatch -w $Node $Script $OutDirectory N2B-input-a
sbatch -w $Node $Script $OutDirectory N2B-input-b
sbatch -w $Node $Script $OutDirectory N2B-input-c

# output
sbatch -w $Node $Script $OutDirectory N2A-output-a
sbatch -w $Node $Script $OutDirectory N2A-output-b
sbatch -w $Node $Script $OutDirectory N2A-output-c

sbatch -w $Node $Script $OutDirectory N2B-output-a
sbatch -w $Node $Script $OutDirectory N2B-output-b
sbatch -w $Node $Script $OutDirectory N2B-output-c
