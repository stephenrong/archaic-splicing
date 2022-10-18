#!/bin/bash

# Process updated sequencing data.

sh star_all_standard_Neanderthal_updated_exonic.sh
sh star_all_standard_process_Neanderthal_updated_exonic.sh
multiqc ../../data/Neanderthal-HEK_293T-final_expanded_exonic/star_align/N* -o ../../data/Neanderthal-HEK_293T-final_expanded_exonic/star_align/multiqc/ -v -f
