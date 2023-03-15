#!/bin/sh

# Create archaic mask BED file
Rscript preprocess_mask.R

# Preprocess 1KGP SNP data (get AC, get B, get mask, convert to hub format, save all SNPs, save AF >= 0.9 SNPs [useful for retrieving lineage-specific SNPs])
Rscript preprocess_1KGP_SNPs.R

# Preprocess archaic data (get B, get mask, convert to hub format, join archaics together, add 1KGP SNP info, add gnomAD SNP info [convert to BED, tabix gnomAD, join gnomAD])
Rscript preprocess_archaics.R
# 	Tabix gnomAD SNPs based on archaic SNP to BED
sh preprocess_archaics_part1_gnomAD.sh
# 	Join gnomAD SNP info to archaic data
Rscript preprocess_archaics_part2_gnomAD.R

# Preprocess old library data (get B, get mask, add 1 KGP SNP info, add archaic SNP info, add gnomAD SNP info [convert to BED, tabix gnomAD, join gnomAD])
Rscript preprocess_library.R
# 	Tabix gnomAD SNPs based on archaic SNP to BED
sh preprocess_library_part1_gnomAD.sh
# 	Join gnomAD SNP info to old library data
Rscript preprocess_library_part2_gnomAD.R

# Preprocess REFisDER data [where ref allele is derived not ancestral] (get B, get mask, convert to hub format, add 1KGP SNP info, add archaic SNP info, add gnomAD SNP info [convert to BED, tabix gnomAD, join gnomAD])
Rscript preprocess_REFisDER.R
# 	Tabix gnomAD SNPs based on archaic SNP to BED
sh preprocess_REFisDER_part1_gnomAD.sh
# 	Join gnomAD SNP info to REFisDER data
Rscript preprocess_REFisDER_part2_gnomAD.R

# Preprocess 1KGP greater 0.9
# 	Tabix gnomAD SNPs based on archaic SNP to BED
sh preprocess_ALL_1KGP_phase3_greater0.9_part1_gnomAD.sh
# 	Join gnomAD SNP info to archaic data
Rscript preprocess_ALL_1KGP_phase3_greater0.9_part2_gnomAD.R

# Preprocess 1KGP greater 0.5
# 	Tabix gnomAD SNPs based on archaic SNP to BED
sh preprocess_ALL_1KGP_phase3_greater0.5_part1_gnomAD.sh
# 	Join gnomAD SNP info to archaic data
Rscript preprocess_ALL_1KGP_phase3_greater0.5_part2_gnomAD.R

# Reprocess to get updated lineage-specific variants and combine with old library data
Rscript preprocess_lineage_part1.R
# 	Join individual files together
Rscript preprocess_lineage_part2.R

# Don't append Rinker et al NDA, RAA, RHA SNPs to final
Rscript preprocess_lineage_part3.R
cp ../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz ../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.txt.gz
