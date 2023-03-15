#!/bin/sh
sh filter_runall_enrichment_GTEx_QTLs.sh  # run to get hapR2 tables
sh filter_runall_enrichment_SpliceAI.sh  # run annotate SpliceAI first
sh filter_runall_enrichment_VEP.sh  # run hapR2 VEP first

sh enrichment_runall_SpliceAI_v2.sh

sh enrichment_runall_VEP_v2.sh

sh enrichment_runall_GTEx_QTLs_v2.sh
sh enrichment_runall_recmap_access_GTEx_QTLs_v2.sh

sh visualize_GTEx_QTLs_runall_enrichment.sh
sh visualize_GTEx_QTLs_recmap_access_runall_enrichment.sh
