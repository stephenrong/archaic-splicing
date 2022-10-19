#!/bin/sh
sh filter_runall_enrichment_GTEx_QTLs.sh
sh filter_runall_enrichment_SpliceAI.sh
sh filter_runall_enrichment_VEP.sh

sh enrichment_runall_GTEx_QTLs_v2.sh
sh enrichment_runall_SpliceAI_v2.sh
sh enrichment_runall_VEP_v2.sh

sh visualize_GTEx_QTLs_runall_enrichment.sh
