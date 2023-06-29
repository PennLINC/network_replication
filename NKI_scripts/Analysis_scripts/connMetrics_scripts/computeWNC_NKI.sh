#!/bin/bash

# /cbica/projects/network_replication/Rscripts/NKI_scripts/Analysis_scripts/connMetrics_scripts
# qsub -l h_vmem=24G,s_vmem=24G computeWNC_NKI.sh

Rscript --save /cbica/projects/network_replication/Rscripts/functions/main_analyses/computeWNC.R "NKI"
