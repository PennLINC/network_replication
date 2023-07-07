#!/bin/bash

# for submitting this Rscript job onto cubic, the following command was used in terminal:
# cd /cbica/projects/network_replication/Rscripts/NKI_scripts/ConnMatrices_scripts
# qsub -l h_vmem=25G,s_vmem=25G NKI_makeConnMatrices.sh
Rscript --save /cbica/projects/network_replication/Rscripts/NKI_scripts/ConnMatrices_scripts/NKI_makeConnMatrices.R

