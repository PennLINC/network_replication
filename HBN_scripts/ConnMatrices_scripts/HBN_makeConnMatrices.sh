#!/bin/bash


# cd /cbica/projects/network_replication/Rscripts/HBN_scripts/ConnMatrices_scripts/
# qsub -l h_vmem=25G,s_vmem=25G HBN_makeConnMatrices.sh

Rscript --save /cbica/projects/network_replication/Rscripts/HBN_scripts/ConnMatrices_scripts/HBN_makeConnMatrices.R
