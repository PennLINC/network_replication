#!/bin/bash

# cd /cbica/projects/network_replication/Rscripts/PNC_scripts/ConnMatrices_scripts/
# qsub -l h_vmem=25G,s_vmem=25G PNC_makeConnMatrices.sh
Rscript --save /cbica/projects/network_replication/Rscripts/PNC_scripts/ConnMatrices_scripts/PNC_makeConnMatrices.R

