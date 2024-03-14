#!/bin/bash

# cd /cbica/projects/network_replication/revisions/code/PNC/CPAC_scripts/nonGSR
# qsub -l h_vmem=25G,s_vmem=25G PNC_makeConnMatrices_nonGSR.sh
Rscript --save /cbica/projects/network_replication/revisions/code/PNC/CPAC_scripts/nonGSR/PNC_makeConnMatrices_nonGSR.R

