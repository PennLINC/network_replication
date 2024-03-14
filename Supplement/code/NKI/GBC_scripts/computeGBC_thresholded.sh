#!/bin/bash


# cd /cbica/projects/network_replication/revisions/code/NKI/GBC_scripts
# qsub -l h_vmem=25G,s_vmem=25G computeGBC_thresholded.sh
Rscript --save /cbica/projects/network_replication/revisions/code/utils/computeGBC_absvalue_thresholded.R "NKI" "thresholded"