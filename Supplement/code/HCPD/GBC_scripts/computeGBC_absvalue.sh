#!/bin/bash


# cd /cbica/projects/network_replication/revisions/code/HCPD/GBC_scripts
# qsub -l h_vmem=25G,s_vmem=25G computeGBC_absvalue.sh
Rscript --save /cbica/projects/network_replication/revisions/code/utils/computeGBC_absvalue_thresholded.R "HCPD" "absvalue"