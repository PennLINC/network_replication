#!/bin/bash

# cd /cbica/projects/network_replication/Rscripts/HBN_scripts/Analysis_scripts/connMetrics_scripts/
# qsub -l h_vmem=25G,s_vmem=25G computeGBC_HBN.sh
Rscript --save /cbica/projects/network_replication/Rscripts/functions/main_analyses/computeGBC.R "HBN"
