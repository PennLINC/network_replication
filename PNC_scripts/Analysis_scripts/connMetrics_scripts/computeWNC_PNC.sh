#!/bin/bash



# cd /cbica/projects/network_replication/Rscripts/PNC_scripts/Analysis_scripts/connMetrics_scripts/
# qsub -l h_vmem=25G,s_vmem=25G computeWNC_edge.sh
Rscript --save /cbica/projects/network_replication/Rscripts/functions/main_analyses/computeWNC.R "PNC"
