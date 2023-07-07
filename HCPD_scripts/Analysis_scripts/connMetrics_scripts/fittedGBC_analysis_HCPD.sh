#!/bin/bash



# cd /cbica/projects/network_replication/Rscripts/HCPD_scripts/Analysis_scripts/connMetrics_scripts
# qsub -l h_vmem=25G,s_vmem=25G fittedGBC_analysis_HCPD.sh

Rscript --save /cbica/projects/network_replication/Rscripts/functions/main_analyses/fittedGBC_analysis.R "HCPD" "GBC"
