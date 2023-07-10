#!/bin/bash

# cd /cbica/projects/network_replication/Rscripts/HCPD_scripts/Analysis_scripts/connMetrics_scripts/
# qsub -l h_vmem=25G,s_vmem=25G sensitivity_computeGBC_HCPD_restOnly.sh

Rscript --save /cbica/projects/network_replication/Rscripts/functions/sensitivity_analyses/computeGBC_restOnly.R "HCPD"
