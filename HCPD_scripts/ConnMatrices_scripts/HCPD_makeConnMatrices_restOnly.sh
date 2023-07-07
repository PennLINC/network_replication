#!/bin/bash

# cd /cbica/projects/network_replication/Rscripts/HCPD_scripts/ConnMatrices_scripts/
# qsub -l h_vmem=25G,s_vmem=25G HCPD_makeConnMatrices_restOnly.sh

Rscript --save /cbica/projects/network_replication/Rscripts/HCPD_scripts/ConnMatrices_scripts/HCPD_makeConnMatrices_restOnly.R
