#!/bin/bash

#cd /cbica/projects/network_replication/revisions/code/HCPD/ConnMatrices_scripts
#qsub -l h_vmem=64G,s_vmem=64G makeConnMatrices_thresholded.sh
Rscript --save /cbica/projects/network_replication/revisions/code/utils/makeConnMatrices_thresholded_absvalue.R "HCPD" "thresholded"
