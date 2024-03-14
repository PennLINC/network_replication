#!/bin/bash

#cd /cbica/projects/network_replication/revisions/code/HBN/ConnMatrices_scripts
#qsub -l h_vmem=64G,s_vmem=64G makeConnMatrices_absvalue.sh
Rscript --save /cbica/projects/network_replication/revisions/code/utils/makeConnMatrices_thresholded_absvalue.R "HBN" "absvalue"
