#!/bin/bash

#cd /cbica/projects/network_replication/revisions/code/NKI/ConnMatrices_scripts
#qsub -l h_vmem=64G,s_vmem=64G NKI_mean_SM_matrix.sh
Rscript --save /cbica/projects/network_replication/revisions/code/utils/somatomotor_FCmap.R "NKI"  
