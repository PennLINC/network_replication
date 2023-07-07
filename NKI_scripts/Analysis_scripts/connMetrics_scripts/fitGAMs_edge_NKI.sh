#!/bin/bash

# cd /cbica/projects/network_replication/software/docker
# singularity pull docker://audreycluo/r-packages-for-cubic:r_forNetworkReplication0.0.3

# cd /cbica/projects/network_replication/Rscripts/NKI_scripts/Analysis_scripts/connMetrics_scripts/
# qsub -l h_vmem=25G,s_vmem=25G fitGAMs_edge_NKI.sh
singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_r_forNetworkReplication0.0.3.sif \
    Rscript --save /cbica/projects/network_replication/Rscripts/functions/main_analyses/edge_fitGAMs.R "NKI"

 
