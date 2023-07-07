#!/bin/bash

# cd /cbica/projects/network_replication/software/docker
# singularity pull docker://audreycluo/r-packages-for-cubic:r_forNetworkReplication0.0.3

# cd /cbica/projects/network_replication/Rscripts/HCPD_scripts/Analysis_scripts/connMetrics_scripts/
# qsub -l h_vmem=64G,s_vmem=64G covbat_edges_HCPD.sh

singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_r_forNetworkReplication0.0.3.sif \
    Rscript --save /cbica/projects/network_replication/Rscripts/functions/main_analyses/covbat_edges.R "HCPD"

 