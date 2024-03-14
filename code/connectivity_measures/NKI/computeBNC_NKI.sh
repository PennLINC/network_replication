#!/bin/bash


# cd /cbica/projects/network_replication/manuscript/code/7_connectivity_measures/NKI
# qsub -l h_vmem=24G,s_vmem=24G computeBNC_NKI.sh
singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_0.0.4.sif \
    Rscript --save /cbica/projects/network_replication/manuscript/code/scripts/main_analyses/computeBNC.R "NKI"
