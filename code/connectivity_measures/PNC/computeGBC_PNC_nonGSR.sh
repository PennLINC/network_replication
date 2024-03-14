#!/bin/bash


# cd /cbica/projects/network_replication/manuscript/code/7_connectivity_measures/PNC
# qsub -l h_vmem=25G,s_vmem=25G computeGBC_PNC_nonGSR.sh

singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_0.0.4.sif \
    Rscript --save /cbica/projects/network_replication/manuscript/code/scripts/sensitivity_analyses/computeGBC_nonGSR.R