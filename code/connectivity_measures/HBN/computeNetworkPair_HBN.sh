#!/bin/bash


singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_0.0.4.sif \
    Rscript --save /cbica/projects/network_replication/manuscript/code/scripts/sensitivity_analyses/computeNetworkPair_conn.R "HBN" 

 
