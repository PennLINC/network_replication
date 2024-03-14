#!/bin/bash
 

# cd /cbica/projects/network_replication/manuscript/code/covbat_harmonization/HBN
# qsub -l h_vmem=64G,s_vmem=64G covbat_edges_HBN.sh

singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_0.0.4.sif \
    Rscript --save /cbica/projects/network_replication/manuscript/code/scripts/main_analyses/covbat_edges.R "HBN"

 