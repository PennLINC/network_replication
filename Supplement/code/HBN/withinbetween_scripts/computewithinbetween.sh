#!/bin/bash


# cd /cbica/projects/network_replication/revisions/code/HBN/withinbetween_scripts
# qsub -l h_vmem=64G,s_vmem=64G computewithinbetween.sh
singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_r_forNetworkReplication0.0.3.sif \
    Rscript --save /cbica/projects/network_replication/revisions/code/utils/computewithinbetween_schaefer200x7.R "HBN" 

 
