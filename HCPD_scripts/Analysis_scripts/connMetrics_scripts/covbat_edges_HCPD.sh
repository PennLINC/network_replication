#!/bin/bash

singularity run --cleanenv \
    /cbica/projects/network_replication/software/docker/r-packages-for-cubic_r4.1.2forNetworkReplication.sif \
    Rscript --save /cbica/projects/network_replication/Rscripts/functions/main_analyses/covbat_edges.R "HCPD"

 