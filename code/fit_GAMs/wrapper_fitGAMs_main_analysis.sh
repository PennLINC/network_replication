#!/bin/bash

datasets=("PNC" "NKI" "HCPD" "HBN")
outputs_dir_logs="/cbica/projects/network_replication/manuscript/logs"
functions_dir="/cbica/projects/network_replication/manuscript/code/scripts/main_analyses/"

##################################
# Fit GAMs for GBC, BNC, and WNC #
##################################
# Loop through each dataset
for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/network_replication/manuscript/code/config_${dataset}.json"
    qsub -l h_vmem=64G,s_vmem=64G \
                -N fitGAMs_FCmetrics_${dataset} \
                -b y \
                -V \
                -j n \
                -o ${outputs_dir_logs} \
                -e ${outputs_dir_logs} \
                ${functions_dir}/fitGAMs_FCmetrics_singularity.sh ${config_file}
done


######################
# Fit GAMs for edges #
######################
for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/network_replication/manuscript/code/config_${dataset}.json"
    qsub -l h_vmem=64G,s_vmem=64G \
                -N fitGAMs_edges_${dataset} \
                -b y \
                -V \
                -j n \
                -o ${outputs_dir_logs} \
                -e ${outputs_dir_logs} \
                ${functions_dir}/fitGAMs_edges_singularity.sh ${config_file}
done


