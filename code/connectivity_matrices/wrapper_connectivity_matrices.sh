#!/bin/bash

 
# make connectivity matrices - main analysis 
datasets=("PNC" "NKI" "HCPD" "HBN")
script_dir="/cbica/projects/network_replication/manuscript/code/connectivity_matrices"
outputs_dir_logs="/cbica/projects/network_replication/manuscript/logs"

# Loop through each dataset
for dataset in "${datasets[@]}"; do
    qsub -l h_vmem=25G,s_vmem=25G \
        -N makeConnMatrices_${dataset} \
        -b y \
        -V \
        -j n \
        -o ${outputs_dir_logs} \
        -e ${outputs_dir_logs} \
        ${script_dir}/${dataset}/${dataset}_makeConnMatrices.sh
done
 