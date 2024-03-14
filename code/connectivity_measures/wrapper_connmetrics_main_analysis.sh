#!/bin/bash

datasets=("PNC" "NKI" "HCPD" "HBN")
measures=("GBC" "BNC" "WNC" "Edge")

script_dir="/cbica/projects/network_replication/manuscript/code/connectivity_measures"
outputs_dir_logs="/cbica/projects/network_replication/manuscript/logs"

# Loop through each dataset
for dataset in "${datasets[@]}"; do
    # Loop through each measure
    for measure in "${measures[@]}"; do
        # Submit job for the current dataset and measure
        qsub -l h_vmem=16G,s_vmem=16G \
            -N compute${measure}_${dataset} \
            -b y \
            -V \
            -j n \
            -o ${outputs_dir_logs} \
            -e ${outputs_dir_logs} \
            ${script_dir}/${dataset}/compute${measure}_${dataset}.sh
    done
done