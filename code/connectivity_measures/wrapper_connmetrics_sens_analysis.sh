#!/bin/bash


# rest only analyses
datasets=("PNC" "HCPD" "HBN")
measures=("GBC")
script_dir="/cbica/projects/network_replication/manuscript/code/connectivity_measures"
outputs_dir_logs="/cbica/projects/network_replication/manuscript/logs"

# Loop through each dataset
for dataset in "${datasets[@]}"; do
    # Loop through each measure
    for measure in "${measures[@]}"; do
        # Submit job for the current dataset and measure
        qsub -l h_vmem=25G,s_vmem=25G \
            -N compute${measure}_${dataset}_restOnly \
            -b y \
            -V \
            -j n \
            -o ${outputs_dir_logs} \
            -e ${outputs_dir_logs} \
            ${script_dir}/${dataset}/sensitivity_compute${measure}_${dataset}_restOnly.sh
    done
done


# network pair analysis
datasets=("PNC" "NKI" "HCPD" "HBN")
measure=("NetworkPair")
script_dir="/cbica/projects/network_replication/manuscript/code/connectivity_measures"
outputs_dir_logs="/cbica/projects/network_replication/manuscript/logs"

# Loop through each dataset
for dataset in "${datasets[@]}"; do
    # Submit job for the current dataset and measure
    qsub -l h_vmem=25G,s_vmem=25G \
        -N compute${measure}_${dataset} \
        -b y \
        -V \
        -j n \
        -o ${outputs_dir_logs} \
        -e ${outputs_dir_logs} \
        ${script_dir}/${dataset}/compute${measure}_${dataset}.sh
done


# PNC - absolute value, thresholded, nonGSR
dataset="PNC"
measure="GBC"
conn_types=("absvalue" "thresholded" "nonGSR")
for conn_type in "${conn_types[@]}"; do
    # Submit job for the current dataset and measure
    qsub -l h_vmem=25G,s_vmem=25G \
        -N compute${measure}_${dataset}_${conn_type} \
        -b y \
        -V \
        -j n \
        -o ${outputs_dir_logs} \
        -e ${outputs_dir_logs} \
        ${script_dir}/${dataset}/compute${measure}_${dataset}_${conn_type}.sh
done