#!/bin/bash

# submit after jobs from wrapper_connectivity_matrices.sh are finished
dataset="PNC"
script_dir="/cbica/projects/network_replication/manuscript/code/connectivity_matrices/${dataset}"
outputs_dir_logs="/cbica/projects/network_replication/manuscript/logs"
qsub -l h_vmem=25G,s_vmem=25G \
			-N makeConnMatrices_absvalue_${dataset} \
			-b y \
			-V \
			-j n \
			-o ${outputs_dir_logs} \
			-e ${outputs_dir_logs} \
			${script_dir}/${dataset}_makeConnMatrices_absvalue.sh

qsub -l h_vmem=25G,s_vmem=25G \
			-N makeConnMatrices_thresholded_${dataset} \
			-b y \
			-V \
			-j n \
			-o ${outputs_dir_logs} \
			-e ${outputs_dir_logs} \
			${script_dir}/${dataset}_makeConnMatrices_thresholded.sh


