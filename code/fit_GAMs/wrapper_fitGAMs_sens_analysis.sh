#!/bin/bash


outputs_dir_logs="/cbica/projects/network_replication/manuscript/logs"
functions_dir="/cbica/projects/network_replication/manuscript/code/scripts/sensitivity_analyses/"

###############################
# Fit GAMs for GBC, rest only #
############################### 
datasets=("PNC" "HCPD" "HBN")
# Loop through each dataset
for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/network_replication/manuscript/code/config_${dataset}.json"
    qsub -l h_vmem=16G,s_vmem=16G \
			-N fitGAMs_GBC_restOnly_${dataset} \
			-b y \
			-V \
			-j n \
			-o ${outputs_dir_logs} \
			-e ${outputs_dir_logs} \
			${functions_dir}/fitGAMs_GBC_restOnly_singularity.sh ${config_file}
done


##############################
# Fit GAMs for Network Pairs #
##############################
datasets=("PNC" "NKI" "HCPD" "HBN")
# Loop through each dataset
for dataset in "${datasets[@]}"; do
    config_file="/cbica/projects/network_replication/manuscript/code/config_${dataset}.json"
    qsub -l h_vmem=16G,s_vmem=16G \
			-N fitGAMs_networkpair_${dataset} \
			-b y \
			-V \
			-j n \
			-o ${outputs_dir_logs} \
			-e ${outputs_dir_logs} \
			${functions_dir}/fitGAMs_networkpair_singularity.sh ${config_file}
done



####################################################
# Fit GAMs for PNC: abs value, thresholded, nonGSR #
####################################################
dataset="PNC"
config_file="/cbica/projects/network_replication/manuscript/code/config_${dataset}.json"
qsub -l h_vmem=16G,s_vmem=16G \
			-N fitGAMs_GBC_absvalue_thresholded_${dataset} \
			-b y \
			-V \
			-j n \
			-o ${outputs_dir_logs} \
			-e ${outputs_dir_logs} \
			${functions_dir}/fitGAMs_GBC_absvalue_thresholded_singularity.sh ${config_file}

qsub -l h_vmem=16G,s_vmem=16G \
			-N fitGAMs_GBC_nonGSR_${dataset} \
			-b y \
			-V \
			-j n \
			-o ${outputs_dir_logs} \
			-e ${outputs_dir_logs} \
			${functions_dir}/fitGAMs_GBC_nonGSR_singularity.sh ${config_file}



