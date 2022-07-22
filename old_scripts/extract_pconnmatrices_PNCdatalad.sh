cd /cbica/projects/spatiotemp_dev_plasticity/RBC_datalad/PNC/xcp

for file in sub* ; do
sub=${file%_*}

mkdir /cbica/projects/spatiotemp_dev_plasticity/GBC/PNC/pconn_files/$sub

#extract pconn.nii files
unzip -j "$file" "xcp_abcd/$sub/ses-PNC1/func/*ses-PNC1_task-rest_acq-singleband_space-fsLR_atlas*91k_bold.pconn.nii" -d /cbica/projects/spatiotemp_dev_plasticity/GBC/PNC/pconn_files/$sub


done

