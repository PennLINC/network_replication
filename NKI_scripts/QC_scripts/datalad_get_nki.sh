#!/bin/bash

#datalad clone ria+ssh://audluo@bblsub.pmacs.upenn.edu:/static/LINC_NKI#~XCP_zipped datalad_xcp 
# --description "cloned RBC production NKI xcp  output_ria into network_replication/input/NKI/datalad_xcp"
​
cd /cbica/projects/network_replication/input/NKI/datalad_xcp  ###change directory to dir with the datalad zips
for file in sub*zip ; do  ###for every subject zip file
sub=${file%_ses*}
if ! [ -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub ] ###if I have not already extracted the output and put it in this folder....
then

datalad get $file  ###get the zip

#get resting state MRI data
mkdir /cbica/projects/network_replication/input/NKI/nki_xcp/$sub  ###make dir to put extracted task and rest timeseries

###this unzips just the one file you want and puts it in the -d directory
###gets task-rest for all 4 atlases

unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Glasser*_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Gordon*_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Schaefer217_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Schaefer417_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub


#get qc files
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*_task-rest*fsLR*qc*.csv" -d /cbica/projects/network_replication/input/NKI/nki_xcp/qc_files  ###this unzips just the one file you want and puts it in the -d director

datalad drop $file   ###drop the zip
fi
done
 