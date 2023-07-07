#!/bin/bash

#datalad clone ria+ssh://audluo@bblsub.pmacs.upenn.edu:/static/LINC_HCPD#~XCP_zipped /cbica/projects/network_replication/input/HCPD/datalad_xcp
# --description "cloned RBC production HCP-D xcp  output_ria into network_replication/input/HCPD/datalad_xcp"
​
 
cd /cbica/projects/network_replication/input/HCPD/datalad_xcp  ###change directory to dir with the datalad zips
for file in sub*zip ; do  ###for every subject zip file
sub=${file%_*}  ###get just the sub-id (extract from the sub-number from the entire file name)
if ! [ -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub ] ###if I have not already extracted the output and put it in this folder....
then

datalad get $file  ###get the zip


#get fluctuation amplitude data
mkdir /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub  ###make dir to put extracted task and rest timeseries

###this unzips just the one file you want and puts it in the -d directory
###gets task-carit, task-emotion, task-guessing, task-rest for all 4 atlases and in both AP and PA
unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Glasser*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub
unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Gordon*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub
unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Schaefer217*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub
unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Schaefer417*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub


#get qc files
unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task-rest*fsLR*qc.csv" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/qc_files/  ###this unzips just the one file you want and puts it in the -d director

datalad drop $file   ###drop the zip
fi
done
 



 