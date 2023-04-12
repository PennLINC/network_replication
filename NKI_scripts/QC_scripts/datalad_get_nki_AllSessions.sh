#!/bin/bash

 
#datalad clone ria+ssh://audluo@bblsub.pmacs.upenn.edu:/static/LINC_NKI#~XCP_zipped datalad_xcp 
 
# --description "cloned RBC production NKI xcp  output_ria into network_replication/input/NKI/datalad_xcp"
​
cd /cbica/projects/network_replication/input/NKI/datalad_xcp  ###change directory to dir with the datalad zips
for file in sub*zip ; do  ###for every subject zip file
ses=${file#sub-A*_}
ses=${ses%_xcp-0-0-8.zip}
sub=${file%_ses*}

if ! [ -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub/$ses] ###if I have not already extracted the output and put it in this folder....
then
datalad get $file  ###get the zip

#get resting state MRI data
mkdir /cbica/projects/network_replication/input/NKI/nki_xcp/$sub/$ses  ###make dir to put extracted task and rest timeseries

###this unzips just the one file you want and puts it in the -d directory
###gets task-rest for all 4 atlases

unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Glasser*_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Gordon*_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Schaefer217_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Schaefer417_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub


#get processed timeseries residuals
mkdir /cbica/projects/network_replication/input/NKI/nki_xcp/qc_files/$sub/$ses
unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*_task-rest*fsLR*qc*.csv" -d /cbica/projects/network_replication/input/NKI/nki_xcp/qc_files/$sub  ###this unzips just the one file you want and puts it in the -d director

datalad drop $file   ###drop the zip

fi
done
#rm -r /cbica/projects/network_replication/input/NKI/nki_xcp/sub*/ses*
#rm -r /cbica/projects/network_replication/input/NKI/nki_xcp/qc_files/sub*/ses*

#cd /cbica/projects/network_replication/input/NKI/nki_xcp/qc_files --description concatenated all the qc files into a csv
#cat */*.csv > NKI_xcp_qc_20221217.csv
