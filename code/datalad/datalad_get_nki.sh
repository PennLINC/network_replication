#!/bin/bash
 
cd /cbica/projects/network_replication/input/NKI/datalad_xcp  
for file in sub*zip ; do   
    sub=${file%_ses*}
    if ! [ -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub ]; then

        datalad get $file  ###get the zip

        # get resting state MRI data
        mkdir /cbica/projects/network_replication/input/NKI/nki_xcp/$sub  # make dir to put extracted task and rest timeseries

        # this unzips just the one file you want and puts it in the -d directory
        # gets task-rest for all 4 atlases

        unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Glasser*_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
        unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Gordon*_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
        unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Schaefer217_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub
        unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*rest_acq-*_space-fsLR_atlas*Schaefer417_den-91k_bold.ptseries.nii" -d /cbica/projects/network_replication/input/NKI/nki_xcp/$sub

        # get qc files
        unzip -j "$file" "xcp_abcd/$sub/ses-*/func/*_task-rest*fsLR*qc*.csv" -d /cbica/projects/network_replication/input/NKI/nki_xcp/qc_files  ###this unzips just the one file you want and puts it in the -d director

        datalad drop $file   # drop the zip
    fi
done
 