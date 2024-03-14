#!/bin/bash

cd /cbica/projects/network_replication/input/HCPD/datalad_xcp  
for file in sub*zip ; do   
    sub=${file%_*}   
    if ! [ -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub ]; then

        datalad get $file  
 
        mkdir /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub  

        # this unzips just the one file you want and puts it in the -d directory
        # gets task-carit, task-emotion, task-guessing, task-rest for all 4 atlases and in both AP and PA
        unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Glasser*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub
        unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Gordon*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub
        unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Schaefer217*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub
        unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task*space-fsLR_atlas*Schaefer417*ptseries.nii" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/$sub


        # get qc files
        unzip -j "$file" "xcp_d/$sub/ses-V1/func/*ses-V1_task-rest*fsLR*qc.csv" -d /cbica/projects/network_replication/input/HCPD/hcpd_xcp/qc_files/   

        datalad drop $file   # drop the zip
    fi
done
 



 