#!/bin/bash
 
cd /cbica/projects/network_replication/input/HBN/datalad_xcp   

datalad get sub*/*/*/*Glasser*ptseries.nii -J 3 # get Glasser for movies and rest 
datalad get sub*/*/*/*Gordon*ptseries.nii -J 3
datalad get sub*/*/*/*Schaefer217*ptseries.nii -J 3 # get Schaefer200 for movies and rest
datalad get sub*/*/*/*Schaefer417*ptseries.nii -J 3  # get Schaefer400 for movies and rest
datalad get sub*/*/*/*space-fsLR_den-91k_qc.csv -J 3  

for file in sub* ; do  
    sub=${file%_*}   
    if ! [ -d /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub ]; then

        mkdir /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub  
        
        cp $sub/*/*/*Glasser*ptseries.nii /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub
        datalad drop $sub/*/*/*Glasser*ptseries.nii  

        cp $sub/*/*/*Gordon*ptseries.nii /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub
        datalad drop $sub/*/*/*Gordon*ptseries.nii

        cp $sub/*/*/*Schaefer217*ptseries.nii /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub
        datalad drop $sub/*/*/*Schaefer217*ptseries.nii

        cp $sub/*/*/*Schaefer417*ptseries.nii /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub
        datalad drop sub*/*/*/*Schaefer417*ptseries.nii  
        
        cp $sub/*/*/*space-fsLR_den-91k_qc.csv /cbica/projects/network_replication/input/HBN/HBN_xcp/qc_files/
        datalad drop $sub/*/*/*fsLR_den-91k_qc.csv

        echo ${sub}
    fi
done
 