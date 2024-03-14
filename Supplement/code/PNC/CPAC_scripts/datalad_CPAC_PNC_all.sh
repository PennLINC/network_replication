#!/bin/bash

# cd ~/revisions/input/PNC/gitrepo_CPAC
# git clone git@github.com:ReproBrainChart/PNC_CPAC.git
cd /cbica/projects/network_replication/revisions/input/PNC/gitrepo_CPAC/PNC_CPAC/cpac_RBCv0  ###change directory to dir with the datalad folders
  
for file in sub* ; do  ###for every subject zip file
  participant=${file}
  echo "Datalad getting $participant files"
    if ! [ -d /cbica/projects/network_replication/revisions/input/PNC/CPAC_subfiles/$participant ] ###if I have not already extracted the output and put it in this folder....
    then
    mkdir /cbica/projects/network_replication/revisions/input/PNC/CPAC_subfiles/$participant ###make dir to put extracted task and rest timeseries
    datalad get $participant/*/func/*Schaefer2018p200n17*aCompCor_desc-Mean_timeseries.1D # non-GSR timeseries for task and rest
    datalad get $participant/*/func/*Schaefer2018p200n17*36Parameter_desc-PartialNilearn_correlations.tsv # partial correlations for both task and rest scans (just in case)
    datalad get $participant/*/func/*powerParams_motion.txt # qc for task and rest scans 

    cp $participant/*/func/*Schaefer2018p200n17*aCompCor_desc-Mean_timeseries.1D $participant/*/func/*Schaefer2018p200n17*36Parameter_desc-PartialNilearn_correlations.tsv $participant/*/func/*powerParams_motion.txt /cbica/projects/network_replication/revisions/input/PNC/CPAC_subfiles/$participant 
    fi

done 
datalad drop . 
    
 