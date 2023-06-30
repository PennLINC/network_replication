#!/bin/bash

#datalad clone ria+ssh://audluo@bblsub.pmacs.upenn.edu:/static/LINC_HBN#~XCP_unzipped datalad_xcp

# --description cloned from PMACS to cubic
# first datalad get, then did the cp as a second step (since i had to enter password for each datalad get for PMACS transfer)
​
cd /cbica/projects/network_replication/input/HBN/datalad_xcp  ###change directory to dir with the datalad zips

datalad get sub*/*/*/*Glasser*ptseries.nii -J 3 ###get Glasser for movies and rest 
datalad get sub*/*/*/*Gordon*ptseries.nii -J 3
datalad get sub*/*/*/*Schaefer217*ptseries.nii -J 3 ###get Schaefer200 for movies and rest
datalad get sub*/*/*/*Schaefer417*ptseries.nii -J 3  ###get Schaefer400 for movies and rest
datalad get sub*/*/*/*space-fsLR_den-91k_qc.csv -J 3 # double check this file

for file in sub* ; do  ###for every subject zip file
sub=${file%_*}  ###get just the sub-id (extract from the sub-number from the entire file name)
if ! [ -d /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub ] ###if I have not already extracted the output and put it in this folder....
then

#datalad get $file  ###get the subject folder

mkdir /cbica/projects/network_replication/input/HBN/HBN_xcp/$sub  ###make dir to put extracted task and rest timeseries

 
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

cd /cbica/projects/network_replication/input/HBN/HBN_xcp/
rm -rf sub*/*peer* # peer isn't needed  

cd qc_files
rm -f sub*/*peer*
 