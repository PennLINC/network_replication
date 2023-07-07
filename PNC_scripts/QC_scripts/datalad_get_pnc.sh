#!/bin/bash

#datalad clone ria+file:///cbica/projects/RBC/production/PNC/xcp/output_ria#~data /cbica/projects/network_replication/input/PNC/datalad_xcp
# --description "cloned RBC production PNC xcp  output_ria into network_replication/input/PNC/datalad_xcp"
​
# git annex was used in 1/2023 when RBC data still lived on CUBIC

cd /cbica/projects/network_replication/input/PNC/datalad_xcp/xcp_abcd  ###change directory to dir with the datalad zips
for html_file in sub*.html; do  ###for every subject folder ${folder%_*}
sub=${html_file%%[._]*}  ###get just the sub-id (extract from the sub-number from the entire file name)
if ! [ -d /cbica/projects/network_replication/input/PNC/datalad_xcp/$sub ] ###if I have not already extracted the output and put it in this folder....
then

mkdir /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub

git annex get $sub/*/*/*Glasser*ptseries.nii  ###get Glasser for idemo, frac2back, and rest
cp $sub/*/*/*Glasser*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
git annex drop $sub/*/*/*Glasser*ptseries.nii

git annex get $sub/*/*/*Gordon*ptseries.nii  ###get Glasser for idemo, frac2back, and rest
cp $sub/*/*/*Gordon*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
git annex drop $sub/*/*/*Gordon*ptseries.nii

git annex get $sub/*/*/*Schaefer217*ptseries.nii  ###get Glasser for idemo, frac2back, and rest
cp $sub/*/*/*Schaefer217*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
git annex drop $sub/*/*/*Schaefer217*ptseries.nii

git annex get $sub/*/*/*Schaefer417*ptseries.nii  ###get Glasser for idemo, frac2back, and rest
cp $sub/*/*/*Schaefer417*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
git annex drop $sub/*/*/*Schaefer417*ptseries.nii

git annex get $sub/*/*/*framewisedisplacement_bold.tsv
cp $sub/*/*/*framewisedisplacement_bold.tsv /cbica/projects/network_replication/input/PNC/pnc_xcp/qc_files/
git annex drop $sub/*/*/*framewisedisplacement_bold.tsv

git annex get $sub/*/*/*fsLR_desc-qc_bold.csv
cp $sub/*/*/*fsLR_desc-qc_bold.csv /cbica/projects/network_replication/input/PNC/pnc_xcp/qc_files
git annex drop $sub/*/*/*fsLR_desc-qc_bold.csv

fi
done

