#!/bin/bash


cd /cbica/projects/network_replication/input/PNC/datalad_xcp/xcp_abcd  
for html_file in sub*.html; do   
    sub=${html_file%%[._]*}   
    if ! [ -d /cbica/projects/network_replication/input/PNC/datalad_xcp/$sub ]; then

        mkdir /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub

        git annex get $sub/*/*/*Glasser*ptseries.nii  # get idemo, frac2back, and rest
        cp $sub/*/*/*Glasser*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
        git annex drop $sub/*/*/*Glasser*ptseries.nii

        git annex get $sub/*/*/*Gordon*ptseries.nii   
        cp $sub/*/*/*Gordon*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
        git annex drop $sub/*/*/*Gordon*ptseries.nii

        git annex get $sub/*/*/*Schaefer217*ptseries.nii   
        cp $sub/*/*/*Schaefer217*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
        git annex drop $sub/*/*/*Schaefer217*ptseries.nii

        git annex get $sub/*/*/*Schaefer417*ptseries.nii   
        cp $sub/*/*/*Schaefer417*ptseries.nii /cbica/projects/network_replication/input/PNC/pnc_xcp/$sub
        git annex drop $sub/*/*/*Schaefer417*ptseries.nii

        git annex get $sub/*/*/*framewisedisplacement_bold.tsv
        cp $sub/*/*/*framewisedisplacement_bold.tsv /cbica/projects/network_replication/input/PNC/pnc_xcp/qc_files/
        git annex drop $sub/*/*/*framewisedisplacement_bold.tsv

         # get qc files
        git annex get $sub/*/*/*fsLR_desc-qc_bold.csv
        cp $sub/*/*/*fsLR_desc-qc_bold.csv /cbica/projects/network_replication/input/PNC/pnc_xcp/qc_files
        git annex drop $sub/*/*/*fsLR_desc-qc_bold.csv

    fi
done

