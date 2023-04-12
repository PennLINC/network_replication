#!/bin/bash

#datalad clone ria+ssh://audluo@bblsub.pmacs.upenn.edu:/static/LINC_HBN#~FMRIPREP-FUNC_zipped datalad_fmriprep
 

# --description cloned from PMACS to cubic
# first datalad get, then did the cp as a second step (since i had to enter password for each datalad get for PMACS transfer)
​

cd /cbica/projects/network_replication/input/HBN/datalad_fmriprep  ###change directory to dir with the datalad zips

datalad get sub*.zip -J 3  ###get the zip


for file in sub*zip ; do  ###for every subject zip file
sub=${file%_*}

if ! [ -d /cbica/projects/network_replication/input/HBN/HBN_fmriprep/$sub] ###if I have not already extracted the output and put it in this folder....
then


#get fmri brain mask for each subject
mkdir /cbica/projects/network_replication/input/HBN/HBN_fmriprep/$sub  ###make dir to put extracted fmri brain mask

###this unzips just the one file you want and puts it in the -d directory
unzip -j "$file" "fmriprep/$sub/ses-*/func/*brain_mask*" -d /cbica/projects/network_replication/input/HBN/HBN_fmriprep/$sub
 

echo ${sub}
fi
done

cd /cbica/projects/network_replication/input/HBN/HBN_fmriprep/
rm -rf sub*/*peer* #didn't realize peer wasn't needed until after all files were copied over!

cd qc_files
rm -f sub*/*peer*

# cd /cbica/projects/network_replication/input/HBN/qc_files --description concatenated all the qc files into a csv
# cat sub*/*space-fsLR_den-91k_qc.csv > /cbica/projects/network_replication/input/HBN/sample_selection/HBN_fmriprep_qc.csv

 