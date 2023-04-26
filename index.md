<br>
<br>

# Functional Connectivity Development Aligns with the Sensorimotor-Association Cortical Axis in Four Independent Datasets 

Cortical maturation has been posited to be organized along the sensorimotor-association axis, a hierarchical axis of brain organization that spans from unimodal sensorimotor cortices to transmodal association cortices. In this preregistered study, we used four large-scale datasets to investigate whether the development of functional connectivity reliably varies along the sensorimotor-association axis during childhood through adolescence. To do so, we examined four datasets that included youth ages 5-22: the Philadelphia Neurodevelopmental Cohort (N=1207), Nathan Kline Institute-Rockland (N=381), Human Connectome Project: Development (N=625), and Healthy Brain Network (N=1379). Across all datasets, the spatial patterning of connectivity at region aligned with the S-A axis through development. Specifically, global connectivity in sensorimotor regions increased, whereas global connectivity declined in association cortices. Convergent findings across four independent datasets robustly establish that the sensorimotor-association axis is not only a major axis of brain organization, but also encodes the dominant pattern of functional connectivity development.

### Project Lead
Audrey Luo

### Faculty Lead
Theodore D. Satterthwaite

### Analytic Replicator
Valerie J. Sydnor

### Collaborators 
Valerie J. Sydnor, Arielle S. Keller, Aaron F. Alexander-Bloch, Matthew Cieslak, Sydney Covitz, Andrew Chen, Eric Feczko, Alexandre R. Franco, Raquel E. Gur, Ruben C. Gur, Audrey Houghton, Fengling Hu, Gregory Kiar, Bart Larsen, Adam Pines, Giovanni Salum, Tinashe Tapera, Ting Xu, Chenying Zhao, Damien A. Fair, Michael P. Milham, Theodore D. Satterthwaite

### Project Start Date
December 2021

### Current Project Status
In preparation

### Datasets
RBC PNC (Health excude), NKI, HCP-D, and HBN

### Github Repository
<https://github.com/PennLINC/network_replication>

### Conference Presentations
 
- Poster presented at Flux Congress, September 2022.
- Poster to be presented at Organization for Human Brain Mapping, July 2023

### Cubic Project Directory
`/cbica/projects/network_replication`
<br>

`software/`: project software 

`Rscripts/`: contains functions for main and sensitivity analyses as well as Rmd files for each step of the analytic workflow (see below) for each dataset. 

`atlases/dlabel/*.nii` : Cortical parcellations for Schaefer 200 (7 and 17 network), Schaefer 400 (7 and 17 network), Gordon, and HCP-MMP 

`atlases/parcellations/*regionlist_final.csv`: parcel labels for each cortical parcellation (see below) used for study

`atlases/edge/*_edge.csv`: edge names used for edge-level analysis for each cortical parcellation

`input/<dataset>/datalad_xcp/`: MRI data for each dataset pulled via datalad get

`input/<dataset>/<dataset>_xcp/`: only fMRI data in fsLR space and qc data for each subject


`input/<dataset>/connMatricesData/connectivity_matrices`: connectivity matrices derived from concatenated task and rest fMRI scans 
 
<br>

Demographics .csv's all live in `input/<dataset>/sample_selection` but have different file names: 
* PNC: `pnc_participants.tsv` 
* NKI: `nki_participants.tsv`
* HCP-D: `hcpd_demographics.csv`
* HBN: `participants.tsv`

<br>

Final sample lists for each dataset all live in `input/<dataset>/sample_selection` but have different file names:
* PNC: `PNC_demographics_finalsample_20230103.csv` 
* NKI: `NKI_demographics_finalsample_20221219.csv`
* HCP-D: `HCPD_demographics_finalsample_20221226.csv`
* HBN: `HBN_demographics_finalsample_202230226.csv`

<br>

`output/<dataset>/<functional_connectivity_metric>/GAM`: GAM results. Includes effect sizes, p-values, fitted- values, smooth estimates. Outputs for HBN and HCP-D include covbat harmonized outputs. 
<br>
<br>

# CODE DOCUMENTATION  

All project analyses are described below along with the corresponding code on Github. The following outline describes the order of the analytic workflow:

1. Parcellating the sensorimotor-association (S-A) axis  
2. Sample selection for each dataset: PNC (discovery), NKI, HCP-D, and HBN (replication)  
3. Constructing connectivity matrices for each dataset 

* **Main analysis:**
    + *Data*: concatenated task and resting-state fMRI
    + *Cortical parcellation*: [Schaefer 200 atlas](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/Schaefer2018_200Parcels_17Networks_order.dlabel.nii)
         + *Network solution*: [7 Network](https://github.com/ThomasYeoLab/CBIG/blob/6d1400a2d643261246f6b042e7ef5fbe417506cd/utilities/matlab/FC/CBIG_ReorderParcelIndex.m) 

* **Sensitivity analysis:** 
    + *Data*: resting-state fMRI only 
    + *Cortical parcellation*: [Schaefer 400](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/Schaefer2018_400Parcels_17Networks_order.dlabel.nii), [HCP multimodal](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/glasser_space-fsLR_den-32k_desc-atlas.dlabel.nii), [Gordon](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/gordon_space-fsLR_den-32k_desc-atlas.dlabel.nii)
        + *Network solution*: 7 Network and 17 Network (Schaefer atlases)

4. Image harmonization: applying [covbat-gam](https://github.com/andy1764/ComBatFamily) to multi-site data (HCP-D and HBN)
5. Quantification of functional connectivity metrics: global brain connectivity, between- and within-network connectivity
6. Fitting generalized additive models (GAMs) 
7. Characterization of relationships between functional connectivity metrics, age, and the S-A axis

 
<br>

### 1. Parcellating the sensorimotor-association (S-A) axis  
We parcellated the fslr/cifti [Sensorimotor-Association Axis](https://github.com/PennLINC/S-A_ArchetypalAxis/blob/main/FSLRVertex/SensorimotorAssociation_Axis.dscalar.nii) with cortical atlases (Schaefer 200, Schaefer 400, Gordon, and HCP-MMP) using [/Rscripts/Generate_input/1_parcellate_SAaxis.Rmd](https://github.com/PennLINC/network_replication/blob/main/Generate_input/1_parcellate_SAaxis.Rmd).


### 2. Sample selection for each dataset: PNC (discovery), NKI, HCP-D, and HBN (replication)  
The final samples for each dataset were consructed using /Rscripts/<dataset>/QC_scripts/<dataset>_SampleSelection.Rmd. Links to the corresponding github code and descriptions of each final sample are as follows: 

* PNC: [PNC_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/PNC_scripts/QC_scripts/PNC_SampleSelection.Rmd)

        1) Original sample: N=1559, 4546 scans total, ages 8-22
        2) Exclude participants with medical conditions affecting brain function, gross neurological abnormalities, and psychoactive medical medications: N=1413, 4136 scans 
        3) Include passing T1 QC: N=1374, 4041 scans
        4) Include meanFD < 0.3: N=1262, 3365 scans
        5) Include scans with at least 7 minutes of scan time: N= 1207,  3310 scans (final sample), 646 females. Age: mean = 15, SD = 3.3. 
* NKI: [NKI_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/NKI_scripts/QC_scripts/NKI_SampleSelection.Rmd)

        1) Original sample: N=1268, ages 6-85, 6226 scans total
        2) Include ages 6-22: N=424, 2570 scans
        3) Include passing T1 QC: N=402, 2281 scans. 
        - delete all scans from a given session if fail T1 
        4) Include meanFD < 0.3: N=386, 1816 scans
        5) Choose the session that has the most scans surviving the head motion exclusion: N=386, 998 scans
        6) Include scans with at least 7 minutes of scan time: N=381, 993 scans (final sample),  177 females. Age: mean=14.5, SD=4.4

* HCP-D: [HCPD_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/HCPD_scripts/QC_scripts/HCPD_SampleSelection.Rmd)

        1) Original sample: N=652, ages 5-21, 5716 scans total
        2) Exclude participants with medical conditions affecting brain function, gross neurological abnormalities: N=631, 5527 scans
        3) Include passing T1 QC: all scans in dataset have survived T1 QC already 
        4) Include meanFD < 0.3: N=629, 5165 scans 
        5) Include scans with at least 7 minutes of scan time: N=625, 5159 scans (final sample), 337 females. Age: mean=14.5, SD=4.1

* HBN: [HBN_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/HBN_scripts/QC_scripts/HBN_SampleSelection.Rmd)

        1) Original sample: N= 2255, ages 5-21, 6915 scans total
        2) Exclude participants with medical conditions affecting brain function, gross neurological abnormalities: no medical exclusion in HBN
        3) Include passing T1 QC: pending T1 QC from RBC
        4) Include meanFD < 0.3: N= 1649, 3964 scans 
        5) Include scans with at least 7 minutes of scan time: N= 1438, 3939 scans (final sample), 546 females
        6) Exclude participants with missing demographics data (i.e. age and sex): N= 1380, 3767 scans, 546 females. Age: mean=11.6, sd= 3.7

### 3. Constructing connectivity matrices for each dataset

fMRIPrep 20.2.3 (PNC and NKI) and 22.0.2 (HCP-D and HBN) were run with the following parameters:

```bash
$  /usr/local/miniconda/bin/fmriprep inputs/data prep participant -w /scratch/rbc/SGE_820171/$subid/ds/.git/tmp/wkdir --n_cpus 1 --stop-on-first-crash --fs-license-file code/license.txt --skip-bids-validation --bids-filter-file /scratch/rbc/SGE_820171/$subid.json --output-spaces MNI152NLin6Asym:res-2 --participant-label $subid --force-bbr --cifti-output 91k -v -v
```

xcp-d 0.0.8 (NKI) and 0.3.2 (PNC, HCP-D, HBN) were run with the following parameters: 

```bash
$ /usr/local/miniconda/bin/xcp_abcd inputs/data/fmriprep xcp participant --despike --nthreads 1 --lower-bpf 0.01 --upper-bpf 0.08 --participant_label $subid -p 36P -f 10 --cifti

```


### 4. Image harmonization: applying [covbat-gam](https://github.com/andy1764/ComBatFamily) to multi-site data (HCP-D and HBN)
### 5. Quantification of functional connectivity metrics 
### 6. Fitting generalized additive models (GAMs) 
### 7. Characterization of relationships between functional connectivity metrics, age, and the S-A axis

 



### Data Interpretation and Visualization 


### Sensitivity Analyses