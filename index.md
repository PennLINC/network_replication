<br>
<br>

# Functional Connectivity Development Aligns with the Sensorimotor-Association Cortical Axis in Four Independent Datasets 

*Cortical maturation has been posited to be organized along the sensorimotor-association axis, a hierarchical axis of brain organization that spans from unimodal sensorimotor cortices to transmodal association cortices. In this preregistered study, we used four large-scale datasets to investigate whether the development of functional connectivity reliably varies along the sensorimotor-association axis during childhood through adolescence. To do so, we examined four datasets that included youth ages 5-22: the Philadelphia Neurodevelopmental Cohort (N=1207), Nathan Kline Institute-Rockland (N=381), Human Connectome Project: Development (N=625), and Healthy Brain Network (N=1379). Across all datasets, the spatial patterning of connectivity at region aligned with the S-A axis through development. Specifically, global connectivity in sensorimotor regions increased, whereas global connectivity declined in association cortices. Convergent findings across four independent datasets robustly establish that the sensorimotor-association axis is not only a major axis of brain organization, but also encodes the dominant pattern of functional connectivity development.* 

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

`atlases/parcellations/*regionlist_final.csv`: 
 parcel labels for each cortical parcellation (see below) used for study

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
### 2. Sample selection for each dataset: PNC (discovery), NKI, HCP-D, and HBN (replication)  
### 3. Constructing connectivity matrices for each dataset
### 4. Image harmonization: applying [covbat-gam](https://github.com/andy1764/ComBatFamily) to multi-site data (HCP-D and HBN)
### 5. Quantification of functional connectivity metrics 
### 6. Fitting generalized additive models (GAMs) 
### 7. Characterization of relationships between functional connectivity metrics, age, and the S-A axis

 
