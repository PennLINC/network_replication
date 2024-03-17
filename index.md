<br>
<br>

# Development of Functional Connectivity Aligns with the Cortical Hierarchy in Four Datasets 

Human cortical maturation has been posited to be organized along the sensorimotor-association axis, a hierarchical axis of brain organization that spans from unimodal sensorimotor cortices to transmodal association cortices. Here, we investigate the hypothesis that the development of functional connectivity during childhood through adolescence conforms to the cortical hierarchy defined by the sensorimotor-association axis. We tested this pre-registered hypothesis in four large-scale, independent datasets (total n = 3,355; ages 5-23 years): the Philadelphia Neurodevelopmental Cohort (n = 1,207), Nathan Kline Institute-Rockland Sample (n = 397), Human Connectome Project: Development (n = 625), and Healthy Brain Network (n = 1,126). Across datasets, the development of functional connectivity systematically varied along the sensorimotor-association axis. Connectivity in sensorimotor regions increased, whereas connectivity in association cortices declined, refining and reinforcing the cortical hierarchy. These consistent and generalizable results establish that the sensorimotor-association axis of cortical organization encodes the dominant pattern of functional connectivity development.   

*Note*: Global brain connectivity (GBC) and functional connectivity (FC) strength are used interchangeably in this documentation. 

### Project Lead
Audrey C. Luo

### Faculty Lead
Theodore D. Satterthwaite

### Analytic Replicator
Valerie J. Sydnor

### Collaborators 
Valerie J. Sydnor, Adam Pines, Bart Larsen, Aaron F. Alexander-Bloch, Matthew Cieslak, Sydney Covitz, Andrew Chen, Nathalia Bianchini Esper, Eric Feczko, Alexandre R. Franco, Raquel E. Gur, Ruben C. Gur, Audrey Houghton, Fengling Hu,l, Arielle S. Keller, Gregory Kiar, Kahini Mehta, Giovanni A. Salum, Tinashe Tapera, Ting Xu, Chenying Zhao, Damien A. Fair, Taylor Salo, Russell T. Shinohara, Michael P. Milham, Theodore D. Satterthwaite

### Project Start Date
December 2021

### Current Project Status
Accepted in principle at Nature Communications

### Datasets
RBC PNC (Health excude), NKI, HCP-D, and HBN

### Github Repository
<https://github.com/PennLINC/network_replication>
### Slack Channel:
#network_replication 

### Conference Presentations
 
- Poster presented at Flux Congress, September 2022.
- Poster presented at Organization for Human Brain Mapping, July 2023

### Cubic Project Directory Overview
`/cbica/projects/network_replication`
<br>

`/software/`: project software 

`/manuscript/`: final code, input, output, and figures for manuscript. __Note that the github repo contains the content in this directory.__

`/atlases/dlabel/*.nii` : Cortical parcellations for Schaefer 200 (7 and 17 network), Schaefer 400 (7 and 17 network), Gordon, and HCP-MMP 

`/atlases/parcellations/*regionlist_final.csv`: parcel labels for each cortical parcellation (see below) used for study

`/atlases/edge/*_edge.csv`: edge names used for edge-level analysis for each cortical parcellation

`/input/<dataset>/datalad_xcp/`: MRI data for each dataset pulled via datalad get

`/input/<dataset>/<dataset>_xcp/`: only fMRI data in fsLR space and qc data for each subject


`/manuscript/input/<dataset>/connMatricesData/connectivity_matrices`: connectivity matrices derived from concatenated task and rest fMRI scans 
 
<br>

Demographics .csv's all live in `/input/<dataset>/sample_selection` but have different file names: 
* PNC: `pnc_participants.tsv` 
* NKI: `nki_participants.tsv`
* HCP-D: `hcpd_demographics.csv`
* HBN: `HBN_Demographics.tsv`

<br>

Final sample lists for each dataset all live in `/input/<dataset>/sample_selection` but have different file names:
* PNC: `PNC_demographics_finalsample.csv` 
* NKI: `NKI_demographics_finalsample.csv`
* HCP-D: `HCPD_demographics_finalsample.csv`
* HBN: `HBN_demographics_finalsample.csv`

`/manuscript/output/<dataset>/<functional_connectivity_metric>/GAM`: GAM results. Includes effect sizes, p-values, fitted- values, smooth estimates. Outputs for HBN and HCP-D include covbat harmonized outputs. 
<br>
<br>

# CODE DOCUMENTATION  

All project analyses are described below along with the corresponding code on Github. The following outline describes the order of the analytic workflow:

*0.* Get static data from PMACS  
*1.* Parcellating the sensorimotor-association (S-A) axis   
*2.* Formatting parcel labels for different cortical parcellations  
*3.* Creating the spin test parcel rotation matrix for significance testing  
*4.* Constructing the sample for each dataset:  
* *Discovery: PNC* 
* *Replication: NKI, HCP-D, and HBN*    

*5.* Constructing connectivity matrices for each dataset  
*6.* Quantifying functional connectivity metrics: global brain connectivity, between- and within-network connectivity  
*7.* Image harmonization: applying [covbat-gam](https://github.com/andy1764/ComBatFamily) to multi-site data (HCP-D and HBN)  
*8.* Fitting generalized additive models (GAMs) and doing age-resolved analysis  
*9.* Characterizing relationships between functional connectivity metrics, age, and the S-A axis  
*10.* Visualizing results!
 

 
<br>

### 0. Getting static data from PMACS 
Scripts in [/manuscript/code/datalad](https://github.com/PennLINC/network_replication/tree/main/code/datalad) were used to download the static data for PNC, NKI, HCP-D, and HBN used for this project from PMACS using datalad get.

### 1. Parcellating the sensorimotor-association (S-A) axis  
We parcellated the fslr/cifti [Sensorimotor-Association Axis](https://github.com/PennLINC/S-A_ArchetypalAxis/blob/main/FSLRVertex/SensorimotorAssociation_Axis.dscalar.nii) with cortical atlases (Schaefer 200, Schaefer 400, Gordon, and HCP-MMP) using [/manuscript/code/parcellations/parcellate_SAaxis.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/parcellations/parcellate_SAaxis.Rmd).

### 2. Formatting labels for different cortical parcellations

Our primary parcellation utilized the Schaefer 200 atlas; secondary atlases included the Schaefer 400 atlas, the Gordon atlas, and HCP-MMP atlas. For analyses of secondary outcome measures that require community structure, namely average between- and within-network connectivity, we evaluated both the Yeo 7 and 17-network partitions associated with the Schaefer atlas. Results from sensitivity analyses are presented in the supplement.

+ *Cortical parcellations*: [Schaefer 200 atlas](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/Schaefer2018_200Parcels_17Networks_order.dlabel.nii), [Schaefer 400](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/Schaefer2018_400Parcels_17Networks_order.dlabel.nii), [HCP multimodal](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/glasser_space-fsLR_den-32k_desc-atlas.dlabel.nii), [Gordon](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/gordon_space-fsLR_den-32k_desc-atlas.dlabel.nii) (note that some of these files have changed since they were originally downloaded from the xcp-d git repo)

Importantly, the timeseries data for all datasets was available in Schaefer 200x17 rather than Schaefer 200x7. While the two parcellations have the same number of parcels, the ordering is slightly different. Since Schaefer 200x7 was our primary parcellation and network solution, we determined the ordering between Schaefer 200x7 and 200x17 using [this function from the Yeo group](https://github.com/ThomasYeoLab/CBIG/blob/6d1400a2d643261246f6b042e7ef5fbe417506cd/utilities/matlab/FC/CBIG_ReorderParcelIndex.m). This ordering can be found at `/software/schaefer7to17_reordering/index_schaefer200x7to17.csv`. 
 
Parcel labels are formatted and reordered using:
* [/manuscript/code/parcellations
/format_parcel_labels.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/parcellations/format_parcel_labels.Rmd)

For edge-wise analyses, we also created edge labels for each edge, named after the two cortical endpoints using:
* [/manuscript/code/parcellations/edge_labels.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/parcellations/edge_labels.Rmd)


### 3. Creating the spin test parcel rotation matrix for significance testing
For spin tests used in step 9, we first needed to generate rotated permutations of a parcellation. 

This was done using
[/manuscript/code/spin_test_matrix/SpinTest_RotationMatrix.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/spin_test_matrix/SpinTest_RotationMatrix.Rmd)


This Rmd file used [rotate_parcellation](https://github.com/frantisekvasa/rotate_parcellation) to generate 10k rotated permutations of a parcellation, given the coordinates of left and right hemispheres of this parcellation on the sphere.

The output is a vector of 1:Number_of_parcels for each permutation which tries to conserve the relative position of each parcel.

- Schaefer .annot files were downloaded from [https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3). Centroid coordinates were extracted from .annot files. 


- glasser_spherical_coords.csv was downloaded from [https://github.com/frantisekvasa/rotate_parcellation/blob/master/sphere_HCP.txt](https://github.com/frantisekvasa/rotate_parcellation/blob/master/sphere_HCP.txt )

### 4. Sample selection for each dataset: PNC (discovery), NKI, HCP-D, and HBN (replication)  
The final samples for each dataset were constructed using 
```
/manuscript/code/sample_selection/<dataset>_SampleSelection.Rmd
```
Note that `PNC_SampleSelection_nonGSR_CPAC.Rmd` is for supplementary analyses. 

Links to the corresponding github code and descriptions of each final sample are as follows: 

* PNC: [/manuscript/code/sample_selection/PNC_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/sample_selection/PNC_SampleSelection.Rmd)

        1) Original sample: N=1559, 4546 scans total, ages 8-23
        2) Exclude participants with medical conditions affecting brain function, gross neurological abnormalities, and psychoactive medical medications: N=1413, 4136 scans 
        3) Include passing T1 QC: N=1374, 4041 scans
        4) Include meanFD < 0.3: N=1262, 3365 scans
        5) Include scans with at least 7 minutes of scan time: N= 1207,  3310 scans (final sample), 646 females.  
                - range scan time: 3.60 to 33.25 min
                - median: 28.25 min
        Age: mean = 15.4, SD = 3.5
 
        Race
        - white = 551; 45.7%
        - black = 513; 42.5%
        - other (native american + hawaiian pi) = 4 + 1; 0.41%
        - asian = 11; 0.9%
        - missing = 0
        - mixed = 127; 10.5%


* NKI: [/manuscript/code/sample_selection/NKI_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/sample_selection/NKI_SampleSelection.Rmd)

        1) Original sample: N=1288, ages 6-85, 6346 scans total
        2) Include ages 6-22: N=426, 2584 scans
        3) Include passing T1 QC: N=421, 2474 scans. 
        - delete all scans from a given session if fail T1 
        4) Include meanFD < 0.3: N=403, 1918 scans
        5) Choose the session that has the most scans surviving the head motion exclusion: N=403, 1037 scans
        6) Include scans with at least 7 minutes of scan time: N=397, 1031 scans (final sample), 186 females    
                range scan time: 7.75075 to 24.10167 min
                median: 24.10167 min
         Age: mean=14.5, SD=4.4

         Race 
        - White = 258 - 65.0%
        - Asian = 34 - 8.5%
        - Black = 82 - 20.7%
        - Other = 11 - 2.8%
        - Missing = 12 - 3.0%



* HCP-D: [/manuscript/code/sample_selection/HCPD_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/sample_selection/HCPD_SampleSelection.Rmd)

        1) Original sample: N=652, ages 5-22, 5716 scans total
        2) Exclude participants with medical conditions affecting brain function, gross neurological abnormalities: N=631, 5527 scans
        3) Include passing T1 QC: all scans in dataset have survived T1 QC already 
        4) Include meanFD < 0.3: N=629, 5165 scans 
        5) Include scans with at least 7 minutes of scan time: N=625, 5159 scans (final sample), 337 females. Age: mean=14.5, SD=4.1
                - range scan time: 7.466667 to 42.66667 min
                - median scan time: 42.666667 min
        Age: mean=14.5, SD=4.1

        Race:
        - white = 395; 63.2%
        - asian = 48; 7.7%%
        - black = 69; 11%
        - other (native american, hawaiin or pacific islander) = 3; 0.48%
        - missing = 16; 2.6%
        - mixed = 94; 15%

* HBN: [/manuscript/code/sample_selection/HBN_SampleSelection.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/sample_selection/HBN_SampleSelection.Rmd)

        1) original sample: N=2255, ages 5-22, 6915 scans total
        2) exclude participants with medical conditions affecting brain function, gross neurological abnormalities: no medical exclusion in HBN
        3) include passing T1 QC: N=1669, 5111 scans
        4) include meanFD < 0.3: N=1315, 3283 scans
        5) include scans with at least 7 minutes of scan time: N=1126, 3094 scans (final sample), 439 females
                - scan time range: 8.33 to 23.33 min 
                - median: 18.33 min 
        Age: mean = 11.6, SD= 3.6

        Race
        - White = 498; 44.2%
        - Asian = 33; 2.93%
        - Black = 139; 12.34%
        - Other (native american, hawaiin pacific islander) = 29 - 2.6%
        - missing = 6 + 143 = 149; 13.2%  
        - Mixed = 169; 15.0%
        - Hispanic = 109; 9.7%



### 5. Constructing connectivity matrices for each dataset

fMRIPrep 20.2.3 (PNC and NKI) and 22.0.2 (HCP-D and HBN) were run with the following parameters:

```bash
$  /usr/local/miniconda/bin/fmriprep inputs/data prep participant -w /scratch/rbc/SGE_820171/$subid/ds/.git/tmp/wkdir --n_cpus 1 --stop-on-first-crash --fs-license-file code/license.txt --skip-bids-validation --bids-filter-file /scratch/rbc/SGE_820171/$subid.json --output-spaces MNI152NLin6Asym:res-2 --participant-label $subid --force-bbr --cifti-output 91k -v -v
```

xcp-d 0.0.8 (NKI) and 0.3.2 (PNC, HCP-D, HBN) were run with the following parameters: 

```bash
$ /usr/local/miniconda/bin/xcp_abcd inputs/data/fmriprep xcp participant --despike --nthreads 1 --lower-bpf 0.01 --upper-bpf 0.08 --participant_label $subid -p 36P -f 10 --cifti

```


**3a. Main analysis:**  
Vertex-level fMRI timeseries were parcellated with fsLR surface atlases utilizing Connectome Workbench 1.5.0.19. This produced fMRI timeseries within individual cortical regions. The Schaefer200 atlas was used as the primary atlas and Schaefer 400, HCP-MMP, and Gordon atlases were used in sensitivity analyses. Next, parcellated rest and task fMRI timeseries were concatenated and the Pearson correlation between concatenated timeseries was computed for every pair of cortical regions. 

*Data*: task and resting-state fMRI were concatenated using the following code:
+ PNC: [/manuscript/code/connectivity_matrices/PNC/PNC_makeConnMatrices.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/PNC/PNC_makeConnMatrices.R)
+ NKI: [/manuscript/code/connectivity_matrices/NKI/NKI_makeConnMatrices.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/NKI/NKI_makeConnMatrices.R)
+ HCP-D: [/manuscript/code/connectivity_matrices/HCPD/HCPD_makeConnMatrices.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/HCPD/HCPD_makeConnMatrices.R)
+ HBN:  [/manuscript/code/connectivity_matrices/HBN/HBN_makeConnMatrices.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/HBN/HBN_makeConnMatrices.R)

The above Rscripts were submitted on CUBIC as jobs via the following wrapper script: 
[/manuscript/code/connectivity_matrices/wrapper_connectivity_matrices.sh](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/wrapper_connectivity_matrices.sh)

This wrapper script calls a separate script that runs a singularity image containing all the required R packages. The wrapper above submits this script as a job on CUBIC. You can check out the scripts running singularity images at `/manuscript/code/connectivity_matrices/<dataset>/<dataset>_makeConnMatrices.sh`.
 


**3b. Sensitivity analysis:** 

We wanted to investigate whether our findings were potentially driven by confounding factors including use of task and rest scans and atlas used for cortical parcellation. 

First, sensitivity analyses were performed with only resting state data while excluding fMRI acquired during task conditions. Datasets used in this rest-only sensitivity analysis include PNC, HCP-D, and HBN. Analyses for NKI were completed with resting-state fMRI only due to absence of task scans. We included participants with at least 6 minutes of resting-state fMRI. We analyzed data from 998 participants (549 females) from PNC, 611 participants (328 females) from HCP-D, and 842 participants (1477 females) from HBN. The total scan time for fully acquired resting-state scans was 11 minutes and 12 seconds (224 volumes) for PNC, 25 minutes and 30 seconds (1912 volumes) for HCP-D, and 10 minutes and 9 seconds (750 volumes) for HBN. 

Second, analyses were also evaluated using additional cortical parcellations. Our primary parcellation utilized the Schaefer 200 atlas; secondary atlases included the Schaefer 400 atlas, the Gordon atlas, and HCP-MMP atlas. For analyses of secondary outcome measures that require community structure, namely average between- and within-network connectivity, we evaluated both the Yeo 7 and 17-network partitions associated with the Schaefer atlas. Results from sensitivity analyses are presented in the supplement.

Lastly, we made connectivity matrices for PNC using absolute correlation coefficient, thresholded matrices, and without use of GSR (from CPAC). 

*Data*: resting-state fMRI only was used to construct connectivity matrices with the following code: 

* PNC: [/manuscript/code/connectivity_matrices/PNC/PNC_makeConnMatrices_restOnly.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/PNC/PNC_makeConnMatrices_restOnly.R)
* NKI: no sensitivity analyses done with resting-state only since NKI only has resting-state to begin with! 
* HCP-D: [/manuscript/code/connectivity_matrices/HCPD/HCPD_makeConnMatrices_restOnly.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/HCPD/HCPD_makeConnMatrices_restOnly.R)
* HBN:  [/manuscript/code/connectivity_matrices/HBN/HBN_makeConnMatrices_restOnly.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/HBN/HBN_makeConnMatrices_restOnly.R)


Furthermore, absolute correlation coefficient, thresholded matrices, and non-GSR matrices in PNC were made using the following: 
* [/manuscript/code/scripts/makeConnMatrices_thresholded_absvalue.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/makeConnMatrices_thresholded_absvalue.R) 
* [/manuscript/code/connectivity_matrices/PNC/PNC_makeConnMatrices_nonGSR.R](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/PNC/PNC_makeConnMatrices_nonGSR.R)

Sensitivity analyses were submitted using the following two wrapper scripts:
* [/manuscript/code/connectivity_matrices/wrapper_connectivity_matrices_sens_analyses.sh](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/wrapper_connectivity_matrices_sens_analyses.sh): creates rest only matrices for relevant datasets and non-GSR for PNC  
* [/manuscript/code/connectivity_matrices/wrapper_PNC_sens_connectivity_matrices.sh](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_matrices/wrapper_PNC_sens_connectivity_matrices.sh): makes abs value and thresholded matrices in PNC. This script was submitted after the main analysis connectivity matrices were made, since it depends on those to generate absolute correlation coefficient and thresholded matrices.

 


### 5. Quantifying functional connectivity metrics: global brain connectivity (aka FC strength), between- and within-network connectivity  
Global brain connectivity (GBC), which we call functional connectivity (FC) strength in the manuscript, was calculated for each cortical parcel by averaging its timeseries correlation with all other parcels. Hence, global brain connectivity represents the mean edge strength of a given region with all other regions, without thresholding. Average between-network connectivity (BNC) was defined as the mean edge strength (Pearson correlation) of a given region and all other regions not in that region’s network. Average within-network connectivity (WNC) was defined as the mean edge strength (Pearson correlation) of a given region and all other regions within that region’s network. We also examined functional connectivity at the edge level by extracting the Pearson correlation between timeseries for each pair of regions. 

**5a. Main analyses:**  
[/manuscript/code/scripts/main_analyses/computeConnectivityMetrics.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/computeConnectivityMetrics.R): contains the functions needed for computing the conn measures

GBC, BNC, WNC, and edge-level connectivity are computed using the following scripts:
* [/manuscript/code/scripts/main_analyses/computeGBC.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/computeGBC.R)
* [/manuscript/code/scripts/main_analyses/computeBNC.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/computeBNC.R)
* [/manuscript/code/scripts/main_analyses/computeWNC.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/computeWNC.R)
* [/manuscript/code/scripts/main_analyses/computeEdge.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/computeEdge.R)

 
The above Rscripts were submitted on CUBIC via the following wrapper script: 
[/manuscript/code/connectivity_measures/wrapper_connmetrics_main_analysis.sh](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_measures/wrapper_connmetrics_main_analysis.sh)

This wrapper script calls a separate, dataset-specific script that runs a singularity image containing all the required R packages. The wrapper above submits this script as a job on CUBIC. You can look through the scripts running singularity images at `/manuscript/code/connectivity_measures/<dataset>/compute<metric>_<dataset>.sh`.

**5b. Sensitivity analyses:**  
[/manuscript/code/scripts/sensitivity_analyses/computeConnectivityMetrics_restOnly.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/computeConnectivityMetrics_restOnly.R): contains the functions needed for computing the conn measures using rest-only data

The following describes the scripts used for computing each metric:
* GBC based on rest-only data was computed using [/manuscript/code/scripts/sensitivity_analyses/computeGBC_restOnly.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/computeGBC_restOnly.R)

* GBC based on absolute correlation coefficient and thresholded matrices in the PNC were computed using [/manuscript/code/scripts/sensitivity_analyses/computeGBC_absvalue_thresholded.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/computeGBC_absvalue_thresholded.R)

* GBC based on data that was not pre-processed using GSR in the PNC were computed using [/manuscript/code/scripts/sensitivity_analyses/computeGBC_nonGSR.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/computeGBC_nonGSR.R) 

* Connectivity between pairs of networks was computed using [/manuscript/code/scripts/sensitivity_analyses/computeNetworkPair_conn.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/computeNetworkPair_conn.R)

 
The above Rscripts were submitted on CUBIC via the following wrapper script: 
[/manuscript/code/connectivity_measures/wrapper_connmetrics_sens_analysis.sh](https://github.com/PennLINC/network_replication/blob/main/code/connectivity_measures/wrapper_connmetrics_sens_analysis.sh)

This wrapper script calls a separate, dataset-specific script that runs a singularity image containing all the required R packages. The wrapper above submits this script as a job on CUBIC. You can look through the scripts running singularity images at `/manuscript/code/connectivity_measures/<dataset>/<sensitivity_analysis_name>.sh`

### 6. Image harmonization: applying [covbat-gam](https://github.com/andy1764/ComBatFamily) to multi-site data (HCP-D and HBN)

[Correcting Covariance Batch Effects (CovBat)](https://onlinelibrary.wiley.com/doi/10.1002/hbm.25688) was used to harmonize multi-site MRI data to ensure the imaging measures were comparable across sites. CovBat was applied to functional connectivity metrics for multi-site dataset (HCP-D and HBN) using the [ComBatFamily R package](https://github.com/andy1764/ComBatFamily). Sex and in-scanner motion were included as covariates with age modeled as a smooth term via a [generalized additive model](https://www.sciencedirect.com/science/article/pii/S1053811919310419?via%3Dihub) using the ‘covfam’ function, which enables flexible covariate modeling with penalized splines.  

Covbat was applied to GBC, BNC, and WNC in the following Rmd files:
 
+ HCP-D: [/manuscript/code/covbat_harmonization/HCPD/HCPD_covbat_FCmetrics.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/covbat_harmonization/HCPD/HCPD_covbat_FCmetrics.Rmd)
+ HBN: [/manuscript/code/covbat_harmonization/HBN/HBN_covbat_FCmetrics.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/covbat_harmonization/HBN/HBN_covbat_FCmetrics.Rmd)

Harmonization of edge-level data was completed using:
[/manuscript/code/scripts/main_analyses/covbat_edges.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/covbat_edges.R)

* Jobs for HCP-D and HBN were submitted using [/manuscript/code/covbat_harmonization/HCPD/covbat_edges_HCPD.sh](https://github.com/PennLINC/network_replication/blob/main/code/covbat_harmonization/HBN/covbat_edges_HBN.sh) and [/manuscript/code/covbat_harmonization/HBN/covbat_edges_HBN.sh](https://github.com/PennLINC/network_replication/blob/main/code/covbat_harmonization/HBN/covbat_edges_HBN.sh) respectively


Covbat was applied to metrics used in sensitivity analyses using the following Rmd files:
 
+ HCP-D: [/manuscript/code/covbat_harmonization/HCPD/HCPD_covbat_sensitivity.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/covbat_harmonization/HCPD/HCPD_covbat_sensitivity.Rmd)
+ HBN: [/manuscript/code/covbat_harmonization/HBN/HBN_covbat_sensitivity.Rmd](https://github.com/PennLINC/network_replication/blob/main/code/covbat_harmonization/HBN/HBN_covbat_sensitivity.Rmd)


### 8. Fitting generalized additive models (GAMs) 

**8a. Main analyses:**  
GAMs were fit for GBC, BNC, WNC, and edge-level connectivity for each cortical region to examine age-dependent changes in each of these functional connectivity metrics. 

[/manuscript/code/scripts/main_analyses/GAM_functions.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/GAM_functions.R) includes the set of functions to fit GAM models. This script includes:  
* *gam.fit*: Function to fit GAM `measure ~ s(smooth_var) + covariates)` per each region in atlas and save out statistics and derivative-based characteristics
* *gam.predsmooth*: Function to fit GAM smooths based on model-predicted data
* *gam.smooth.predict_posterior*: Function to predict fitted values of a measure based on a fitted GAM smooth `measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)` and a prediction df and for individual draws from the simulated posterior distribution

<br>

The following describes the scripts used for fitting GAMs for each main analysis:

* To fit GAMs for FC metrics, estimate GAM smooths, predict fitted GBC values, and compute alignment of fitted GBC values with the S-A axis across age, we used [/manuscript/code/scripts/main_analyses/fitGAMs_FCmetrics.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/fitGAMs_FCmetrics.R). 


> Note: the age-resolved analysis (predicting fitted GBC values and computing alignment of values with S-A axis across age) was done to examine how the spatial distribution of FC strength aligns with the S-A axis across the broad age range studied. This analysis was performed to gain insight into whether the spatial patterning of functional connectivity across the cortical mantle becomes increasingly hierarchical through development. 

* To fit GAMs for edges, we used [/manuscript/code/scripts/main_analyses/fitGAMs_edges.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/main_analyses/fitGAMs_edges.R). 


These two scripts for fitting GAMs were run for all datasets using [/manuscript/code/fit_GAMs/wrapper_fitGAMs_main_analysis.sh](https://github.com/PennLINC/network_replication/blob/main/code/fit_GAMs/wrapper_fitGAMs_main_analysis.sh). 
 
This wrapper script calls separate scripts that runs a singularity image containing all the required R packages. The wrapper submits these scripts as jobs on CUBIC. You can look through the scripts running singularity images at `/manuscript/code/scripts/main_analyses/fitGAMs_<main_analysis_name>_singularity.sh`.


**8b. Sensitivity analyses:**  

The following describes the scripts used for fitting GAMs for each sensitivity analysis:
* GAMs for GBC based on rest-only data were fit using [/manuscript/code/scripts/sensitivity_analyses/fitGAMs_GBC_restOnly.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/fitGAMs_GBC_restOnly.R)

* GAMs for GBC based on absolute correlation coefficient and thresholded matrices in the PNC were computed using [/manuscript/code/scripts/sensitivity_analyses/fitGAMs_GBC_absvalue_thresholded.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/fitGAMs_GBC_absvalue_thresholded.R)

* GAMs for GBC based on data that was not pre-processed using GSR in the PNC were computed using [/manuscript/code/scripts/sensitivity_analyses/fitGAMs_GBC_nonGSR.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/fitGAMs_GBC_nonGSR.R) 

* GAMs for connectivity between pairs of networks were computed using [/manuscript/code/scripts/sensitivity_analyses/fitGAMs_networkpair.R](https://github.com/PennLINC/network_replication/blob/main/code/scripts/sensitivity_analyses/fitGAMs_networkpair.R)

 
The above Rscripts were submitted on CUBIC via the following wrapper script: 
[/manuscript/code/connectivity_measures/wrapper_fitGAMs_sens_analysis.sh](https://github.com/PennLINC/network_replication/blob/main/code/fit_GAMs/wrapper_fitGAMs_sens_analysis.sh)

This wrapper script calls separate scripts that runs a singularity image containing all the required R packages. The wrapper submits these scripts as jobs on CUBIC. You can look through the scripts running singularity images at `/manuscript/code/scripts/sensitivity_analyses/fitGAMs_<sensitivity_analysis_name>_singularity.sh`.


### 9. Characterization of relationships between functional connectivity metrics, age, and the S-A axis

We used Spearman’s rank correlations to quantify the association between S-A axis ranks and observed developmental effects. This was performed within 
[/manuscript/results/main_figures.Rmd](https://github.com/PennLINC/network_replication/blob/main/results/main_figures.Rmd).

To investigate how the development of edge-level connectivity differs across sensorimotor-association axis, we examined age-related changes in connectivity across edges by fitting a bivariate smooth interaction. The effect of S-A axis rank on edge-level age effects was modeled using a tensor product smooth. This analysis can be found in [lines 885-919](https://github.com/PennLINC/network_replication/blob/main/results/main_figures.Rmd#L885-L919) from the Rmd above. 


<br>

### 10. Visualizing the results
Results from main analyses are visualized using this Rmd file: [/manuscript/results/main_figures.Rmd](https://github.com/PennLINC/network_replication/blob/main/results/main_figures.Rmd). The knitted Rmd file displaying main figures 1-6 can be downloaded at [/manuscript/results/main_figures.html](https://github.com/PennLINC/network_replication/blob/main/results/main_figures.html) and viewed on your browswer.

Results from sensitivity analyses are visualized using this Rmd file: [/manuscript/results/supp_figures.Rmd](https://github.com/PennLINC/network_replication/blob/main/results/supp_figures.Rmd). The knitted Rmd file displaying supplementary figures can be downloaded at [/manuscript/results/supp_figures.html](https://github.com/PennLINC/network_replication/blob/main/results/supp_figures.html) and viewed on your browswer.
