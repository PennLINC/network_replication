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
Valerie J. Sydnor, Arielle S. Keller, Aaron F. Alexander-Bloch, Matthew Cieslak, Sydney Covitz, Andrew Chen, Eric Feczko, Alexandre R. Franco, Raquel E. Gur, Ruben C. Gur, Audrey Houghton, Fengling Hu, Gregory Kiar, Bart Larsen, Adam Pinesa, Giovanni Salumg, Tinashe Taperaa, Ting Xu, Damien A. Fair, Michael P. Milham, Theodore D. Satterthwaite

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
/cbica/projects/network_replication

```
subdirectory/: description
```



<br>
<br>
# CODE DOCUMENTATION  

All project analyses are described below along with the corresponding code on Github. The analytic workflow was executed in the following order:  

1. Parcellating the sensorimotor-association (S-A) axis  
2. Sample selection for each dataset: PNC (discovery), NKI, HCP-D, and HBN (replication)  
3. Constructing connectivity matrices for each dataset
* *Main analysis:* 
    + Data: concatenated task and resting-state fMRI
    + Cortical parcellation: [Schaefer 200](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/Schaefer2018_200Parcels_17Networks_order.dlabel.nii)
    + Network sollution: [7 Network](https://github.com/ThomasYeoLab/CBIG/blob/6d1400a2d643261246f6b042e7ef5fbe417506cd/utilities/matlab/FC/CBIG_ReorderParcelIndex.m)
* *Sensitivity analysis:* resting-state fMRI only 
4. Image harmonization: applying [covbat-gam](https://github.com/andy1764/ComBatFamily) to multi-site data (HCP-D and HBN)
5. Quantification of global brain connectivity, between- and within-network connectivity
6. Fitting generalized additive models (GAMs) 
7. Characterization of relationships between functional connectivity metrics, age, and the S-A axis

his workflow includes quantification of regional fluctuation amplitude, PNC sample selection, fitting of generalized additive models (GAMs), and characterization of relationships between fluctuation amplitude, age, environmental variability, and the sensorimotor-association axis. Scripts were implemented in the order outlined below.
<br>
### Fluctuation Amplitude Quantification
Resting state functional MRI data were processed with [fmriprep 20.2.3](https://hub.docker.com/layers/fmriprep/nipreps/fmriprep/20.2.3/images/sha256-102db5fe8b0a34298f2eb2fd5962ad99ff0a948d258cbf08736fcc1b845cad9f?context=explore) and [xcp-d 0.0.4](https://hub.docker.com/layers/xcp_abcd/pennlinc/xcp_abcd/0.0.4/images/sha256-317160b8078cf7978eaf9db6fef32df78864232cb8a8759a354832813d1faf02?context=explore) to quantify fluctuation amplitude at each vertex on the fslr 32k cortical surface. 

fmriprep 20.2.3 was run with the following parameters:

```bash
$ singularity run pennlinc-containers/.datalad/environments/fmriprep-20-2-3/image inputs/data prep participant --output-spaces MNI152NLin6Asym:res-2 --participant-label "$subid" --force-bbr --cifti-output 91k -v -v
```

xcp-d  0.0.4 was run with the following parameters: 

```bash
$ singularity run pennlinc-containers/.datalad/environments/xcp-abcd-0-0-4/image inputs/data/fmriprep xcp participant --despike --lower-bpf 0.01 --upper-bpf 0.08 --participant_label $subid -p 36P -f 10 –cifti
```

Vertex-wise fluctuation amplitude maps were then parcellated with [/fluctuation_amplitude/parcellate.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/fluctuation_amplitude/parcellate.Rmd) to quantify mean fluctuation amplitude in each cortical region. Regions were defined with the [HCP multimodal atlas](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/glasser_space-fsLR_den-32k_desc-atlas.dlabel.nii) (i.e. Glasser 360, primary analyses) and with the [Schaefer 400 atlas](https://github.com/PennLINC/xcp_d/blob/main/xcp_d/data/ciftiatlas/Schaefer2018_400Parcels_17Networks_order.dlabel.nii) (sensitivity analysis). 

Fluctuation amplitude analyses were only conducted in brain regions that reliably exhibited high signal to noise ratio (SNR) in PNC functional MRI data. The vertex level SNR map generated in Cui et al., 2020, Neuron was parcellated with Glasser 360 and Schaefer 400 atlases with the script [/fluctuation_amplitude/SNR_mask.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/fluctuation_amplitude/SNR_mask.R) for use in study analyses. Regions wherein >= 25% of vertices had attenuated signal were excluded from analyses.


### Sample Construction
The final study sample was constructed with [/sample_construction/finalsample.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/sample_construction/finalsample.Rmd). The final sample was generated from the 1374 PNC participants with dominant group ses-PNC1_task-rest_acq-singleband scans (non-variant CuBIDS acquisitions). The following exclusions were applied to generate the final sample of N = 1033 participants:

> *Health exclude*: 120 participants with medical problems that could impact brain function or incidentally-encountered brain structure abnormalities were excluded from the sample  
> *T1 QA exclude*: 23 participants with T1-weighted scans that failed visual QC were excluded from the sample  
> *fMRI motion exclude*: 179 participants with a mean relative RMS motion value > 0.2 during the resting state fMRI scan were excluded from the sample  
> *Fluctuation amplitude outlier exclude*: from the remaining 1052 participants, 19 individuals that had outlier (+- 4 SD from the mean) fluctuation amplitude data in more than 5% of Glasser 360 parcels were excluded from the sample   
<br>

### Model Fitting 
***GAM Functions***  
To characterize age-dependent changes in spontaneous activity fluctuations across the developing cortex as well as associations between fluctuation amplitude and environmental factors, generalized additive models were fit in each cortical region. GAM models and associated statistics, fitted values, smooths, and derivatives were quantified with the set of functions included in [/gam_models/GAM_functions.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/GAM_functions.R). This script includes the following functions:
-	*gam.fit.smooth*: A function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and save out statistics (smooth partial R squared and p-value) and derivative-based characteristics
-	*gam.smooth.predict*: A function to predict fitted values of a measure based on a GAM model and a prediction data frame
-	*gam.estimate.smooth*: A function to estimate zero-averaged gam smooth functions
-	*gam.posterior.smooths*: A function to simulate the posterior distribution from a fit GAM model, calculate smooths for individual posterior draws, and determine smooth max and min median values + 95% credible intervals
-	*gam.derivatives*: A function to compute derivatives for the smooth term from a main GAM model and for individual draws from the simulated posterior distribution; can return true model derivatives or posterior derivatives
-	*gam.fit.covariate*: A function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariate of interest + control covariates)) and save out statistics (partial R squared and p-value) for covariate of interest
-	*gam.varyingcoefficients*: A function to estimate how the linear association between a predictor and y varies along a smooth function (varying coefficient interaction GAMs); can calculate associations (varying coefficient slopes) along the smooth function for the true model or for individual draws from the simulated posterior distribution
-	*gam.factorsmooth.interaction*: A function to fit a GAM with a factor-smooth interaction and obtain statistics for the interaction term 

***Model Fitting: Age-Dependent Changes in Regional Fluctuation Amplitude***  
Developmental models were fit with [gam_models/fitGAMs_fluctuationamplitude_age.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/fitGAMs_fluctuationamplitude_age.R), which calls the functions described above. Age-focused GAMs were implemented for all regions included in the Glasser 360 and Schaefer 400 atlases, using the final study sample of N = 1033 generated during the sample construction process. Models examining the impact of participant pubertal stage on regional fluctuation amplitude were fit with [gam_model/fitGAMs_fluctuationamplitude_pubertalstage.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/fitGAMs_fluctuationamplitude_pubertalstage.R).

***Model Fitting: Developmental Environment-Dependent Variation in Regional Fluctuation Amplitude***  
Models examining age-independent (main effect) and age-dependent (interaction) associations between fluctuation amplitude and neighborhood socioeconomic circumstances were run via [gam_models/fitGAMs_fluctuationamplitude_environment.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/fitGAMs_fluctuationamplitude_environment.R) on the N = 1033 study sample. 


### Data Interpretation and Visualization 

Model results were examined and studied within our hierarchical neurodevelopmental plasticity framework in the markdown file [/developmental_effects/hierarchical_development.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd). This markdown executes the following:

- Generates all manuscript Figures
- Examines the cortical distribution and significance of associations between regional fluctuation amplitude and age
- Quantifies associations between developmental changes in fluctuation amplitude and in the T1w/T2w ratio
- Characterizes alignment of fluctuation amplitude age-related changes with the sensorimotor-association axis, including the magnitude and timing of change
- Calculates age smooths for 10 bins of the sensorimotor-association axis
- Performs a temporal sliding window analysis to uncover when developmental change in fluctuation amplitude is maximally aligned (and not aligned) with the sensorimotor-association axis
- Examines the significance and distribution of associations between regional fluctuation amplitude and environmental variability, and how environment effects are stratified by the S-A axis

A rendered html of hierarchical_development.Rmd can be viewed [here](https://htmlpreview.github.io/?https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.html)!

### Sensitivity Analyses

The robustness of our developmental findings was confirmed in a series of sensitivity analyses. Sensitiviy analysis GAMs were fit with [sensitivity_analyses/fitGAMs_sensitivityanalyses.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/sensitivity_analyses/fitGAMs_sensitivityanalyses.R) and results were examined in [sensitivity_analyses/sensitivityresults.Rmd](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/sensitivity_analyses/sensitivityresults.Rmd). The following sensitivity analyses were performed:

- *Low motion*: Findings were assessed in a low motion sample of N = 690 participants with a relative mean RMS < 0.075
- *Psychiatry exclusions*: Findings were assessed in a sample that excluded participants with current psychotropic medication use or a history of psychiatric hospitalization
- *Vascular control*: All fluctuation amplitude models were re-fit while controlling for regional cerebral blood flow (quantified with ASL)
- *T2 signal control*: All fluctuation amplitude models were re-fit while controlling for regional BOLD signal level (T2*) during the fMRI scan as quantified with [sensitivity_analyses/parcellated_meanT2star.R](https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/sensitivity_analyses/parcellated_meanT2star.R)
- *Normalization*: Developmental GAMs were run with mean normalized fluctuation amplitude as the dependent variable
- *Atlas*: Results were reproduced using the Schaefer 400 atlas