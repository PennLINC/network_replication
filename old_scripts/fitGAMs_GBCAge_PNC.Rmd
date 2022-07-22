---
title: "GBC GAMs: PNC"
author: "Valerie Jill Sydnor"
output:
  html_document:
    code_folding: show
    highlight: haddock
    theme: lumen
    toc: yes
    toc_depth: 4
    toc_float: yes
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
source("./GAM_functions.R")
```


**Prepare GBC + Demographics data frame for GAMs**

```{r}
participants <- read.csv("") #final study sample pnc IDs (bblid/scanid) from Adam Pines's multiscale paper
demographics <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/sample_info/PNC/n1601_demographics_go1_20161212.csv") #PNC demographics information
demographics <- demographics %>% mutate(age = (ageAtScan1/12)) #format age variable
demographics$sex <- as.factor(demographics$sex) #format sex variable

demographics <- merge(participants, demographics, by="scanid") #get demographics for final study sample only via ID matching
relRMS <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/RBC_timeseries/PNC/qc_files/RBC_PNC_relMeans.csv") 
demographics <- merge(demographics, relRMS, by="rbcid") #add motion information to demographics df
```

```{r}
gbc.subxparcel.matrix.glasser <- read.csv("GBC_subxparcel_matrix_glasser.csv") #participant x glasser region GBC spreadsheet
gbc.glasser.pnc <- merge(gbc.subxparcel.matrix.glasser, demographics, by="rbcid") #get glasser GBC data and demographics for final study sample via ID matching
write.csv(gbc.glasser.pnc, "GBCglasser_demographics_finalsample.csv", quote = F, row.names = F)

gbc.subxparcel.matrix.schaefer <- read.csv("GBC_subxparcel_matrix_schaefer.csv") #participant x schaefer region GBC spreadsheet
gbc.schaefer.pnc <- merge(gbc.subxparcel.matrix.schaefer, demographics, by="rbcid") #get schaefer GBC data and demographics for final study sample via ID matching
write.csv(gbc.schaefer.pnc, "GBCschaefer_demographics_finalsample.csv")
```

```{r}
glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T) #glasser parcel names in order of surface data
schaefer.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/schaefer400_regionlist.csv", header = T) #schaefer parcel names in order of surface data
```

**Fit GAMs**

Fit GAM (GBC ~ s(age) + sex + RMSmotion)) per each region in atlas and save out statistics and derivative-based characteristics

```{r}
gam.GBC.age.glasser <- matrix(data=NA, nrow=360, ncol=9) #empty matrix to save gam.fit output to

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  region <- glasser.parcel.labels$orig_parcelname[row] 
  GAM.RESULTS <- gam.fit(measure = "gbc", atlas = "glasser", dataset = "pnc", region = region, smooth_var = "age", covariates = "sex + RMSmotion") #run the gam.fit function
  gam.GBC.age.glasser[row,] <- GAM.RESULTS}

gam.GBC.age.glasser <- as.data.frame(gam.GBC.age.glasser)
colnames(gam.GBC.age.glasser) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.GBC.increase","age.maturation")
cols = c(2:9)    
gam.GBC.age.glasser[,cols] = apply(gam.GBC.age.glasser[,cols], 2, function(x) as.numeric(as.character(x)))
gam.GBC.age.glasser <- gam.GBC.age.glasser %>% mutate(significant = (Anova.age.pvalue < 0.05))

write.csv(gam.GBC.age.glasser, "GAMresults.GBC.age.glasser.csv", row.names = F, quote = F)
```

```{r}
gam.GBC.age.schaefer <- matrix(data=NA, nrow=400, ncol=9) #empty matrix to save gam.fit output to

for(row in c(1:nrow(schaefer.parcel.labels))){ #for each schaefer region
  region <- schaefer.parcel.labels$label[row] 
  GAM.RESULTS <- gam.fit(measure = "gbc", atlas = "schaefer", dataset = "pnc", region = region, smooth_var = "age", covariates = "sex + RMSmotion") #run the gam.fit function
  gam.GBC.age.schaefer[row,] <- GAM.RESULTS}

gam.GBC.age.schaefer <- as.data.frame(gam.GBC.age.schaefer)
colnames(gam.GBC.age.schaefer) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.GBC.increase","age.maturation")
cols = c(2:9)    
gam.GBC.age.schaefer[,cols] = apply(gam.GBC.age.schaefer[,cols], 2, function(x) as.numeric(as.character(x)))
gam.GBC.age.schaefer <- gam.GBC.age.schaefer %>% mutate(significant = (Anova.age.pvalue < 0.05))

write.csv(gam.GBC.age.schaefer, "GAMresults.GBC.age.schaefer.csv", row.names = F, quote = F)
```

Estimate GAM smooths based on model-predicted data and save out predicted y data

```{r}
gam.smooths.glasser <- matrix(data=NA, ncol=7) #empty matrix to save gam.predsmooth fits to
colnames(gam.smooths.glasser) <- c("age","fit","se.fit","selo","sehi","index","label")

gam.peaks.glasser <- matrix(data=NA, ncol=2) #empty matrix to save measure peaks to (value of x when y is predicted to be maximal)
colnames(gam.peaks.glasser) <- c("label","age.peak")

for(row in c(1:nrow(glasser.parcel.labels))){ #for each glasser region
  
  region <- glasser.parcel.labels$orig_parcelname[row] #get the region name
  GAM.SMOOTH <- gam.predsmooth(measure = "gbc", atlas = "glasser", dataset = "pnc", region = region, smooth_var = "age", covariates = "sex + RMSmotion") #run the gam.predsmooth function
  
  preddata <- as.data.frame(GAM.SMOOTH[3]) #get predicted.smooth df from function output
  preddata$index <- rep(x=row, 1000) #region index
  preddata$label <- rep(x=GAM.SMOOTH[1], 1000) #label
  gam.smooths.glasser <- rbind(gam.smooths.glasser, preddata)
  
  datapeak <- as.data.frame(cbind(GAM.SMOOTH[1], GAM.SMOOTH[2]))
  colnames(datapeak) <- c("label","age.peak")
  gam.peaks.glasser <- rbind(gam.peaks.glasser, datapeak)
}
gam.smooths.glasser <- gam.smooths.glasser[-1,] #remove empty initilization row
gam.smooths.glasser$label <- as.character(gam.smooths.glasser$label)
gam.peaks.glasser <- gam.peaks.glasser[-1,] #remove empty initilization row
gam.peaks.glasser$label <- as.character(gam.peaks.glasser$label)
gam.peaks.glasser$age.peak <- as.numeric(gam.peaks.glasser$age.peak)

write.csv(gam.smooths.glasser,"GAMsmoothfits.GBC.age.glasser.csv", row.names = F, quote = F)
write.csv(gam.peaks.glasser, "GAMpeaks.GBC.age.glasser.csv", row.names = F, quote = F)
```

```{r}
gam.smooths.schaefer <- matrix(data=NA, ncol=7) #empty matrix to save gam.predsmooth fits to
colnames(gam.smooths.schaefer) <- c("age","fit","se.fit","selo","sehi","index","label")

gam.peaks.schaefer <- matrix(data=NA, ncol=2) #empty matrix to save measure peaks to (value of x when y is predicted to be maximal)
colnames(gam.peaks.schaefer) <- c("label","age.peak")

for(row in c(1:nrow(schaefer.parcel.labels))){ #for each schaefer region
  
  region <- schaefer.parcel.labels$label[row] #get the region name
  GAM.SMOOTH <- gam.predsmooth(measure = "gbc", atlas = "schaefer", dataset = "pnc", region = region, smooth_var = "age", covariates = "sex + RMSmotion") #run the gam.predsmooth function
  
  preddata <- as.data.frame(GAM.SMOOTH[3]) #get predicted.smooth df from function output
  preddata$index <- rep(x=row, 1000) #region index
  preddata$label <- rep(x=GAM.SMOOTH[1], 1000) #label
  gam.smooths.schaefer <- rbind(gam.smooths.schaefer, preddata)
  
  datapeak <- as.data.frame(cbind(GAM.SMOOTH[1], GAM.SMOOTH[2]))
  colnames(datapeak) <- c("label","age.peak")
  gam.peaks.schaefer <- rbind(gam.peaks.schaefer, datapeak)
}
gam.smooths.schaefer <- gam.smooths.schaefer[-1,] #remove empty initilization row
gam.smooths.schaefer$label <- as.character(gam.smooths.schaefer$label)
gam.peaks.schaefer <- gam.peaks.schaefer[-1,] #remove empty initilization row
gam.peaks.schaefer$label <- as.character(gam.peaks.schaefer$label)
gam.peaks.schaefer$age.peak <- as.numeric(gam.peaks.schaefer$age.peak)

write.csv(gam.smooths.schaefer,"GAMsmoothfits.GBC.age.schaefer.csv", row.names = F, quote = F)
write.csv(gam.peaks.schaefer, "GAMpeaks.GBC.age.schaefer.csv", row.names = F, quote = F)
```
