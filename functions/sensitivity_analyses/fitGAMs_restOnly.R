# functions for fitting GAMs 
library(dplyr)
library(janitor)
library(purrr)
source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/GAM_functions.R") 

# load parcellated S-A axis 

schaefer200x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x17_SAaxis.csv")
schaefer200x17_SAaxis$label <- gsub("17Network", "Network", schaefer200x17_SAaxis$label) 
schaefer200_SAaxis <- schaefer200x17_SAaxis 
 

# load parcel labels   
schaefer200x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
schaefer200.parcel.labels <- schaefer200x17.parcel.labels #gbc
 
# Function for making connMetrics_demographics.csv
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs  
# @param dataset A character string, name of dataset
make_connMetricsDemog <- function(atlas_network, metric, demog_csv, subj_id, dataset){
  subxparcel.matrix <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/%2$s_subxparcel_matrix_%3$s.csv", dataset, metric, atlas_network))
  print(paste(atlas_network, metric, "file loaded"))
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  write.csv(connMetricsDemog_df, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/%2$s%3$s_demographics_finalsample.csv", dataset, metric, atlas_network), row.names = F)
  return(connMetricsDemog_df)
}


# Function for making connMetrics_demographics.csv with covbat
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_connMetricsDemog_covbat <- function(atlas_network, metric, demog_csv, subj_id, dataset){
  subxparcel.matrix <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/%2$s_subxparcel_matrix_%3$s.csv", dataset, metric, atlas_network))
  print(paste(atlas_network, metric, "file loaded"))
  subxparcel.matrix$subject <- demographics$subject
  subxparcel.matrix <- subxparcel.matrix %>% dplyr::relocate(subject)
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  write.csv(connMetricsDemog_df, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/%2$s%3$s_demographics_finalsample_covbat.csv", dataset, metric, atlas_network), row.names = F)
  return(connMetricsDemog_df)
}
 


# Function for **Fitting GAMs** for GBC, BNC, and WNC
# Fit GAM (func_conn_metric ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name) {
  parcel.labels <- parcel.labels$label
  gam.age <- matrix(data=NA, nrow=length(parcel.labels), ncol=9) #empty matrix to save gam.fit output to
  
  for(row in c(1:length(parcel.labels))){ #for each region
    region <- parcel.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name, network_parcellation), dataset = dataset_name, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  write.csv(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  
  return(gam.age)
}


# Function for **Fitting GAMs** for GBC, BNC, and WNC - covbat
# Fit GAM (func_conn_metric ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs_covbat <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name) {
  parcel.labels <- parcel.labels$label
  gam.age <- matrix(data=NA, nrow=length(parcel.labels), ncol=9) #empty matrix to save gam.fit output to
  
  for(row in c(1:length(parcel.labels))){ #for each region
    region <- parcel.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name, network_parcellation), dataset = dataset_name, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  write.csv(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s_covbat.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  
  return(gam.age)
}

 
 
# Function for estimating GAM smooths based on model-predicted data and save out predicted y data
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
estimate_GAMsmooths <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name) {
  parcel.labels <- parcel.labels$label
  gam.smooths <- matrix(data=NA, ncol=7) #empty matrix to save gam.predsmooth fits to
  colnames(gam.smooths) <- c("age","fit","se.fit","selo","sehi","index","label")
  
  for(row in c(1:length(parcel.labels))){ #for each region
    
    region <- parcel.labels[row] #get the region name
    GAM.SMOOTH <- gam.predsmooth(measure = metric, atlas = paste0(atlas_name, network_parcellation), dataset = dataset_name, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.predsmooth function
    
    preddata <- as.data.frame(GAM.SMOOTH[3]) #get predicted.smooth df from function output
    preddata$index <- rep(x=row, 1000) #region index
    preddata$label <- rep(x=GAM.SMOOTH[1], 1000) #label
    gam.smooths <- rbind(gam.smooths, preddata)
    
    
  }
  gam.smooths <- gam.smooths[-1,] #remove empty initilization row
  gam.smooths$label <- as.character(gam.smooths$label)
  
  
  write.csv(gam.smooths, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s%4$s.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)

  return(gam.smooths)
  
}
 

# Function for estimating GAM smooths based on model-predicted data and save out predicted y data - covbat
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
estimate_GAMsmooths_covbat <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name) {
  parcel.labels <- parcel.labels$label
  gam.smooths <- matrix(data=NA, ncol=7) #empty matrix to save gam.predsmooth fits to
  colnames(gam.smooths) <- c("age","fit","se.fit","selo","sehi","index","label")
  
  for(row in c(1:length(parcel.labels))){ #for each region
    
    region <- parcel.labels[row] #get the region name
    GAM.SMOOTH <- gam.predsmooth(measure = metric, atlas = paste0(atlas_name, network_parcellation), dataset = dataset_name, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.predsmooth function
    
    preddata <- as.data.frame(GAM.SMOOTH[3]) #get predicted.smooth df from function output
    preddata$index <- rep(x=row, 1000) #region index
    preddata$label <- rep(x=GAM.SMOOTH[1], 1000) #label
    gam.smooths <- rbind(gam.smooths, preddata)
    
    
  }
  gam.smooths <- gam.smooths[-1,] #remove empty initilization row
  gam.smooths$label <- as.character(gam.smooths$label)
  
  
  write.csv(gam.smooths, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s%4$s_covbat.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  
  return(gam.smooths)
}

 


