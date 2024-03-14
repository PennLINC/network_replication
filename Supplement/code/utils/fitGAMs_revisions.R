# functions for fitting GAMs 
library(dplyr)
library(janitor)
library(purrr)
source("/cbica/projects/network_replication/revisions/code/utils/GAM_functions_revisions.R")

# load parcellated S-A axis
schaefer200x7_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x7_SAaxis.csv")
schaefer200x7_SAaxis$label <- gsub("7Network", "Network", schaefer200x7_SAaxis$label) 
 

# load parcel labels  
schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
 


# Function for making connMetrics_demographics dataframe
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_connMetricsDemog <- function(atlas_network, metric, demog_csv, subj_id, dataset, conn_type){
  subxparcel.matrix <- read.csv(sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/%2$s_subxparcel_matrix_%3$s_%4$s.csv", dataset, metric, atlas_network, conn_type))
  print(paste(atlas_network, metric, "file loaded"))
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  write.csv(connMetricsDemog_df, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample_%4$s.csv", dataset, metric, atlas_network, conn_type), row.names = F)
  return(connMetricsDemog_df)
}


# Function for making connMetrics_demographics dataframe (using covbat harmonized output)
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_connMetricsDemog_covbat <- function(atlas_network, metric, demog_csv, subj_id, dataset, conn_type){
  subxparcel.matrix <- read.csv(sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/%2$s_subxparcel_matrix_%3$s_%4$s_covbat.csv", dataset, metric, atlas_network, conn_type))

  print(paste(atlas_network, metric, "file loaded"))
  subxparcel.matrix$subject <- demographics$subject
  subxparcel.matrix <- subxparcel.matrix %>% dplyr::relocate(subject)
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  connMetricsDemog_df <- connMetricsDemog_df[,-2]
  write.csv(connMetricsDemog_df, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample_%4$s_covbat.csv", dataset, metric, atlas_network, conn_type), row.names = F)
  return(connMetricsDemog_df)
}

   

# Function for **Fitting GAMs** for GBC, BNC, and WNC
# Fit GAM (func_conn_metric ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name, conn_type) {
  parcel.labels <- parcel.labels$label
  gam.age <- matrix(data=NA, nrow=length(parcel.labels), ncol=9) #empty matrix to save gam.fit output to
  
  for(row in c(1:length(parcel.labels))){ #for each region
    region <- parcel.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name, network_parcellation), conn_type = conn_type, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  write.csv(gam.age, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s_%5$s.csv", dataset_name, metric, atlas_name, network_parcellation, conn_type), row.names = F, quote = F)
  
  return(gam.age)
}




# Function for **Fitting GAMs** for GBC, BNC, and WNC (using covbat harmonized output)
# Fit GAM (func_conn_metric ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs_covbat <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name, conn_type) {
  parcel.labels <- parcel.labels$label
  gam.age <- matrix(data=NA, nrow=length(parcel.labels), ncol=9) #empty matrix to save gam.fit output to
  
  for(row in c(1:length(parcel.labels))){ #for each region
    region <- parcel.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name, network_parcellation), conn_type = conn_type, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  write.csv(gam.age, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s_%5$s_covbat.csv", dataset_name, metric, atlas_name, network_parcellation, conn_type), row.names = F, quote = F)
  return(gam.age)
}

  

 
# Function for estimating GAM smooths based on model-predicted data and save out predicted y data
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
estimate_GAMsmooths <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name, conn_type) {
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
  
  
  write.csv(gam.smooths, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s%4$s_%5$s.csv", dataset_name, metric, atlas_name, network_parcellation, conn_type), row.names = F, quote = F)

  return(gam.smooths)
  
}
 


# Function for estimating GAM smooths based on model-predicted data and save out predicted y data (using covbat harmonized output)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
estimate_GAMsmooths_covbat <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name, conn_type) {
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
  
  
  write.csv(gam.smooths, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s%4$s_%5$s_covbat.csv", dataset_name, metric, atlas_name, network_parcellation, conn_type), row.names = F, quote = F)
  return(gam.smooths)
}


## Function to calculate Region-wise Fitted GBC values  
# derived from https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd#L699-L709)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")

dev_fitted <- function(parcel.labels, metric, atlas_name, dataset_name, conn_type){
  np <- 200 #number of ages to get the fitted at
  npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)
  
  gam.fitted.atlas <- matrix(data=NA, ncol=9)
  colnames(gam.fitted.atlas) <- c("age","fitted","se","lower","upper","significant","significant.fit","index","label")
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  
  if (atlas_name=="glasser") {
    to_replace_in_glasser <- setdiff(SAaxis$label, unique(parcel.labels))
    glasser_labels_to_change <- gsub("-", '.', to_replace_in_glasser)
    index <- which(SAaxis$label %in% to_replace_in_glasser)
    SAaxis$label[index] <- glasser_labels_to_change
  } 
  
  if(metric == "GBC" & str_detect(atlas_name, "x")){
    atlas_name <- gsub("x(.*)", "", atlas_name)
  } else {
    atlas_name <- atlas_name
  }
  for(row in c(1:nrow(SAaxis))){ #for each region
    region <- SAaxis$label[row] 
    GAM.fitted.output <- gam.smooth.predict_posterior(measure = metric, atlas = atlas_name, dataset = dataset_name, 
                                              region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes", 
                                              knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_fits = FALSE) #run the gam.fitted.atlas function to get true model fitted
    GAM.fitted.output$index <- as.numeric(rep(x=row, np)) #region index
    GAM.fitted.output$label <- as.character(rep(x=region, np)) #region label
    gam.fitted.atlas <- rbind(gam.fitted.atlas, GAM.fitted.output)
  }
  
  gam.fitted.atlas <- gam.fitted.atlas[-1,] #remove empty initialization row
  gam.fitted.atlas <- left_join(gam.fitted.atlas, SAaxis, by="label", sort = F)
  filename <- paste0("gam.", metric, ".fitted.", atlas_name)
  write.csv(gam.fitted.atlas, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/GAM/%3$s_%4$s.csv", dataset_name, metric, filename, conn_type), row.names = F, quote = F)
  return(gam.fitted.atlas)
}
 

## Function to calculate Region-wise Fitted GBC values (using covbat harmonized output)
# derived from https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd#L699-L709)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")
dev_fitted_covbat <- function(parcel.labels, metric, atlas_name, dataset_name, conn_type){
  np <- 200 #number of ages to get the fitted at
  npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)
  
  gam.fitted.atlas <- matrix(data=NA, ncol=9)
  colnames(gam.fitted.atlas) <- c("age","fitted","se","lower","upper","significant","significant.fit","index","label")
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  
  if (atlas_name=="glasser") {
    to_replace_in_glasser <- setdiff(SAaxis$label, unique(parcel.labels))
    glasser_labels_to_change <- gsub("-", '.', to_replace_in_glasser)
    index <- which(SAaxis$label %in% to_replace_in_glasser)
    SAaxis$label[index] <- glasser_labels_to_change
  } 
  
  if(metric == "GBC" & str_detect(atlas_name, "x")){
    atlas_name <- gsub("x(.*)", "", atlas_name)
  } else {
    atlas_name <- atlas_name
  }
  for(row in c(1:nrow(SAaxis))){ #for each region
    region <- SAaxis$label[row] 
    GAM.fitted.output <- gam.smooth.predict_posterior(measure = metric, atlas = atlas_name, dataset = dataset_name, 
                                                      region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes", 
                                                      knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_fits = FALSE) #run the gam.fitted.atlas function to get true model fitted
    GAM.fitted.output$index <- as.numeric(rep(x=row, np)) #region index
    GAM.fitted.output$label <- as.character(rep(x=region, np)) #region label
    gam.fitted.atlas <- rbind(gam.fitted.atlas, GAM.fitted.output)
  }
  
  gam.fitted.atlas <- gam.fitted.atlas[-1,] #remove empty initialization row
  gam.fitted.atlas <- left_join(gam.fitted.atlas, SAaxis, by="label", sort = F)
  filename <- paste0("gam.", metric, ".fitted.", atlas_name)
  write.csv(gam.fitted.atlas, sprintf("/cbica/projects/network_replication/revisions/output/%1$s/%2$s/GAM/%3$s_%4$s_covbat.csv", dataset_name, metric, filename, conn_type), row.names = F, quote = F)
  return(gam.fitted.atlas)
}





