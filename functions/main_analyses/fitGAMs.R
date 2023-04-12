# functions for fitting GAMs 
library(dplyr)
library(janitor)
library(purrr)
source("/cbica/projects/network_replication/adapted_Rscripts/GAM_functions.R")

# load parcellated S-A axis
glasser_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/glasser_SAaxis.csv")
gordon_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/gordon_SAaxis.csv")

schaefer200x7_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x7_SAaxis.csv")
schaefer200x7_SAaxis$label <- gsub("7Network", "Network", schaefer200x7_SAaxis$label) 

schaefer200x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x17_SAaxis.csv")
schaefer200x17_SAaxis$label <- gsub("17Network", "Network", schaefer200x17_SAaxis$label) 
schaefer200_SAaxis <- schaefer200x17_SAaxis 

schaefer400x7_SAaxis  <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer400x7_SAaxis.csv")
schaefer400x7_SAaxis$label <-  gsub("7Networks", "Networks", schaefer400x7_SAaxis$label)

schaefer400x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer400x17_SAaxis.csv")
schaefer400x17_SAaxis$label <- gsub("17Network", "Network", schaefer400x17_SAaxis$label)
schaefer400_SAaxis <- schaefer400x17_SAaxis


# load parcel labels  
glasser.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/glasser360_regionlist_final.csv")

gordon.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/gordon_regionlist_final.csv")

schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")

schaefer200x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
schaefer200.parcel.labels <- schaefer200x17.parcel.labels #gbc

schaefer400x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x7_regionlist_final.csv")

schaefer400x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x17_regionlist_final.csv")
schaefer400.parcel.labels <- schaefer400x17.parcel.labels #gbc




# Function for making connMetrics_demographics.csv
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_connMetricsDemog <- function(atlas_network, metric, demog_csv, subj_id, dataset){
  subxparcel.matrix <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s_subxparcel_matrix_%3$s.csv", dataset, metric, atlas_network))
  print(paste(atlas_network, metric, "file loaded"))
  #subxparcel.matrix <- remove_empty(subxparcel.matrix, which = "cols")
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  #write.csv(connMetricsDemog_df, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample.csv", dataset, metric, atlas_network), row.names = F)
  return(connMetricsDemog_df)
}


# Function for making connMetrics_demographics.csv with covbat
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_connMetricsDemog_covbat <- function(atlas_network, metric, demog_csv, subj_id, dataset){
  subxparcel.matrix <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s_subxparcel_matrix_%3$s_covbat.csv", dataset, metric, atlas_network))

  print(paste(atlas_network, metric, "file loaded"))
  subxparcel.matrix$subject <- demographics$subject
  subxparcel.matrix <- subxparcel.matrix %>% dplyr::relocate(subject)
  #subxparcel.matrix <- remove_empty(subxparcel.matrix, which = "cols")
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  
  #write.csv(connMetricsDemog_df, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample_covbat.csv", dataset, metric, atlas_network), row.names = F)
  return(connMetricsDemog_df)
}


# Function for making edges_demographics.csv
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_EdgeDemog <- function(atlas_network, demog_csv, subj_id, dataset){
  subxedge.matrix <- readRDS(sprintf("/cbica/projects/network_replication/output/%1$s/edge/subxedge_%2$s.RData", dataset, atlas_network))
  print("file loaded")
  subxedge.matrix[subxedge.matrix == 1] <- NA
  subxedge.matrix <- remove_empty(subxedge.matrix, which = "cols") # remove columns that are all ==1 (connection to self)
  print("rm columns of self-connections")
  edgeDemog_df <- merge(subxedge.matrix, demog_csv, by=subj_id)
  saveRDS(edgeDemog_df, sprintf("/cbica/projects/network_replication/output/%1$s/edge/%2$s_demographics_finalsample.RData", dataset, atlas_network))
  print(paste(atlas_network, "done"))
  return(edgeDemog_df)
}


# Function for making edges_demographics.csv
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_EdgeDemog_covbat <- function(atlas_network, demog_csv, subj_id, dataset){
  subxedge.matrix <- readRDS(sprintf("/cbica/projects/network_replication/output/%1$s/edge/subxedge_%2$s_covbat.RData", dataset, atlas_network))
  print("file loaded")
  subxedge.matrix[subxedge.matrix == 1] <- NA
  subxedge.matrix <- remove_empty(subxedge.matrix, which = "cols") # remove columns that are all ==1 (connection to self)
  print("rm columns of self-connections")
  subxedge.matrix$subject <- rownames(subxedge.matrix) 
  subxedge.matrix <- subxedge.matrix %>% relocate(subject)
  edgeDemog_df <- merge(subxedge.matrix, demog_csv, by=subj_id)
  saveRDS(edgeDemog_df, sprintf("/cbica/projects/network_replication/output/%1$s/edge/%2$s_demographics_finalsample_covbat.RData", dataset, atlas_network))
  print(paste(atlas_network, "done"))
  return(edgeDemog_df)
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
  ##write.csv(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  
  return(gam.age)
}


# Function for **Fitting GAMs** for GBC, BNC, and WNC
# Fit GAM (func_conn_metric ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs_intercept <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name) {
  parcel.labels <- parcel.labels$label
  gam.age <- matrix(data=NA, nrow=length(parcel.labels), ncol=10) #empty matrix to save gam.fit output to
  
  for(row in c(1:length(parcel.labels))){ #for each region
    region <- parcel.labels[row] 
    GAM.RESULTS <- gam.fit_intercept(measure = metric, atlas = paste0(atlas_name, network_parcellation), dataset = dataset_name, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","Intercept", "GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:10)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  ##write.csv(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s_intercept.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  
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
  #write.csv(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s_covbat.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  return(gam.age)
}


# Function for **Fitting GAMs** for parcel-by-parcel edges
# Fit GAM (edge ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param edge.labels A vector of edge.labels
# @param metric A character string of connectivity metric (e.g. "edge")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs_edge <- function(edge.labels, metric, atlas_name, dataset_name) {
  gam.age <- matrix(data=NA, nrow=length(edge.labels), ncol=9) #empty matrix to save gam.fit output to
  for(row in c(1:length(edge.labels))){ #for each region
    edge_name <- edge.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name), dataset = dataset_name, region = edge_name, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS
    print(paste(row, "/", length(edge.labels)))}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  saveRDS(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/edge/GAM/GAMresults.%2$s.age.%3$s.RData", dataset_name, metric, atlas_name))
  #saveRDS(gam.age, sprintf("/Users/audluo/Desktop/%1$s_GAMresults.%2$s.age.%3$s.RData", dataset_name, metric, atlas_name))
  return(gam.age)
}


# Function for **Fitting GAMs** for parcel-by-parcel edges - with covbat
# Fit GAM (edge ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param edge.labels A vector of edge.labels
# @param metric A character string of connectivity metric (e.g. "edge")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs_edge_covbat <- function(edge.labels, metric, atlas_name, dataset_name) {
  gam.age <- matrix(data=NA, nrow=length(edge.labels), ncol=9) #empty matrix to save gam.fit output to
  for(row in c(1:length(edge.labels))){ #for each region
    edge_name <- edge.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name), dataset = dataset_name, region = edge_name, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS
    print(paste(row, "/", length(edge.labels)))}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  saveRDS(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/edge/GAM/GAMresults.%2$s.age.%3$s_covbat.RData", dataset_name, metric, atlas_name))
  return(gam.age)
}


# Function for **Fitting GAMs** for parcel-by-parcel edges -- using bam()
# Fit GAM (edge ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param edge.labels A vector of edge.labels
# @param metric A character string of connectivity metric (e.g. "edge")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs_edge_bam <- function(edge.labels, metric, atlas_name, dataset_name) {
  gam.age <- matrix(data=NA, nrow=length(edge.labels), ncol=9) #empty matrix to save gam.fit output to
  for(row in c(1:length(edge.labels))){ #for each region
    edge_name <- edge.labels[row] 
    GAM.RESULTS <- bam.fit(measure = metric, atlas = paste0(atlas_name), dataset = dataset_name, region = edge_name, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS
    print(paste(row, "/", length(edge.labels)))}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  #saveRDS(gam.age, sprintf("/cbica/projects/network_replication/output/%1$s/edge/GAM/GAMresults.%2$s.age.%3$s.RData", dataset_name, metric, atlas_name))
   saveRDS(gam.age, sprintf("/Users/audluo/Desktop/%1$s_GAMresults.%2$s.age.%3$s.RData", dataset_name, metric, atlas_name))
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
  
  
  #write.csv(gam.smooths, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s%4$s.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)

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
  
  
  #write.csv(gam.smooths, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s%4$s_covbat.csv", dataset_name, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  return(gam.smooths)
}



## Function to calculate Region-wise Developmental Derivatives 
# derived from https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd#L699-L709)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")


dev_fitted <- function(parcel.labels, metric, atlas_name, dataset_name){
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
  #write.csv(gam.fitted.atlas, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/%3$s.csv", dataset_name, metric, filename), row.names = F, quote = F)
  return(gam.fitted.atlas)
}
 

## Function to calculate Region-wise Developmental Derivatives 
# derived from https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd#L699-L709)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")


dev_fitted_covbat <- function(parcel.labels, metric, atlas_name, dataset_name){
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
  #write.csv(gam.fitted.atlas, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/%3$s_covbat.csv", dataset_name, metric, filename), row.names = F, quote = F)
  return(gam.fitted.atlas)
}



## Function to calculate Region-wise Developmental Derivatives 
# derived from https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd#L699-L709)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")


dev_derivatives <- function(parcel.labels, metric, atlas_name, dataset_name){
  np <- 200 #number of ages to get the derivative at
  npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)
  
  gam.derivatives.atlas <- matrix(data=NA, ncol=9)
  colnames(gam.derivatives.atlas) <- c("age","derivative","se","lower","upper","significant","significant.derivative","index","label")
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
    GAM.DERIVATIVES.output <- gam.derivatives(measure = metric, atlas = atlas_name, dataset = dataset_name, 
                                              region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes", 
                                              knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_derivatives = FALSE) #run the gam.derivatives.atlas function to get true model derivatives
    GAM.DERIVATIVES.output$index <- as.numeric(rep(x=row, np)) #region index
    GAM.DERIVATIVES.output$label <- as.character(rep(x=region, np)) #region label
    gam.derivatives.atlas <- rbind(gam.derivatives.atlas, GAM.DERIVATIVES.output)
  }
  
  gam.derivatives.atlas <- gam.derivatives.atlas[-1,] #remove empty initialization row
  gam.derivatives.atlas <- left_join(gam.derivatives.atlas, SAaxis, by="label", sort = F)
  filename <- paste0("gam.", metric, ".derivatives.", atlas_name)
  #write.csv(gam.derivatives.atlas, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/%3$s.csv", dataset_name, metric, filename), row.names = F, quote = F)
  return(gam.derivatives.atlas)
}




## Function to calculate Region-wise Developmental Derivatives- covbat
# derived from https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd#L699-L709)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")


dev_derivatives_covbat <- function(parcel.labels, metric, atlas_name, dataset_name){
  np <- 200 #number of ages to get the derivative at
  npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)
  
  gam.derivatives.atlas <- matrix(data=NA, ncol=9)
  colnames(gam.derivatives.atlas) <- c("age","derivative","se","lower","upper","significant","significant.derivative","index","label")
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
    GAM.DERIVATIVES.output <- gam.derivatives(measure = metric, atlas = atlas_name, dataset = dataset_name, 
                                              region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes", 
                                              knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_derivatives = FALSE) #run the gam.derivatives.atlas function to get true model derivatives
    GAM.DERIVATIVES.output$index <- as.numeric(rep(x=row, np)) #region index
    GAM.DERIVATIVES.output$label <- as.character(rep(x=region, np)) #region label
    gam.derivatives.atlas <- rbind(gam.derivatives.atlas, GAM.DERIVATIVES.output)
  }
  
  gam.derivatives.atlas <- gam.derivatives.atlas[-1,] #remove empty initialization row
  gam.derivatives.atlas <- left_join(gam.derivatives.atlas, SAaxis, by="label", sort = F)
  filename <- paste0("gam.", metric, ".derivatives.", atlas_name)
  #write.csv(gam.derivatives.atlas, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/%3$s_covbat.csv", dataset_name, metric, filename), row.names = F, quote = F)
  return(gam.derivatives.atlas)
}


## Function to calculate Region-wise Posterior Smooth Derivatives Correlation with SA Axis 
#function to estimate npd posterior draw derivatives for each region 
#and compute the correlation between regional derivative and regional S-A axis rank at each age, for each draw
#executed via rerun below

np <- 200 #number of ages to get the derivative at
npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)

 
compute_axis_correlation <- function(metric, atlas_name, dataset_name){ 
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  parcel.labels <- get(paste0(atlas_name, ".parcel.labels"))
  print("Get posterior derivatives")
  gam.derivatives.atlas <- map_dfr(parcel.labels, function(x){
    gam.derivatives(measure = metric, atlas = atlas_name, dataset = dataset_name, region = as.character(x), smooth_var = "age", covariates = "sex + meanFD_avgSes", knots = 3, set_fx = TRUE, 
                    draws = npd, increments = np, return_posterior_derivatives = TRUE)}) #run gam.derivatives to get simulated derivatives
  # Compute and save correlations with S-A axis
  gam.derivatives.atlas <- left_join(SAaxis, gam.derivatives.atlas, by="label", sort=F) #assign axis rank to each label
  corr_values <- gam.derivatives.atlas %>%
    group_by(draw,age) %>%
    do(SAcorrelation = cor(as.numeric(.$SA.axis_rank), as.numeric(.$posterior.derivative), method=c("spearman"))) %>% #correlation between parcel rank and derivative at each age for each draw
    unnest(cols = c(SAcorrelation))
  corr_values.wide <- corr_values %>% pivot_wider(names_from = "draw", values_from = "SAcorrelation", names_sort = FALSE)
  corr_values.wide <- corr_values.wide %>% select(contains("draw"))
  print(paste(metric, atlas_name, "finished"))
  return(corr_values.wide)
}

make_deriv.SAaxis.posteriorcor <- function(metric, atlas_name, dataset_name){
  deriv.SAaxis.posteriorcor <- rerun(10, compute_axis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; each column is a draw and has the correlation between parcel derivative and parcel axis ranking at each age (row)
  #write.table(deriv.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorderivative_correlation_byage_%3$s.csv", dataset_name, metric, atlas_name), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
  rm(deriv.SAaxis.posteriorcor)
  gc()
}





## Function to calculate Region-wise Posterior Smooth Derivatives Correlation with SA Axis 
#function to estimate npd posterior draw derivatives for each region 
#and compute the correlation between regional derivative and regional S-A axis rank at each age, for each draw
#executed via rerun below

np <- 200 #number of ages to get the derivative at
npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)


compute_axis_correlation <- function(metric, atlas_name, dataset_name){ 
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  parcel.labels <- get(paste0(atlas_name, ".parcel.labels"))
  print("Get posterior derivatives")
  gam.derivatives.atlas <- map_dfr(parcel.labels, function(x){
    gam.derivatives(measure = metric, atlas = atlas_name, dataset = dataset_name, region = as.character(x), smooth_var = "age", covariates = "sex + meanFD_avgSes", knots = 3, set_fx = TRUE, 
                    draws = npd, increments = np, return_posterior_derivatives = TRUE)}) #run gam.derivatives to get simulated derivatives
  # Compute and save correlations with S-A axis
  gam.derivatives.atlas <- left_join(SAaxis, gam.derivatives.atlas, by="label", sort=F) #assign axis rank to each label
  corr_values <- gam.derivatives.atlas %>%
    group_by(draw,age) %>%
    do(SAcorrelation = cor(as.numeric(.$SA.axis_rank), as.numeric(.$posterior.derivative), method=c("spearman"))) %>% #correlation between parcel rank and derivative at each age for each draw
    unnest(cols = c(SAcorrelation))
  corr_values.wide <- corr_values %>% pivot_wider(names_from = "draw", values_from = "SAcorrelation", names_sort = FALSE)
  corr_values.wide <- corr_values.wide %>% select(contains("draw"))
  print(paste(metric, atlas_name, "finished"))
  return(corr_values.wide)
}

make_deriv.SAaxis.posteriorcor<- function(metric, atlas_name, dataset_name){
  deriv.SAaxis.posteriorcor <- rerun(10, compute_axis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; each column is a draw and has the correlation between parcel derivative and parcel axis ranking at each age (row)
  #write.table(deriv.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorderivative_correlation_byage_%3$s.csv", dataset_name, metric, atlas_name), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
  rm(deriv.SAaxis.posteriorcor)
  gc()
}
make_deriv.SAaxis.posteriorcor_covbat <- function(metric, atlas_name, dataset_name){
  deriv.SAaxis.posteriorcor <- rerun(10, compute_axis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; each column is a draw and has the correlation between parcel derivative and parcel axis ranking at each age (row)
  #write.table(deriv.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorderivative_correlation_byage_%3$s_covbat.csv", dataset_name, metric, atlas_name), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
  rm(deriv.SAaxis.posteriorcor)
  gc()
}




## Function to calculate GAM Fitted Value Predictions alignment with S-A Axis
#function to estimate npd posterior draw derivs for each region 
#and compute the correlation between regional GBC fitted value and regional S-A axis rank at each age, for each draw
#executed via rerun below

np <- 200 #number of ages to get the deriv at
npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)


compute_GBCaxis_correlation <- function(metric, atlas_name, dataset_name){ 
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  parcel.labels <- get(paste0(atlas_name, ".parcel.labels"))
  parcel.labels <- parcel.labels$label
  print("Get posterior predicted fitted values")
  gam.fits.atlas <- map_dfr(parcel.labels, function(x){
    gam.smooth.predict_posterior(measure = metric, atlas = atlas_name, dataset = dataset_name, region = as.character(x), smooth_var = "age", covariates = "sex + meanFD_avgSes", knots = 3, set_fx = TRUE, 
                                 draws = npd, increments = np, return_posterior_fits = TRUE)}) #run gam.fits to get simulated fits
  # Compute and save correlations with S-A axis
  gam.fits.atlas <- left_join(SAaxis, gam.fits.atlas, by="label", sort=F) #assign axis rank to each label
  corr_values <- gam.fits.atlas %>%
    group_by(draw,age) %>%  
    do(SAcorrelation = cor(as.numeric(.$SA.axis_rank), as.numeric(.$posterior.fits), method=c("spearman"))) %>% #correlation between parcel rank and fit at each age for each draw 
    # SAcorrelation is computed  by taking the correlation between parcel rank and parcel posterior fit at a given draw and age. 
    # So for draw 1 of age 8, take the SAcorrelation for all the parcels (a correlation of length(atlas) parcels). 
    # Then for draw 1 of age 8.067, again take the correlation of all the parcels. Etc.
    unnest(cols = c(SAcorrelation))
  corr_values.wide <- corr_values %>% pivot_wider(names_from = "draw", values_from = "SAcorrelation", names_sort = FALSE)
  corr_values.wide <- corr_values.wide %>% select(contains("draw"))
  print(paste(metric, atlas_name, "finished"))
  rm(gam.fits.atlas)
  gc()
  return(corr_values.wide)
}


make_fits.SAaxis.posteriorcor <- function(metric, atlas_name, dataset_name){
  fits.SAaxis.posteriorcor <- rerun(10, compute_GBCaxis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; 
  # each column is a draw and has the correlation between parcel fits and parcel axis ranking at each age (row)
  # basically is a df of SA correlations across all parcels, for a given draw and age
  saveRDS(fits.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorfits_correlation_byage_%3$s.RData", dataset_name, metric, atlas_name))
  rm(fits.SAaxis.posteriorcor)
  gc()
}

make_fits.SAaxis.posteriorcor_covbat <- function(metric, atlas_name, dataset_name){
  fits.SAaxis.posteriorcor <- rerun(10, compute_GBCaxis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; each column is a draw and has the correlation between parcel fits and parcel axis ranking at each age (row)
  saveRDS(fits.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorfits_correlation_byage_%3$s_covbat.RData", dataset_name, metric, atlas_name))
  rm(fits.SAaxis.posteriorcor)
  gc()
}
 

