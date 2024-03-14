# functions for computing global brain connectivity (GBC) for rest-only sensitivity analyses

library(cifti)
library(dplyr)
library(magrittr)
library(rjson)
library(stringr)

# Function for Computing Global Brain Connectivity  
# @param subject rbcid of subject of interest
# @param atlas A character string, name of atlas of interest (e.g. schaefer200)
# @param dataset A character string, name of dataset
computeGBC_restOnly <- function(subject, atlas, dataset){
  
  #read in connectivity matrix and mask
  connect.matrix <- readRDS(sprintf("/cbica/projects/network_replication/manuscript/input/%1$s/connMatricesData/connectivity_matrices_restOnly/%2$s_ConnMatrices.RData",dataset, subject))
  if(atlas=="schaefer200" | atlas=="schaefer400") {
    atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
  } else {
    atlas_name <- atlas
  } 
  connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
 
   
  #compute average connectivity 
  GBC <- as.array(rowMeans(connect.matrix))
  return(GBC)
}
 
  
 
# needed for following function 
glasser.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/glasser360_regionlist_final.csv")

gordon.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/gordon_regionlist_final.csv")

schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")

schaefer200x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
 
schaefer400x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x7_regionlist_final.csv")
 
schaefer400x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x17_regionlist_final.csv")
 
# Function for making empty dataframes for connectivity metric outputs
# @param subject_list list of participant ID's
# @param atlas A character string, name of atlas of interest (gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", or "schaefer400x17")
make_output_dfs <- function(subject_list, atlas){
  if(atlas == "glasser"){
    n <- 361
    regionheaders <- as.character(glasser.parcel.labels$label)
  } else if (atlas =="gordon"){
    n <- 334
    regionheaders <- as.character(gordon.parcel.labels$label)
  } else if (atlas =="schaefer200x7"| atlas=="schaefer200"){
    n <- 201
    regionheaders <- as.character(schaefer200x7.parcel.labels$label)
  } else if (atlas =="schaefer200x17" ) {
    n <- 201
    regionheaders <- as.character(schaefer200x17.parcel.labels$label)
  } else if (atlas =="schaefer400x7" | atlas=="schaefer400") {
    n <- 401
    regionheaders <- as.character(schaefer400x7.parcel.labels$label)
  } else if (atlas =="schaefer400x17") {
    n <- 401
    regionheaders <- as.character(schaefer400x17.parcel.labels$label)
  }
  subxparcel.matrix  <- matrix(data = NA, nrow = length(subject_list), ncol = n)
  demoheaders <- c("subject")
  colheaders <- as.matrix(c(demoheaders,regionheaders))
  colnames(subxparcel.matrix) <- colheaders
  return(subxparcel.matrix)
  
}
 



 
 