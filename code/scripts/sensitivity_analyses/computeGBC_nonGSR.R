library(dplyr)
library(magrittr)
library(rjson)
library(stringr)

atlas = "schaefer200x7"

################## 
# Set Directories 
################## 
config_data <- fromJSON(file="/cbica/projects/network_replication/manuscript/code/config_PNC.json")
data_root <- config_data$data_root
conn_matrices_dir <- paste0(data_root, "CPAC_connMatricesData/connectivity_matrices_nonGSR/")
outputs_root <- config_data$outputs_root
sample_selection_dir <- config_data$sample_selection_data_dir
metric_output_dir <- paste0(outputs_root, "sensitivity_analyses/GBC/")

if (!dir.exists(outputs_root)) {
  # If directory doesn't exist, create it
  dir.create(outputs_root, recursive = TRUE)
  print(paste("Directory", outputs_root, "created."))
} else {
  print(paste("Directory",outputs_root, "already exists."))
}

if (!dir.exists(metric_output_dir)) {
  # If directory doesn't exist, create it
  dir.create(metric_output_dir, recursive = TRUE)
  print(paste("Directory", metric_output_dir, "created."))
} else {
  print(paste("Directory", metric_output_dir, "already exists."))
}

 
################## 
# Define Functions
################## 
# Function for making empty dataframes for connectivity metric outputs (schaefer 200 only)
# @param subject_list list of participant ID's
make_output_dfs <- function(subject_list){
  n <- 201
  regionheaders <- as.character(schaefer200x7.parcel.labels$label)
  subxparcel.matrix  <- matrix(data = NA, nrow = length(subject_list), ncol = n)
  demoheaders <- c("subject")
  colheaders <- as.matrix(c(demoheaders,regionheaders))
  colnames(subxparcel.matrix) <- colheaders
  return(subxparcel.matrix)
  
}


# Function for Computing Global Brain Connectivity  
# @param subject id of subject of interest
computeGBC <- function(subject){
  
  #read in connectivity matrix 
  connect.matrix <- readRDS(paste0(conn_matrices_dir, sprintf("%1$s_ConnMatrices.RData", subject)))
  
  #compute average connectivity 
  GBC <- as.array(rowMeans(connect.matrix, na.rm=TRUE))
  return(GBC)
}


################## 
# Read files 
################## 
schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
# nonGSR participants 
participants <- read.csv("/cbica/projects/network_replication/revisions/input/PNC/CPAC_sample_selection/CPAC_nonGSR_PNC_demographics_finalsample_20231130.csv")
participants <- participants$sub
 

################## 
# Compute Metric
################## 
# make output df: 
##GBC.subxparcel.matrix.schaefer200x7 
GBC.subxparcel.matrix.schaefer200x7 <- make_output_dfs(participants) 

# compute GBC
for(sub in c(1:length(participants))){
  subjectID=as.character(participants[sub])
  id.data <- computeGBC(subjectID)
  GBC.subxparcel.matrix.schaefer200x7[sub,] <- cbind(subjectID, t(id.data)) #update the name of the df every iteration to GBC.subxparcel.matrix.[atlas]
  print(paste(sub, "/", length(participants), "-", subjectID))
  
}

write.csv(GBC.subxparcel.matrix.schaefer200x7, sprintf("%1$sGBC_subxparcel_matrix_nonGSR.csv", metric_output_dir), row.names=F, quote=F)
