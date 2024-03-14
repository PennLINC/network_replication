# note that lib.loc has to be specified since this file is linked to in Rscripts submitted as jobs onto cubic (i.e. /Rscripts/functions/main_analyses/compute_GBC.R)
library(cifti, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")

args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
conn_type = args[2] # thresholded or absvalue

indir <- sprintf("/cbica/projects/network_replication/revisions/input/%1$s/connMatricesData/connectivity_matrices/", dataset)
outdir <- sprintf("/cbica/projects/network_replication/revisions/output/%1$s/GBC/", dataset)

schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")

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
  connect.matrix <- readRDS(paste0(indir, sprintf("%1$s_ConnMatrices_%2$s.RData", subject, conn_type)))
  
  #compute average connectivity 
  GBC <- as.array(rowMeans(connect.matrix, na.rm=TRUE))
  return(GBC)
}




if (dataset == "NKI"){
  participants <- read.csv("/cbica/projects/network_replication/input/NKI/sample_selection/NKI_demographics_finalsample_20230629.csv")
  participants <- gsub("sub-", "", participants$sub) 
} else if (dataset=="HCPD") {
  participants <- read.csv("/cbica/projects/network_replication/input/HCPD/sample_selection/HCPD_demographics_finalsample_20221226.csv")
  participants <- gsub("HCD", "sub-", participants$src_subject_id)
  
} else if (dataset=="PNC") {
  participants <- read.csv("/cbica/projects/network_replication/input/PNC/sample_selection/PNC_demographics_finalsample_20230629.csv")
  participants <- participants$sub
} else if (dataset == "HBN") {
  participants <- read.csv("/cbica/projects/network_replication/input/HBN/sample_selection/HBN_demographics_finalsample_20230629.csv")
  participants <- participants$sub
} else {
  participants <- NA
  print("Provide valid dataset")
}



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

write.csv(GBC.subxparcel.matrix.schaefer200x7, sprintf("%1$sGBC_subxparcel_matrix_schaefer200_%2$s.csv", outdir, conn_type), row.names=F, quote=F)




