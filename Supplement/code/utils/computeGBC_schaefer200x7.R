
# note that lib.loc has to be specified since this file is linked to in Rscripts submitted as jobs onto cubic (i.e. /Rscripts/functions/main_analyses/compute_GBC.R)
library(cifti, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")

args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]

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
  atlas <- "schaefer200"
  #read in connectivity matrix
  if(dataset == "PNC" | dataset == "HCPD" | dataset == "HBN") {
    connect.matrix <- readRDS(sprintf("/cbica/projects/network_replication/input/%1$s/connMatricesData/connectivity_matrices/%2$s_ConnMatrices.RData", dataset, subject))
    
    if(atlas=="schaefer200" | atlas=="schaefer400") {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
  } else if(dataset == "NKI") {  
    connect.matrix <- readRDS(sprintf("/cbica/projects/network_replication/input/NKI/connMatricesData/connectivity_matrices/%1$s_ConnMatrices.RData", subject))
    ses_name <- str_extract(names(connect.matrix), "[A-Z]{3}1")[1]
    if(atlas=="schaefer200" | atlas=="schaefer400") {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(ses_name, "_", atlas_name, "_conn")]]
  } else {
    print("Provide valid dataset")
    
  }
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

write.csv(GBC.subxparcel.matrix.schaefer200x7, sprintf("%1$sGBC_subxparcel_matrix_schaefer200x7_orig.csv", outdir), row.names=F, quote=F)




