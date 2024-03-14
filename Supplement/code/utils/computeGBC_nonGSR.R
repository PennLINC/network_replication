library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")

atlas = "schaefer200x7"
dataset="PNC"


indir <- sprintf("/cbica/projects/network_replication/revisions/input/PNC/CPAC_connMatricesData/connectivity_matrices/", dataset)
outdir <- sprintf("/cbica/projects/network_replication/revisions/output/PNC/GBC/", dataset)

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


# nonGSR participants 
participants <- read.csv("/cbica/projects/network_replication/revisions/input/PNC/CPAC_sample_selection/CPAC_nonGSR_PNC_demographics_finalsample_20231130.csv")
participants <- participants$sub
 

# Function for Computing Global Brain Connectivity  
# @param subject id of subject of interest
computeGBC <- function(subject){
  
  #read in connectivity matrix 
  connect.matrix <- readRDS(paste0(indir, sprintf("%1$s_ConnMatrices.RData", subject)))
  
  #compute average connectivity 
  GBC <- as.array(rowMeans(connect.matrix, na.rm=TRUE))
  return(GBC)
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

write.csv(GBC.subxparcel.matrix.schaefer200x7, sprintf("%1$sGBC_subxparcel_matrix_nonGSR.csv", outdir), row.names=F, quote=F)
