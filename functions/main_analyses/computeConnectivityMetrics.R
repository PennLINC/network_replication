# functions for computing global brain connectivity (GBC), between-network coupling (BNC), and within-network coupling (WNC)
# use computeGBC.R, computeBNC.R, computeWNC.R, and computeEdge.R to run jobs.
# This file just consolidates all the code

#library(ciftiTools)
#ciftiTools.setOption('wb_path', '/Applications/workbench/')
#library(gifti)
library(cifti, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")


# Function for Computing Global Brain Connectivity - updated
# @param subject rbcid of subject of interest
# @param atlas A character string, name of atlas of interest (e.g. glasser, gordon, schaefer200, or schaefer400)
# @param dataset A character string, name of dataset
computeGBC <- function(subject, atlas, dataset){
  
  #read in connectivity matrix and mask   
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
  GBC <- as.array(rowMeans(connect.matrix, na.rm = TRUE))
  return(GBC)
}
 

# Function for Computing Between-Network or Within-Network Connectivity
# @param subject rbcid of subject of interest
# @param atlas A character string, name of atlas of interest (gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", or "schaefer400x17")
# @param metric A character string, name of connectivity metric (either "BNC" or "WNC")
# @param dataset A character string, name of dataset
computeBNC_WNC <- function(subject, atlas, metric, dataset){
  print(paste("loading", atlas, "connectivity matrices"))
  #read in connectivity matrix  
  if (dataset == "PNC" | dataset == "HCPD" | dataset == "HBN") {
    connect.matrix <- readRDS(sprintf("/cbica/projects/network_replication/input/%1$s/connMatricesData/connectivity_matrices/%2$s_ConnMatrices.RData",dataset, subject))
    if (str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
  } else if(dataset == "NKI") {  
    connect.matrix <- readRDS(sprintf("/cbica/projects/network_replication/input/NKI/connMatricesData/connectivity_matrices/%1$s_ConnMatrices.RData", subject))
    ses_name <- str_extract(names(connect.matrix), "[A-Z]{3}1")[1]
    if(str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(ses_name, "_", atlas_name, "_conn")]]
  } else {
    print("Provide valid dataset")
  }
  if (atlas == "gordon"){
    parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/gordon_regionlist_final.csv")
    communityAffil <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Gordon333_Network/gordon333CommunityAffiliation_final.csv")
    parcel.labels <- as.data.frame(cbind(parcel.labels$label, communityAffil$CommAffil))
    names(parcel.labels) <- c("label", "network")
    ## make dataframe where rownames are the commAffil number, and then need to transpose V1 to be rows 
    rows <- as.data.frame(t(parcel.labels$network))
    commAffil_mat <- as.data.frame(rows[rep(seq_len(nrow(rows)), each = 333), ])
    rownames(commAffil_mat) <- NULL
  } else if(str_detect(atlas, "schaefer200x7")){
    parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
    communityAffil <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Schaefer_7Network/schaefer200x7CommunityAffiliation_final.csv")
    parcel.labels <- data.frame(cbind(parcel.labels$label, communityAffil$CommAffil))
    names(parcel.labels) <- c("label", "network")
    rows <- as.data.frame(t(parcel.labels$network))
    commAffil_mat <- as.data.frame(rows[rep(seq_len(nrow(rows)), each = 200), ])
    rownames(commAffil_mat) <- NULL
  } else if(str_detect(atlas, "schaefer200x17")){
    parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
    communityAffil <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Schaefer_17Network/schaefer200x17CommunityAffiliation_final.csv")
    parcel.labels <- data.frame(cbind(parcel.labels$label, communityAffil$CommAffil))
    names(parcel.labels) <- c("label", "network")
    rows <- as.data.frame(t(parcel.labels$network))
    commAffil_mat <- as.data.frame(rows[rep(seq_len(nrow(rows)), each = 200), ])
    rownames(commAffil_mat) <- NULL
  } else if(str_detect(atlas, "schaefer400x7")){
    parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x7_regionlist_final.csv")
    communityAffil <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Schaefer_7Network/schaefer400x7CommunityAffiliation_final.csv")
    parcel.labels <- data.frame(cbind(parcel.labels$label, communityAffil$CommAffil))
    names(parcel.labels) <- c("label", "network")
    rows <- as.data.frame(t(parcel.labels$network))
    commAffil_mat <- as.data.frame(rows[rep(seq_len(nrow(rows)), each = 400), ])
    rownames(commAffil_mat) <- NULL
  } else if(str_detect(atlas, "schaefer400x17")){
    parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x17_regionlist_final.csv")
    communityAffil <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Schaefer_17Network/schaefer400x17CommunityAffiliation_final.csv")
    parcel.labels <- data.frame(cbind(parcel.labels$label, communityAffil$CommAffil))
    names(parcel.labels) <- c("label", "network")
    rows <- as.data.frame(t(parcel.labels$network))
    commAffil_mat <- as.data.frame(rows[rep(seq_len(nrow(rows)), each = 400), ])
    rownames(commAffil_mat) <- NULL
  } else {
    print("Please provide valid atlas gordon, schaefer200x7, schaefer200x17, schaefer400x7, or schaefer400x17")
  } 
  
  
  #compute average between-network connectivity OR average within-network connectivity
  if(metric == "BNC"){
    BNC <- c()
    for (i in 1:nrow(connect.matrix)) {  
      between_network_values <- c() 
      BNC_per_node <- c()
      for (j in 1:length(connect.matrix[i,])) {  # entire  row
        if(is.na(connect.matrix[i,j]) || connect.matrix[i,j] == 1) { # add NA if identity
          between_network_values <- append(between_network_values, NA);
        } else if (parcel.labels$network[i] != commAffil_mat[i,j]) { 
          between_network_values <- append(between_network_values, connect.matrix[i,j])
        } else if (parcel.labels$network[i] == commAffil_mat[i,j]) {  
          between_network_values <- append(between_network_values, NA)
        } else {
          print("Bug in BNC calculation ")
        }
      }
      #print(between_network_values)  
      BNC_per_node <- append(BNC_per_node, mean(between_network_values, na.rm=TRUE))  
      #print(BNC_per_node) # 1 value
      BNC <- append(BNC, BNC_per_node)  
      metric_toReturn <- BNC
    }
  } else if (metric == "WNC") {
    WNC <- c()
    for (i in 1:nrow(connect.matrix)) {  
      within_network_values <- c() 
      WNC_per_node <- c()
      for (j in 1:length(connect.matrix[i,])) {  # entire  row
        if(is.na(connect.matrix[i,j]) || connect.matrix[i,j] == 1) { # add NA if identity
          within_network_values <- append(within_network_values, NA);
        } else if (parcel.labels$network[i] == commAffil_mat[i,j]) { 
          within_network_values <- append(within_network_values, connect.matrix[i,j])
        } else if (parcel.labels$network[i] != commAffil_mat[i,j]) {  
          within_network_values <- append(within_network_values, NA)
        } else {
          print("Bug in WNC calculation ")
        }
      }
      #print(within_network_values)  
      WNC_per_node <- append(WNC_per_node, mean(within_network_values, na.rm=TRUE))  
      #print(WNC_per_node) # 1 value
      WNC <- append(WNC, WNC_per_node)  
      metric_toReturn <- WNC
    }
  } else {
    print("Please provide valid connectivity metric: WNC or BNC")
    metric_toReturn <- NA
  }
  return(metric_toReturn)
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
  } else if (atlas =="schaefer200x7"){
    n <- 201
    regionheaders <- as.character(schaefer200x7.parcel.labels$label)
  } else if (atlas =="schaefer200x17"| atlas=="schaefer200") {
    n <- 201
    regionheaders <- as.character(schaefer200x17.parcel.labels$label)
  } else if (atlas =="schaefer400x7") {
    n <- 401
    regionheaders <- as.character(schaefer400x7.parcel.labels$label)
  } else if (atlas =="schaefer400x17"| atlas=="schaefer400") {
    n <- 401
    regionheaders <- as.character(schaefer400x17.parcel.labels$label)
  }
  subxparcel.matrix  <- matrix(data = NA, nrow = length(subject_list), ncol = n)
  demoheaders <- c("subject")
  colheaders <- as.matrix(c(demoheaders,regionheaders))
  colnames(subxparcel.matrix) <- colheaders
  return(subxparcel.matrix)
  
}
 
  
# Function for Extracting Parcel-Parcel Connectivity
# @param subject rbcid of subject of interest
# @param atlas A character string, name of atlas of interest (e.g. glasser, gordon, schaefer200x7, schaefer200x17, schaefer400x7, or schaefer400x17)
# @param dataset A character string, name of dataset
extractParcel2ParcelConn <- function(subject, atlas, dataset){
  
  #read in connectivity matrix  # need to fix to make more efficient  
  if(dataset == "PNC" | dataset=="HCPD" | dataset=="HBN") {
    connect.matrix <- readRDS(sprintf("/cbica/projects/network_replication/input/%1$s/connMatricesData/connectivity_matrices/%2$s_ConnMatrices.RData", dataset, subject))
    if(str_detect(atlas, "schaefer200") | str_detect(atlas, "schaefer400")) {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
  } else if(dataset == "NKI") { #check if BAS1 exists for subject
    connect.matrix <- readRDS(sprintf("/cbica/projects/network_replication/input/NKI/connMatricesData/connectivity_matrices/%1$s_ConnMatrices.RData", subject))
    ses_name <- str_extract(names(connect.matrix), "[A-Z]{3}1")[1]
    if(str_detect(atlas, "schaefer")) {
      atlas_name <- str_extract(atlas, "schaefer[0-9]")
      atlas_name <- paste0(atlas_name, "17")
    } else {
      atlas_name <- atlas
    }
    matrixName <- paste0(ses_name, "_", atlas_name, "_conn")
    connect.matrix <- connect.matrix[[matrixName]]
  }  
    
  
  if (atlas == "glasser"){
    edge <- read.csv("/cbica/projects/network_replication/atlases/edge/glasser_edge.csv")
  } else if (atlas == "gordon"){
    edge <- read.csv("/cbica/projects/network_replication/atlases/edge/gordon_edge.csv")
  }  else if(atlas == "schaefer200x7"){
    edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer200x7_edge.csv")
  } else if(atlas == "schaefer400x7"){
    edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer400x7_edge.csv")
  } else if(atlas == "schaefer200x17"){
    edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer200x17_edge.csv")
  } else if(atlas == "schaefer400x17"){
    edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer400x17_edge.csv")
  } else {
    print("Please provide valid atlas (glasser, gordon, schaefer200, or schaefer400")
  } 
   
  # set upper triangle of matrix to 9999 and vectorize the matrix
  mat <- connect.matrix
  mat[upper.tri(mat)] <- 9999
  connect.vector <- c(mat)
  connect.vector <- connect.vector[-c(which(connect.vector==9999))] # this doesn't remove the diagonal
  connect.vector <- connect.vector[-c(which(connect.vector==1))] # remove edges that connect the same parcel to itself 
  print(paste(which(participants %in% subject), "/", length(participants), "-", atlas))
  return(connect.vector)
  
}






 






 
 