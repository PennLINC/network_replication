source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/computeConnectivityMetrics.R")


args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]



if (dataset == "NKI"){
  participants <- read.csv("/cbica/projects/network_replication/input/NKI/sample_selection/NKI_demographics_finalsample_20230629.csv")
  participants <- participants$subject 
} else if (dataset=="HCPD") {
  participants <- read.csv("/cbica/projects/network_replication/input/HCPD/sample_selection/HCPD_demographics_finalsample_20221226.csv")
  participants <- gsub("HCD", "sub-", participants$src_subject_id)
  
} else if (dataset=="PNC") {
  participants <- read.csv("/cbica/projects/network_replication/input/PNC/sample_selection/PNC_demographics_finalsample_20230629.csv")
  participants <- paste0("sub-", participants$rbcid)
} else if (dataset == "HBN") {
  participants <- read.csv("/cbica/projects/network_replication/input/HBN/sample_selection/HBN_demographics_finalsample_20230629.csv")
  participants <- participants$sub
} else {
  participants <- NA
  print("Provide valid dataset")
}






# make a dataframe in which rows are subjects, and columns are the names of the edges (parcel1_to_parcel2)
# the entries represent the parcel-by-parcel correlation
 
# for each subject, turn the connectivity matrix into a vector

extractEdge.gordon <- lapply(participants, extractParcel2ParcelConn, "gordon", dataset)
extractEdge.glasser <- lapply(participants, extractParcel2ParcelConn, "glasser", dataset)
extractEdge.schaefer200x7 <- lapply(participants, extractParcel2ParcelConn, "schaefer200x7", dataset)
extractEdge.schaefer200x17 <- lapply(participants, extractParcel2ParcelConn, "schaefer200x17", dataset)
extractEdge.schaefer400x7 <- lapply(participants, extractParcel2ParcelConn, "schaefer400x7", dataset)
extractEdge.schaefer400x17 <- lapply(participants, extractParcel2ParcelConn, "schaefer400x17", dataset)



atlases <- c("glasser", "gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", "schaefer400x17")
 
for(i in c(1:length(atlases))){
  extractEdge <- get(paste0("extractEdge.", atlases[i]))
  names(extractEdge) <- participants 
  columns_extractEdge <- bind_rows(extractEdge) # turn list into a dataframe
  subxedge <- as.data.frame(t(columns_extractEdge)) # transpose dataframe
  
  edge <- read.csv(sprintf("/cbica/projects/network_replication/atlases/edge/%1$s_edge.csv", atlases[i]))
  edge <- edge[,2]
  names(subxedge) <- edge # columns = names of edges
  
  subxedge <- subxedge %>% mutate(subject = rownames(subxedge))
  subxedge <- subxedge %>% relocate(subject)
  rownames(subxedge) <- NULL
  
  saveRDS(subxedge, sprintf("/cbica/projects/network_replication/output/%1$s/edge/subxedge_%2$s.RData", dataset, atlases[i])) 
  print(atlases[i])
} 