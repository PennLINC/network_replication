source("/cbica/projects/network_replication/manuscript/code/scripts/main_analyses/computeConnectivityMetrics.R")

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]


################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/manuscript/code/config_%1$s.json", dataset))
outputs_root <- config_data$outputs_root
sample_selection_dir <- config_data$sample_selection_data_dir
metric_output_dir <- paste0(outputs_root, "edge/")

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
# Read files 
################## 
participants <- read.table(sprintf("%1$s%2$s_final_subjectlist.txt", sample_selection_dir, dataset))
if (dataset=="NKI") {
  participants <- gsub("sub-", "", participants$V1)
} else {
  participants <- participants$V1
}



################## 
# Compute metric 
################## 
# make a dataframe in which rows are subjects, and columns are the names of the edges (parcel1_to_parcel2)
# the entries represent the parcel-by-parcel correlation
 
# for each subject, turn the connectivity matrix into a vector

#extractEdge.gordon <- lapply(participants, extractParcel2ParcelConn, "gordon", dataset)
#extractEdge.glasser <- lapply(participants, extractParcel2ParcelConn, "glasser", dataset)
extractEdge.schaefer200x7 <- lapply(participants, extractParcel2ParcelConn, "schaefer200x7", dataset)
#extractEdge.schaefer200x17 <- lapply(participants, extractParcel2ParcelConn, "schaefer200x17", dataset)
#extractEdge.schaefer400x7 <- lapply(participants, extractParcel2ParcelConn, "schaefer400x7", dataset)
#extractEdge.schaefer400x17 <- lapply(participants, extractParcel2ParcelConn, "schaefer400x17", dataset)



atlases <- c("schaefer200x7") # not doing other atlases in manuscript
 
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
  
  saveRDS(subxedge, sprintf("%1$ssubxedge_%2$s.RData", metric_output_dir, atlases[i])) 
  print(atlases[i])
} 