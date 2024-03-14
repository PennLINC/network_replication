source("/cbica/projects/network_replication/manuscript/code/scripts/sensitivity_analyses/computeConnectivityMetrics_restOnly.R")

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]

################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/manuscript/code/config_%1$s.json", dataset))
data_root <- config_data$data_root
conn_matrices_dir <- paste0(data_root, "connMatricesData/connectivity_matrices_restOnly/")
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
# Read files 
################## 
participants <- read.table(sprintf("%1$s%2$s_final_subjectlist_restOnly.txt", sample_selection_dir, dataset))
if (dataset=="NKI") {
  participants <- gsub("sub-", "", participants$V1)
} else {
  participants <- participants$V1
}



# make output dfs: 
## GBC.subxparcel.matrix.schaefer200 
atlases <- "schaefer200" # only include schaefer200 for sensitivity analysis with resting-only data in supplement
metric = "GBC"

for(i in c(1:length(atlases))){
  df_name <- paste0(metric, ".subxparcel.matrix.", atlases[i])
  assign(df_name, make_output_dfs(participants, atlases[i]))
}


# compute GBC
for(sub in c(1:length(participants))){
  for(i in c(1:length(atlases))){
    subjectID=as.character(participants[sub])
    id.data <- computeGBC_restOnly(subjectID, atlases[i], dataset)
    dataname <- sprintf("%s.%s.%s", metric, "subxparcel.matrix", atlases[i]) 
    df_toUpdate <- get(dataname)
    df_toUpdate[sub,] <- cbind(subjectID, t(id.data)) #update the name of the df every iteration to GBC.subxparcel.matrix.[atlas]
    assign(dataname, df_toUpdate)  
    print(paste(sub, "/", length(participants), "-", subjectID, atlases[i]))
  }
}

 
   
write.csv(GBC.subxparcel.matrix.schaefer200, sprintf("%1$sGBC_subxparcel_matrix_schaefer200.csv", metric_output_dir), row.names=F, quote=F)
 