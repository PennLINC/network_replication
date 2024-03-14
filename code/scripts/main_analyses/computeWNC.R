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
metric_output_dir <- paste0(outputs_root, "WNC/")

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
# make output dfs: 
##  WNC.subxparcel.matrix.schaefer200x7 
atlases <- c("schaefer200x7")
metric = "WNC"

for(i in c(1:length(atlases))){
  df_name <- paste0(metric, ".subxparcel.matrix.", atlases[i])
  assign(df_name, make_output_dfs(participants, atlases[i]))
}


# compute WNC
for(sub in c(1:length(participants))){
  for(i in c(1:length(atlases))){
    subjectID=as.character(participants[sub])
    print(paste(subjectID, atlases[i], metric, dataset))
    id.data <- computeBNC_WNC(subjectID, atlases[i], metric, dataset)
    
    dataname <- sprintf("%s.%s.%s", metric, "subxparcel.matrix", atlases[i]) 
    df_toUpdate <- get(dataname)
    df_toUpdate[sub,] <- cbind(subjectID, t(id.data)) #update the name of the df every iteration to WNC.subxparcel.matrix.[atlas]
    assign(dataname, df_toUpdate)  
    print(paste(sub, "/", length(participants), "-", subjectID, atlases[i]))
  }
}


#write.csv(WNC.subxparcel.matrix.gordon, sprintf("%1$sWNC_subxparcel_matrix_gordon.csv", metric_output_dir), row.names=F, quote=F)
write.csv(WNC.subxparcel.matrix.schaefer200x7, sprintf("%1$sWNC_subxparcel_matrix_schaefer200x7.csv", metric_output_dir), row.names=F, quote=F)
#write.csv(WNC.subxparcel.matrix.schaefer400x7, sprintf("%1$sWNC_subxparcel_matrix_schaefer400x7.csv", metric_output_dir), row.names=F, quote=F)
#write.csv(WNC.subxparcel.matrix.schaefer200x17, sprintf("%1$sWNC_subxparcel_matrix_schaefer200x17.csv", metric_output_dir), row.names=F, quote=F)
#write.csv(WNC.subxparcel.matrix.schaefer400x17, sprintf("%1$sWNC_subxparcel_matrix_schaefer400x17.csv", metric_output_dir), row.names=F, quote=F)