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
  participants <- read.csv("/cbica/projects/network_replication/input/PNC/sample_selection/PNC_demographics_finalsample_20230103.csv")
  participants <- paste0("sub-", participants$rbcid)
} else if (dataset == "HBN") {
  participants <- read.csv("/cbica/projects/network_replication/input/HBN/sample_selection/HBN_demographics_finalsample_202230226.csv")
  participants <- participants$sub
} else {
  participants <- NA
  print("Provide valid dataset")
}



# make output dfs: 
## WNC.subxparcel.matrix.gordon, WNC.subxparcel.matrix.schaefer200x7, WNC.subxparcel.matrix.schaefer400x7
atlases <- c("gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", "schaefer400x17")
#atlases <- c("schaefer200x7", "schaefer200x17", "schaefer400x7", "schaefer400x17")

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

 
write.csv(WNC.subxparcel.matrix.gordon, sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNC_subxparcel_matrix_gordon.csv", dataset), row.names=F, quote=F)
write.csv(WNC.subxparcel.matrix.schaefer200x7, sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNC_subxparcel_matrix_schaefer200x7.csv", dataset), row.names=F, quote=F)
write.csv(WNC.subxparcel.matrix.schaefer400x7, sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNC_subxparcel_matrix_schaefer400x7.csv", dataset), row.names=F, quote=F)
write.csv(WNC.subxparcel.matrix.schaefer200x17, sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNC_subxparcel_matrix_schaefer200x17.csv", dataset), row.names=F, quote=F)
write.csv(WNC.subxparcel.matrix.schaefer400x17, sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNC_subxparcel_matrix_schaefer400x17.csv", dataset), row.names=F, quote=F)