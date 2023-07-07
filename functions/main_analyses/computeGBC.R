source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/computeConnectivityMetrics.R")

args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]


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




# make output dfs: 
## GBC.subxparcel.matrix.glasser, GBC.subxparcel.matrix.gordon, GBC.subxparcel.matrix.schaefer200x7, GBC.subxparcel.matrix.schaefer400x7
atlases <- c("glasser", "gordon", "schaefer200", "schaefer400")
metric = "GBC"

for(i in c(1:length(atlases))){
  df_name <- paste0(metric, ".subxparcel.matrix.", atlases[i])
  assign(df_name, make_output_dfs(participants, atlases[i]))
}


# compute GBC
for(sub in c(1:length(participants))){
  for(i in c(1:length(atlases))){
    subjectID=as.character(participants[sub])
    id.data <- computeGBC(subjectID, atlases[i], dataset)
    dataname <- sprintf("%s.%s.%s", metric, "subxparcel.matrix", atlases[i]) 
    df_toUpdate <- get(dataname)
    df_toUpdate[sub,] <- cbind(subjectID, t(id.data)) #update the name of the df every iteration to GBC.subxparcel.matrix.[atlas]
    assign(dataname, df_toUpdate)  
    print(paste(sub, "/", length(participants), "-", subjectID, atlases[i]))
  }
}


write.csv(GBC.subxparcel.matrix.glasser, sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBC_subxparcel_matrix_glasser.csv", dataset), row.names=F, quote=F)
write.csv(GBC.subxparcel.matrix.gordon, sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBC_subxparcel_matrix_gordon.csv", dataset), row.names=F, quote=F)
write.csv(GBC.subxparcel.matrix.schaefer200, sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBC_subxparcel_matrix_schaefer200.csv", dataset), row.names=F, quote=F)
write.csv(GBC.subxparcel.matrix.schaefer400, sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBC_subxparcel_matrix_schaefer400.csv", dataset), row.names=F, quote=F)
