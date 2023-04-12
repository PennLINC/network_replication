# perform covbat-gam on edges

library(ComBatFamily)
library(dplyr)
library(mgcv)

args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]



# function for performing CovBat on connectivity metrics output for all atlases
apply_covbat_edge <- function(atlas, dataset) {
  subxedge <- readRDS(sprintf("/cbica/projects/network_replication/output/%1$s/edge/subxedge_%2$s.RData", dataset, atlas))
  print("Edge RData loaded")
  subjects <- subxedge$subject
  subxedge <- data.frame(subxedge[,-1])
  row.names(subxedge) <- subjects
  print("subxedge loaded")
  if(dataset=="HCPD") {
    participants <- read.csv("/cbica/projects/network_replication/input/HCPD/sample_selection/HCPD_demographics_finalsample_20221226.csv")
    participants <- participants[match(rownames(subxedge), participants$sub),]
    age_vec <- participants$interview_age/12
    sex_vec <- as.factor(participants$sex)
    meanFD_avgSes_vec <- participants$meanFD_avgSes
    batch <- participants$site
    
    covar_df <- bind_cols(participants$sub, as.numeric(age_vec), as.factor(sex_vec), as.numeric(meanFD_avgSes_vec))
    covar_df <- dplyr::rename(covar_df, sub=...1,
                              age = ...2,
                              sex = ...3,
                              meanFD_avgSes = ...4)
    print("Harmonizing data")
    rm(participants)
    data.harmonized <- covfam(data=subxedge, bat = as.factor(batch), covar = covar_df, gam, y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(meanFD_avgSes))
  } else if (dataset=="HBN") { 
    participants <- read.csv("/cbica/projects/network_replication/input/HBN/sample_selection/HBN_demographics_finalsample_202230223.csv")
    participants_noAge_index <- c(which(is.na(participants$age))) #age not collected on 58 participants 
    participants <- participants[-c(which(is.na(participants$age))),] 
    participants <- participants[-c(which(participants$sub=="sub-NDARCE721YB5")),]
    age_vec <- participants$age 
    sex_vec <- as.factor(participants$sex)
    meanFD_avgSes_vec <- participants$meanFD_avgSes
    batch <- participants$ses
    
    covar_df <- bind_cols(participants$sub, as.numeric(age_vec), as.factor(sex_vec), as.numeric(meanFD_avgSes_vec))
    covar_df <- dplyr::rename(covar_df, sub=...1,
                              age = ...2,
                              sex = ...3,
                              meanFD_avgSes = ...4)
    
    print("Harmonizing data")
    rm(participants)
    data.harmonized <- covfam(data=subxedge, bat = as.factor(batch), covar = covar_df, gam, y ~ s(age, k=3, fx=T) + as.factor(sex) + as.numeric(meanFD_avgSes))
    
  } else {
    print("provide proper participants df")
  }
  data.harmonized_covbat <- data.frame(data.harmonized$dat.covbat)
  print(paste(atlas, "edge", dataset, "done"))
  saveRDS(data.harmonized_covbat, sprintf("/cbica/projects/network_replication/output/%1$s/edge/subxedge_%2$s_covbat.RData", dataset, atlas))
  rm(data.harmonized_covbat)
  gc()
  return(covar_df)
}
 

atlases <- c("gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", "schaefer400x17")

#edge_harmonized <- lapply(atlases, apply_covbat_edge, dataset)
apply_covbat_edge("gordon", dataset)
apply_covbat_edge("glasser", dataset)
apply_covbat_edge("schaefer200x7", dataset)
apply_covbat_edge("schaefer200x17", dataset)
apply_covbat_edge("schaefer400x7", dataset)
apply_covbat_edge("schaefer400x17", dataset)
#names(edge_harmonized) <- atlases

#for(i in c(1:length(edge_harmonized))){
#  saveRDS(edge_harmonized[[atlases[i]]], sprintf("/cbica/projects/network_replication/output/%1$s/edge/subxedge_%2$s_covbat.RData", dataset, atlases[i]))
#  print(atlases[i])
#} 