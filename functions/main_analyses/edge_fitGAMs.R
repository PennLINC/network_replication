library(dplyr)
library(janitor)
library(purrr)
library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
library(magrittr)

source("/cbica/projects/network_replication/adapted_Rscripts/GAM_functions.R")
source("/cbica/projects/network_replication/adapted_Rscripts/GAM_functions_Val.R")
source("/cbica/projects/network_replication/adapted_Rscripts/fitGAMs.R")


args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
print(dataset)

if(dataset=="PNC"){
  
  demographics <- read.csv("/cbica/projects/network_replication/input/PNC/sample_selection/PNC_demographics_finalsample_20230103.csv")
  demographics$subject <- paste0("sub-", demographics$rbcid)
} else if(dataset=="HCPD") {
  demographics <- read.csv(sprintf("/cbica/projects/network_replication/input/%1$s/sample_selection/%1$s_demographics_finalsample_20221226.csv", dataset))
  names(demographics)[c(which(names(demographics) == "sub"))] <- "subject"
  demographics$age <- demographics$interview_age/12
} else if(dataset=="NKI") {
  demographics <- read.csv(sprintf("/cbica/projects/network_replication/input/%1$s/sample_selection/%1$s_demographics_finalsample_20221219.csv", dataset))
  demographics$sex <- demographics$gender
} else if(dataset=="HBN") {
  demographics <- read.csv(sprintf("/cbica/projects/network_replication/input/%1$s/sample_selection/%1$s_demographics_finalsample_202230226.csv", dataset))
  demographics <- demographics[-c(which(demographics$sub=="sub-NDARCE721YB5")),]
  names(demographics)[c(which(names(demographics) == "sub"))] <- "subject"
} else {
  print("Provide valid dataset")
}
demographics <- demographics[,c(which(names(demographics) == "subject"), which(names(demographics) == "sex"), which(names(demographics)=="age"), which(names(demographics) =="meanFD_avgSes"))] 
  
  
  
if (dataset=="PNC" | dataset =="NKI") {
  # edge 
  ## edge.glasser.dataset, edge.gordon.dataset, edge.schaefer200x7.dataset, edge.schaefer200x17.dataset, edge.schaefer400x7.dataset, edge.schaefer400x17.dataset
  atlases <- c("glasser","gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", "schaefer400x17")
  #make_EdgeDemog("gordon", demographics, "subject", dataset)
  for(i in c(1:length(atlases))){
    df_name <- paste0("edge.", atlases[i], ".", dataset)
    assign(df_name, make_EdgeDemog(atlases[i],  demographics, "subject", dataset))
  } 
  print("make_EdgeDemog complete")
} else if (dataset=="HBN" | dataset=="HCPD") {
  atlases <- c("schaefer200x7")
  #atlases <- c("glasser","gordon", "schaefer200x7", "schaefer200x17", "schaefer400x7", "schaefer400x17")
  for(i in c(1:length(atlases))){
    df_name <- paste0("edge.", atlases[i], ".", dataset)
    assign(df_name, make_EdgeDemog_covbat(atlases[i],  demographics, "subject", dataset))
  } 
  print("make_EdgeDemog_covbat complete")
}


# load edge labels
gordon_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/gordon_edge.csv")
gordon_edge <- gordon_edge[,2]

glasser_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/glasser_edge.csv")
glasser_edge <- glasser_edge[,2]
schaefer200x7_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer200x7_edge.csv")
schaefer200x7_edge <- schaefer200x7_edge[,2]
schaefer200x17_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer200x17_edge.csv")
schaefer200x17_edge <- schaefer200x17_edge[,2]

schaefer400x7_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer400x7_edge.csv")
schaefer400x7_edge <- schaefer400x7_edge[,2]
schaefer400x17_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer400x17_edge.csv")
schaefer400x17_edge <- schaefer400x17_edge[,2]



glasser_edge <- gsub("-", ".", glasser_edge)

schaefer200x7_edge <- gsub("7Networks", "Networks", schaefer200x7_edge)
schaefer200x17_edge <- gsub("17Networks", "Networks", schaefer200x17_edge)

schaefer400x7_edge <- gsub("7Networks", "Networks", schaefer400x7_edge)

schaefer400x17_edge <- gsub("17Networks", "Networks", schaefer400x17_edge)

  



if(dataset=="PNC"){
  
  edge.schaefer200x7.PNC <- readRDS("/cbica/projects/network_replication/output/PNC/edge/schaefer200x7_demographics_finalsample.RData")
  edge.schaefer200x7.PNC$sex <- as.factor(edge.schaefer200x7.PNC$sex)
  names(edge.schaefer200x7.PNC) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.PNC))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.PNC)
  
  
  edge.gordon.PNC <- readRDS("/cbica/projects/network_replication/output/PNC/edge/gordon_demographics_finalsample.RData")
  edge.gordon.PNC$sex <- as.factor(edge.gordon.PNC$sex)
  gam.edge.age.gordon <- fitGAMs_edge(gordon_edge,"edge","gordon", dataset)
  rm(gam.edge.age.gordon)
  rm(edge.gordon.PNC)
  
  edge.glasser.PNC <- readRDS("/cbica/projects/network_replication/output/PNC/edge/glasser_demographics_finalsample.RData")
  edge.glasser.PNC$sex <- as.factor(edge.glasser.PNC$sex)
  names(edge.glasser.PNC) <- gsub("-", ".",names(edge.glasser.PNC))
  gam.edge.age.glasser <- fitGAMs_edge(glasser_edge,"edge","glasser", dataset) 
  rm(gam.edge.age.glasser)
  rm(edge.glasser.PNC)
  
  
  
  edge.schaefer200x17.PNC <- readRDS("/cbica/projects/network_replication/output/PNC/edge/schaefer200x17_demographics_finalsample.RData")
  edge.schaefer200x17.PNC$sex <- as.factor(edge.schaefer200x17.PNC$sex)
  names(edge.schaefer200x17.PNC) <- gsub("17Networks", "Networks", names(edge.schaefer200x17.PNC))
  gam.edge.age.schaefer200x17 <- fitGAMs_edge(schaefer200x17_edge,"edge","schaefer200x17", dataset)
  rm(gam.edge.age.schaefer200x17)
  rm(edge.schaefer200x17.PNC)
  
  edge.schaefer400x7.PNC <- readRDS("/cbica/projects/network_replication/output/PNC/edge/schaefer400x7_demographics_finalsample.RData")
  edge.schaefer400x7.PNC$sex <- as.factor(edge.schaefer400x7.PNC$sex)
  names(edge.schaefer400x7.PNC) <- gsub("7Networks", "Networks", names(edge.schaefer400x7.PNC))
  gam.edge.age.schaefer400x7 <- fitGAMs_edge(schaefer400x7_edge,"edge","schaefer400x7", dataset)
  rm(gam.edge.age.schaefer400x7)
  rm(edge.schaefer400x7.PNC)
  
  edge.schaefer400x17.PNC <- readRDS("/cbica/projects/network_replication/output/PNC/edge/schaefer400x17_demographics_finalsample.RData")
  edge.schaefer400x17.PNC$sex <- as.factor(edge.schaefer400x17.PNC$sex)
  names(edge.schaefer400x17.PNC) <- gsub("17Networks", "Networks", names(edge.schaefer400x17.PNC))
  gam.edge.age.schaefer400x17 <- fitGAMs_edge(schaefer400x17_edge,"edge","schaefer400x17", dataset)
  rm(gam.edge.age.schaefer400x17)
  rm(edge.schaefer400x17.PNC)
  
   
  
} else if(dataset=="NKI") {
  edge.schaefer200x7.NKI <- readRDS("/cbica/projects/network_replication/output/NKI/edge/schaefer200x7_demographics_finalsample.RData")
  edge.schaefer200x7.NKI$sex <- as.factor(edge.schaefer200x7.NKI$sex)
  names(edge.schaefer200x7.NKI) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.NKI))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.NKI)
  
  
  edge.gordon.NKI <- readRDS("/cbica/projects/network_replication/output/NKI/edge/gordon_demographics_finalsample.RData")
  edge.gordon.NKI$sex <- as.factor(edge.gordon.NKI$sex)
  gam.edge.age.gordon <- fitGAMs_edge(gordon_edge,"edge","gordon", dataset)
  rm(gam.edge.age.gordon)
  rm(edge.gordon.NKI)
  
  
  edge.glasser.NKI <- readRDS("/cbica/projects/network_replication/output/NKI/edge/glasser_demographics_finalsample.RData")
  edge.glasser.NKI$sex <- as.factor(edge.glasser.NKI$sex)
  names(edge.glasser.NKI) <- gsub("-", ".",names(edge.glasser.NKI))
  gam.edge.age.glasser <- fitGAMs_edge(glasser_edge,"edge","glasser", dataset) 
  rm(gam.edge.age.glasser)
  rm(edge.glasser.NKI)

  
  
  
  edge.schaefer200x17.NKI <- readRDS("/cbica/projects/network_replication/output/NKI/edge/schaefer200x17_demographics_finalsample.RData")
  edge.schaefer200x17.NKI$sex <- as.factor(edge.schaefer200x17.NKI$sex)
  names(edge.schaefer200x17.NKI) <- gsub("17Networks", "Networks", names(edge.schaefer200x17.NKI))
  gam.edge.age.schaefer200x17 <- fitGAMs_edge(schaefer200x17_edge,"edge","schaefer200x17", dataset)
  rm(gam.edge.age.schaefer200x17)
  rm(edge.schaefer200x17.NKI)
  
  edge.schaefer400x7.NKI <- readRDS("/cbica/projects/network_replication/output/NKI/edge/schaefer400x7_demographics_finalsample.RData")
  edge.schaefer400x7.NKI$sex <- as.factor(edge.schaefer400x7.NKI$sex)
  names(edge.schaefer400x7.NKI) <- gsub("7Networks", "Networks", names(edge.schaefer400x7.NKI))
  gam.edge.age.schaefer400x7 <- fitGAMs_edge(schaefer400x7_edge,"edge","schaefer400x7", dataset)
  rm(gam.edge.age.schaefer400x7)
  rm(edge.schaefer400x7.NKI)
  
  edge.schaefer400x17.NKI <- readRDS("/cbica/projects/network_replication/output/NKI/edge/schaefer400x17_demographics_finalsample.RData")
  edge.schaefer400x17.NKI$sex <- as.factor(edge.schaefer400x17.NKI$sex)
  names(edge.schaefer400x17.NKI) <- gsub("17Networks", "Networks", names(edge.schaefer400x17.NKI))
  gam.edge.age.schaefer400x17 <- fitGAMs_edge(schaefer400x17_edge,"edge","schaefer400x17", dataset)
  rm(gam.edge.age.schaefer400x17)
  rm(edge.schaefer400x17.NKI)
  
  
   
} else if(dataset=="HCPD") {
  edge.schaefer200x7.HCPD <- readRDS("/cbica/projects/network_replication/output/HCPD/edge/schaefer200x7_demographics_finalsample_covbat.RData")
  edge.schaefer200x7.HCPD$sex <- as.factor(edge.schaefer200x7.HCPD$sex)
  names(edge.schaefer200x7.HCPD) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.HCPD))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge_covbat(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.HCPD)
  
  edge.gordon.HCPD <- readRDS("/cbica/projects/network_replication/output/HCPD/edge/gordon_demographics_finalsample_covbat.RData")
  edge.gordon.HCPD$sex <- as.factor(edge.gordon.HCPD$sex)
  gam.edge.age.gordon <- fitGAMs_edge_covbat(gordon_edge,"edge","gordon", dataset)
  rm(gam.edge.age.gordon)
  rm(edge.gordon.HCPD)
  
  edge.glasser.HCPD <- readRDS("/cbica/projects/network_replication/output/HCPD/edge/glasser_demographics_finalsample_covbat.RData")
  edge.glasser.HCPD$sex <- as.factor(edge.glasser.HCPD$sex)
  names(edge.glasser.HCPD) <- gsub("-", ".",names(edge.glasser.HCPD))
  gam.edge.age.glasser <- fitGAMs_edge_covbat(glasser_edge,"edge","glasser", dataset) 
  rm(gam.edge.age.glasser)
  rm(edge.glasser.HCPD)
  
  
  
  edge.schaefer200x17.HCPD <- readRDS("/cbica/projects/network_replication/output/HCPD/edge/schaefer200x17_demographics_finalsample_covbat.RData")
  edge.schaefer200x17.HCPD$sex <- as.factor(edge.schaefer200x17.HCPD$sex)
  names(edge.schaefer200x17.HCPD) <- gsub("17Networks", "Networks", names(edge.schaefer200x17.HCPD))
  gam.edge.age.schaefer200x17 <- fitGAMs_edge_covbat(schaefer200x17_edge,"edge","schaefer200x17", dataset)
  rm(gam.edge.age.schaefer200x17)
  rm(edge.schaefer200x17.HCPD)
  
  edge.schaefer400x7.HCPD <- readRDS("/cbica/projects/network_replication/output/HCPD/edge/schaefer400x7_demographics_finalsample_covbat.RData")
  edge.schaefer400x7.HCPD$sex <- as.factor(edge.schaefer400x7.HCPD$sex)
  names(edge.schaefer400x7.HCPD) <- gsub("7Networks", "Networks", names(edge.schaefer400x7.HCPD))
  gam.edge.age.schaefer400x7 <- fitGAMs_edge_covbat(schaefer400x7_edge,"edge","schaefer400x7", dataset)
  rm(gam.edge.age.schaefer400x7)
  rm(edge.schaefer400x7.HCPD)
  
  edge.schaefer400x17.HCPD <- readRDS("/cbica/projects/network_replication/output/HCPD/edge/schaefer400x17_demographics_finalsample_covbat.RData")
  edge.schaefer400x17.HCPD$sex <- as.factor(edge.schaefer400x17.HCPD$sex)
  names(edge.schaefer400x17.HCPD) <- gsub("17Networks", "Networks", names(edge.schaefer400x17.HCPD))
  gam.edge.age.schaefer400x17 <- fitGAMs_edge_covbat(schaefer400x17_edge,"edge","schaefer400x17", dataset)
  rm(gam.edge.age.schaefer400x17)
  rm(edge.schaefer400x17.HCPD)
  
  
} else if(dataset=="HBN") {
  edge.schaefer200x7.HBN <- readRDS("/cbica/projects/network_replication/output/HBN/edge/schaefer200x7_demographics_finalsample_covbat.RData")
  edge.schaefer200x7.HBN$sex <- as.factor(edge.schaefer200x7.HBN$sex)
  names(edge.schaefer200x7.HBN) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.HBN))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge_covbat(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.HBN)
  
  #edge.gordon.HBN <- readRDS("/cbica/projects/network_replication/output/HBN/edge/gordon_demographics_finalsample_covbat.RData")
  #edge.gordon.HBN$sex <- as.factor(edge.gordon.HBN$sex)
  #gam.edge.age.gordon <- fitGAMs_edge_covbat(gordon_edge,"edge","gordon", dataset)
  #rm(gam.edge.age.gordon)
  #rm(edge.gordon.HBN)
  
  #edge.glasser.HBN <- readRDS("/cbica/projects/network_replication/output/HBN/edge/glasser_demographics_finalsample_covbat.RData")
  #edge.glasser.HBN$sex <- as.factor(edge.glasser.HBN$sex)
  #names(edge.glasser.HBN) <- gsub("-", ".",names(edge.glasser.HBN))
  #gam.edge.age.glasser <- fitGAMs_edge_covbat(glasser_edge,"edge","glasser", dataset) 
  #rm(gam.edge.age.glasser)
  #rm(edge.glasser.HBN)
  
  
  
  #edge.schaefer200x17.HBN <- readRDS("/cbica/projects/network_replication/output/HBN/edge/schaefer200x17_demographics_finalsample_covbat.RData")
  #edge.schaefer200x17.HBN$sex <- as.factor(edge.schaefer200x17.HBN$sex)
  #names(edge.schaefer200x17.HBN) <- gsub("17Networks", "Networks", names(edge.schaefer200x17.HBN))
  #gam.edge.age.schaefer200x17 <- fitGAMs_edge_covbat(schaefer200x17_edge,"edge","schaefer200x17", dataset)
  #rm(gam.edge.age.schaefer200x17)
  #rm(edge.schaefer200x17.HBN)
  
  #edge.schaefer400x7.HBN <- readRDS("/cbica/projects/network_replication/output/HBN/edge/schaefer400x7_demographics_finalsample_covbat.RData")
  #edge.schaefer400x7.HBN$sex <- as.factor(edge.schaefer400x7.HBN$sex)
  #names(edge.schaefer400x7.HBN) <- gsub("7Networks", "Networks", names(edge.schaefer400x7.HBN))
  #gam.edge.age.schaefer400x7 <- fitGAMs_edge_covbat(schaefer400x7_edge,"edge","schaefer400x7", dataset)
  #rm(gam.edge.age.schaefer400x7)
  #rm(edge.schaefer400x7.HBN)
  
  #edge.schaefer400x17.HBN <- readRDS("/cbica/projects/network_replication/output/HBN/edge/schaefer400x17_demographics_finalsample_covbat.RData")
  #edge.schaefer400x17.HBN$sex <- as.factor(edge.schaefer400x17.HBN$sex)
  #names(edge.schaefer400x17.HBN) <- gsub("17Networks", "Networks", names(edge.schaefer400x17.HBN))
  #gam.edge.age.schaefer400x17 <- fitGAMs_edge_covbat(schaefer400x17_edge,"edge","schaefer400x17", dataset)
  #rm(gam.edge.age.schaefer400x17)
  #rm(edge.schaefer400x17.HBN)
  
} else {
  print("Provide valid dataset")
}






  