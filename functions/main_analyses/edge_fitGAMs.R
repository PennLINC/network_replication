library(dplyr)
library(janitor)
library(purrr)
library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
library(magrittr)

source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/GAM_functions.R")
source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/fitGAMs.R")


args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
print(dataset)

# note that subject IDs have slightly different formatting in each dataset, which is why there are different manipulations to the "subject" column
if(dataset=="PNC"){
  demographics <- read.csv("/cbica/projects/network_replication/input/PNC/sample_selection/PNC_demographics_finalsample_20230629.csv")
  names(demographics)[c(which(names(demographics) == "sub"))] <- "subject"
} else if(dataset=="HCPD") {
  demographics <- read.csv(sprintf("/cbica/projects/network_replication/input/%1$s/sample_selection/%1$s_demographics_finalsample_20221226.csv", dataset))
  names(demographics)[c(which(names(demographics) == "sub"))] <- "subject"
  demographics$age <- demographics$interview_age/12
} else if(dataset=="NKI") {
  demographics <- read.csv(sprintf("/cbica/projects/network_replication/input/%1$s/sample_selection/%1$s_demographics_finalsample_20230629.csv", dataset))
  names(demographics)[c(which(names(demographics) == "sub"))] <- "subject"
  demographics$subject <- gsub("sub-", "", demographics$subject)
} else if(dataset=="HBN") {
  demographics <- read.csv(sprintf("/cbica/projects/network_replication/input/%1$s/sample_selection/%1$s_demographics_finalsample_20230629.csv", dataset))
  names(demographics)[c(which(names(demographics) == "sub"))] <- "subject"
} else {
  print("Provide valid dataset")
}
demographics <- demographics[,c(which(names(demographics) == "subject"), which(names(demographics) == "sex"), which(names(demographics)=="age"), which(names(demographics) =="meanFD_avgSes"))] 
   
  
if (dataset=="PNC" | dataset =="NKI") {
  atlases <- c("schaefer200x7")
  for(i in c(1:length(atlases))){
    df_name <- paste0("edge.", atlases[i], ".", dataset)
    assign(df_name, make_EdgeDemog(atlases[i],  demographics, "subject", dataset))
  } 
  print("make_EdgeDemog complete")
} else if (dataset=="HBN" | dataset=="HCPD") { 
  atlases <- c("schaefer200x7")
  for(i in c(1:length(atlases))){
    df_name <- paste0("edge.", atlases[i], ".", dataset)
    assign(df_name, make_EdgeDemog_covbat(atlases[i],  demographics, "subject", dataset))
  } 
  print("make_EdgeDemog_covbat complete")
}


# load edge labels 
schaefer200x7_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer200x7_edge.csv")
schaefer200x7_edge <- schaefer200x7_edge[,2] 
schaefer200x7_edge <- gsub("7Networks", "Networks", schaefer200x7_edge)
  
if(dataset=="PNC"){
  
  edge.schaefer200x7.PNC <- readRDS("/cbica/projects/network_replication/output/PNC/edge/schaefer200x7_demographics_finalsample.RData")
  edge.schaefer200x7.PNC$sex <- as.factor(edge.schaefer200x7.PNC$sex)
  names(edge.schaefer200x7.PNC) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.PNC))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.PNC)
  
} else if(dataset=="NKI") {
  edge.schaefer200x7.NKI <- readRDS("/cbica/projects/network_replication/output/NKI/edge/schaefer200x7_demographics_finalsample.RData")
  edge.schaefer200x7.NKI$sex <- as.factor(edge.schaefer200x7.NKI$sex)
  names(edge.schaefer200x7.NKI) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.NKI))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.NKI)
  
  
} else if(dataset=="HCPD") {
  edge.schaefer200x7.HCPD <- readRDS("/cbica/projects/network_replication/output/HCPD/edge/schaefer200x7_demographics_finalsample_covbat.RData")
  edge.schaefer200x7.HCPD$sex <- as.factor(edge.schaefer200x7.HCPD$sex)
  names(edge.schaefer200x7.HCPD) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.HCPD))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge_covbat(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.HCPD)
  
  
} else if(dataset=="HBN") {
  edge.schaefer200x7.HBN <- readRDS("/cbica/projects/network_replication/output/HBN/edge/schaefer200x7_demographics_finalsample_covbat.RData")
  edge.schaefer200x7.HBN$sex <- as.factor(edge.schaefer200x7.HBN$sex)
  names(edge.schaefer200x7.HBN) <- gsub("7Networks", "Networks", names(edge.schaefer200x7.HBN))
  gam.edge.age.schaefer200x7 <- fitGAMs_edge_covbat(schaefer200x7_edge,"edge","schaefer200x7", dataset)
  rm(gam.edge.age.schaefer200x7)
  rm(edge.schaefer200x7.HBN)
  
} else {
  print("Provide valid dataset")
}






  