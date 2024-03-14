library(base64enc)
library(cifti)
library(dplyr)
library(gifti)
library(magrittr)
library(rjson)
library(stringr)


print("This script 1) loads a subject's cifti files, and 2) concatenates timeseries and makes thresholded / absvalue connectivity matrices")

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
conn_type = args[2] # thresholded or absvalue


################## 
# Set Directories 
################## 
config_data <- fromJSON(file="/cbica/projects/network_replication/manuscript/code/config_PNC.json")
data_root <- config_data$data_root
cifti_lists_dir <- paste0(data_root, "connMatricesData/cifti_lists/")
conn_matrices_dir <- paste0(data_root, "connMatricesData/connectivity_matrices_", conn_type, "/")
sample_selection_dir <- config_data$sample_selection_data_dir
CUBIDS_csv_dir <- config_data$CUBIDS_csv_dir


if (!dir.exists(conn_matrices_dir)) {
  # If directory doesn't exist, create it
  dir.create(conn_matrices_dir, recursive = TRUE)
  print(paste("Directory", conn_matrices_dir, "created."))
} else {
  print(paste("Directory", conn_matrices_dir, "already exists."))
}

 
################## 
# Define Functions
################## 
## This function extracts the timeseries data from cifti files and concatenates those timeseries for existing acquisitions
### @param name A character string, the name of cifti containing timeseries of interest
### @param cifti A list of ciftis containing the cifti of interest
#class(cifti_list)
loadCiftiTimeseries <- function(name, cifti) {
  timeseries <- cifti[[name]]$data
  timeseries <- t(timeseries)
  return(timeseries)
}

## This function concatenates the timeseries then creates thresholded connectivity matrix
### @param list_timeseries A list of timeseries from each available acquisition for a given subject to concatenate
makeConnMatrix_thresholded <- function(list_timeseries){
  rbound <- do.call(rbind, list_timeseries)
  corr_matrix <- round(cor(rbound, use = "complete.obs"), 5)
  corr_matrix_thresholded <- ifelse(corr_matrix > 0, corr_matrix, NA)
  return(corr_matrix_thresholded)
}

## This function concatenates the timeseries then creates connectivity matrix with absolute correlations 
### @param list_timeseries A list of timeseries from each available acquisition for a given subject to concatenate
makeConnMatrix_absvalue <-function(list_timeseries){
  rbound <- do.call(rbind, list_timeseries)
  corr_matrix_absvalue <- abs(round(cor(rbound, use = "complete.obs"), 5))
  return(corr_matrix_absvalue)
}

## This function concatenates the timeseries then creates connectivity matrix with absolute correlations 
### @param list_timeseries A list of timeseries from each available acquisition for a given subject to concatenate
makeConnMatrix <-function(list_timeseries){
  rbound <- do.call(rbind, list_timeseries)
  corr_matrix <- round(cor(rbound, use = "complete.obs"), 5)
  return(corr_matrix)
}

 
## This wrapper function does the following:
# load appropriate cifti timeseries (schaefer 200 for a given subject) using loadCiftiTimeseries
# then apply makeConnMatrix function
# and saves out! 

wrapperFunc <- function(timeseries_filename){
  timeseries_file <- readRDS(timeseries_filename) 
  if (dataset == "PNC" | dataset == "HCPD") { 
    subject <- str_extract(timeseries_filename, "sub-[0-9]+")
  } else if (dataset == "HBN") {
    subject <- str_extract(timeseries_filename, "sub-[A-Z0-9]+")
  } else if (dataset == "NKI") {
    subject <- str_extract(timeseries_filename, "A[0-9]+")
  }
  
  print(paste(subject, ": concatenating", "timeseries and making connectivity matrices"))
  
  schaefer217_names <- names(timeseries_file)[grep("Schaefer217", names(timeseries_file))]  
  timeseries_schaefer217 <- lapply(schaefer217_names, loadCiftiTimeseries, timeseries_file)  
  if (conn_type=="thresholded") {
    schaefer217_conn <- makeConnMatrix_thresholded(timeseries_schaefer217) 
  } else if (conn_type =="absvalue") {
    schaefer217_conn <- makeConnMatrix_absvalue(timeseries_schaefer217) 
  } else {
    print("Provide conn_type argument (thresholded or absvalue)")
  }
  
  matrix_output <- schaefer217_conn
  
  print(paste(which(all_sub_timeseries == timeseries_filename), "/", length(all_sub_timeseries), subject, "conn matrices DONE")) 
  saveRDS(matrix_output, paste0(conn_matrices_dir, subject, sprintf("_ConnMatrices_%1$s.RData", conn_type)))  
  rm(timeseries_file)
  rm(matrix_output)
  gc()
}
 

############################################# 
# Compute Abs Corr or Thresholded matrices #
############################################# 
all_sub_timeseries <- list.files(cifti_lists_dir, full.names=TRUE) # takes the cifti list from the main analyses and computes abs value or thresholded conn matrices
lapply(all_sub_timeseries, wrapperFunc)


