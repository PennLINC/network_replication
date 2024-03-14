# note that lib.loc needs to be specified since this .R file gets submitted as a job on cubic
library(base64enc, lib.loc="/cbica/home/luoau/Rlibs")
library(gifti, lib.loc="/cbica/home/luoau/Rlibs")
library(cifti, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")

args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
conn_type = args[2] # thresholded or absvalue

outdir <- sprintf("/cbica/projects/network_replication/revisions/input/%1$s/connMatricesData/", dataset)
print("This script 1) loads a subject's cifti files, and 2) concatenates timeseries and makes thresholded connectivity matrices")
 
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
  saveRDS(matrix_output, paste0(outdir, "connectivity_matrices/", subject, sprintf("_ConnMatrices_%1$s.RData", conn_type)))  
  rm(timeseries_file)
  rm(matrix_output)
  gc()
}
 
all_sub_timeseries <- list.files(sprintf("/cbica/projects/network_replication/input/%1$s/connMatricesData/cifti_lists", dataset), full.names=TRUE)
# note that for NKI, the session that has the most scans surviving the head motion exclusion was used. This selection was done during sample selection.
# the cifti_lists directory for NKI contains the final session used in analyses 

lapply(all_sub_timeseries, wrapperFunc)


