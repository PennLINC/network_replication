 
library(base64enc, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")


print("This script 1) loads each subjects' timeseries, and 2) concatenates timeseries and makes connectivity matrices ")

# get subjects  
 
PNC_CIFTI_summary <- readRDS("/cbica/projects/network_replication/revisions/input/PNC/CPAC_sample_selection/CPAC_nonGSR_fMRIinclude_CUBIDS_20231130.RData")
PNC_CIFTI_summary$task <- gsub("rest_", "rest_acq-", PNC_CIFTI_summary$task)
participants <- PNC_CIFTI_summary$sub

# Specify the path to the directory containing the files
CPAC_path <- "/cbica/projects/network_replication/revisions/input/PNC/CPAC_subfiles"

# Get a list of files ending with ".1D" in all subdirectories of the specified path
timeseries_list <- list.files(path = CPAC_path, pattern = "(.*)timeseries.1D",  recursive = TRUE) # Note: rows of .1D files correspond to timepoints, columns correspond to ROIs

# extract from timeseries_list the subject+task timeseries according to PNC_CIFTI_summary
pattern <- mapply(function(sub, task) {
  sprintf("%s/%s(.*)%s(.*)timeseries.1D", sub, sub, task)
}, PNC_CIFTI_summary$sub, PNC_CIFTI_summary$task)

# Use sapply to find filenames matching the pattern above
matching_filenames <- timeseries_list[sapply(pattern, function(p) which(grepl(p, timeseries_list)))]
file_list_participants <- str_extract(matching_filenames, "sub-[0-9]*")
  
## This function extracts the timeseries data from .1D files and concatenates those timeseries for existing acquisitions
### @param participant_id A character string, the subject ID for files containing timeseries of interest
loadTimeseries <- function(subject) {
  file_list <- matching_filenames[which(file_list_participants %in% subject)]  
  file_list <- paste0(CPAC_path, "/", file_list)
  timeseries <- lapply(file_list, read.csv) # note that CPAC is different from xcp in that rows are timepoints and cols are ROIs
  return(timeseries)
}


## This function concatenates the timeseries then creates connectivity matrix
### @param list_timeseries A list of timeseries from each available acquisition for a given subject to concatenate
makeConnMatrix <- function(list_timeseries){
  rbound <- do.call(rbind, list_timeseries)
  corr_matrix <- round(cor(rbound, use = "complete.obs"), 5)
  return(corr_matrix)
}


## This function reads the cifti's for each subject and creates a list of ciftis containing all acquisitions for a given subject
## then creates list of connectivity matrices for each atlas
### @param subject A character string, subject number  
timeseriesWrapper <- function(subject){
    
  timeseries_schaefer217 <- loadTimeseries(subject) 
  schaefer217_conn <- makeConnMatrix(timeseries_schaefer217) 
  print(paste(subject, ":", "schaefer217_conn made"))

  print(paste(which(unique(participants) == subject), "/", length(unique(participants)), subject, "conn matrices DONE")) 
  saveRDS(schaefer217_conn, paste0("/cbica/projects/network_replication/revisions/input/PNC/CPAC_connMatricesData/connectivity_matrices/", subject,"_ConnMatrices.RData"))
}

# make list of ciftis, extract and concatenate timeseries for task and rest scans, compute connectivity matrices for subjects, ages 8-23
lapply(unique(participants), timeseriesWrapper)



 
