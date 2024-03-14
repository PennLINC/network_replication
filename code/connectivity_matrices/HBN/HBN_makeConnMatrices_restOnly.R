library(base64enc)
library(cifti)
library(dplyr)
library(gifti)
library(magrittr) 
library(rjson)
library(stringr)
library(tidyr)

print("This script 1) creates a list of ciftis to extract timeseries data from, and 2) concatenates timeseries and makes connectivity matrices ")


################## 
# Set Directories 
################## 
config_data <- fromJSON(file="/cbica/projects/network_replication/manuscript/code/config_HBN.json")
data_root <- config_data$data_root
cifti_lists_dir <- paste0(data_root, "connMatricesData/cifti_lists_restOnly/")
conn_matrices_dir <- paste0(data_root, "connMatricesData/connectivity_matrices_restOnly/")
sample_selection_dir <- config_data$sample_selection_data_dir
CUBIDS_csv_dir <- config_data$CUBIDS_csv_dir


if (!dir.exists(cifti_lists_dir)) {
  # If directory doesn't exist, create it
  dir.create(cifti_lists_dir, recursive = TRUE)
  print(paste("Directory", cifti_lists_dir, "created."))
} else {
  print(paste("Directory",cifti_lists_dir, "already exists."))
}

if (!dir.exists(conn_matrices_dir)) {
  # If directory doesn't exist, create it
  dir.create(conn_matrices_dir, recursive = TRUE)
  print(paste("Directory", conn_matrices_dir, "created."))
} else {
  print(paste("Directory", conn_matrices_dir, "already exists."))
}



# get subjects and filepaths
HBN_CIFTI_summary  <- readRDS("/cbica/projects/network_replication/input/HBN/sample_selection/HBN_FinalSample_filenames_restOnly.RData") 
participants <- unique(HBN_CIFTI_summary[[1]]$sub)


# paths to all the cifti files (task carit, emotion, guessing, rest for each atlas)
HBN_CIFTI_filepaths <- HBN_CIFTI_summary[[2]]$path


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
### @param list_filepath A list of paths to all cifti files of interest
### @param ciftiFiles_df A dataframe with info on subject numbers, subfolders and cifti filenames
listCifti <- function(subject, list_filepath, ciftiFiles_df){
  index_filepath <- which(str_detect(list_filepath, subject))
  atlas <- gsub("atlas-", "", str_extract(list_filepath[index_filepath], "atlas-[a-zA-Z0-9]+"))
  task <- gsub("task-", "", str_extract(list_filepath[index_filepath], "task-[a-zA-Z]+"))
  run <- str_extract(list_filepath[index_filepath], "run-[1-2]")
  run <- unlist(lapply(run, replace_na, ""))
  run <- gsub("run", "_run", run)
  
  cifti_list <- lapply(list_filepath[index_filepath], read_cifti) # creates list of cifti files in atlas order for each participant
  cifti_names <- as.data.frame(cbind(subject, atlas, task, run))
  cifti_names <- cifti_names %>% mutate(ciftiList_names = paste0(subject, "_", task, run, "_", atlas))
  names(cifti_list) <- cifti_names$ciftiList_names
  print(paste(which(participants == subject), "/", length(participants), "Cifti list done for", subject))
  saveRDS(cifti_list, paste0(cifti_lists_dir, subject, "_timeseries_restOnly.RData"))
  
  
  if(length(cifti_list) == 0){
    print(paste(subject, "missing cifti_list"))
    toReturn <- NA
  } else {
    print(paste(subject, ": concatenating", "timeseries and making connectivity matrices"))
    glasser_names <- names(cifti_list)[grep("Glasser", names(cifti_list))] # get names of timeseries with Glasser
    timeseries_glasser <- lapply(glasser_names, loadCiftiTimeseries, cifti_list) # concatenate all the timeseries (task and rest) for a given atlas
    glasser_conn <- makeConnMatrix(timeseries_glasser) # make connectivity matrix 
    print(paste(subject, ":", "glasser_conn made"))
    
    gordon_names <- names(cifti_list)[grep("Gordon", names(cifti_list))]  
    timeseries_gordon <- lapply(gordon_names, loadCiftiTimeseries, cifti_list)
    gordon_conn <- makeConnMatrix(timeseries_gordon)  
    print(paste(subject, ":", "gordon_conn made"))
    
    schaefer217_names <- names(cifti_list)[grep("Schaefer217", names(cifti_list))]  
    timeseries_schaefer217 <- lapply(schaefer217_names, loadCiftiTimeseries, cifti_list)  
    schaefer217_conn <- makeConnMatrix(timeseries_schaefer217) 
    print(paste(subject, ":", "schaefer217_conn made"))
    
    schaefer417_names <- names(cifti_list)[grep("Schaefer417", names(cifti_list))]  
    timeseries_schaefer417 <- lapply(schaefer417_names, loadCiftiTimeseries, cifti_list)  
    schaefer417_conn <- makeConnMatrix(timeseries_schaefer417)  
    print(paste(subject, ":", "schaefer417_conn made"))
    
    matrix_output <- list(glasser_conn, gordon_conn, schaefer217_conn, schaefer417_conn)
    names(matrix_output) <- c(paste0("glasser_conn"), paste0("gordon_conn"), paste0("schaefer217_conn"), paste0("schaefer417_conn"))
  } 
  print(paste(which(participants == subject), "/", length(participants), subject, "conn matrices DONE")) 
  saveRDS(matrix_output, paste0(conn_matrices_dir, subject,"_ConnMatrices.RData"))
  rm(cifti_list)
  rm(matrix_output)
  gc()
}


######################
# Make conn matrices #
######################
# make list of ciftis, extract and concatenate timeseries for task and rest scans, compute connectivity matrices for subjects, ages 5-22 
lapply(participants, listCifti, list_filepath= HBN_CIFTI_filepaths, ciftiFiles_df = HBN_CIFTI_summary)



