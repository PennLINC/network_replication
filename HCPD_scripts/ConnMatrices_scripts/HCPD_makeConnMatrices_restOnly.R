
#library(ciftiTools, lib.loc="/cbica/home/luoau/Rlibs")
#ciftiTools.setOption('wb_path', '/cbica/share/modules/connectome_workbench/1.4.2')
library(base64enc, lib.loc="/cbica/home/luoau/Rlibs")
library(gifti, lib.loc="/cbica/home/luoau/Rlibs")
library(cifti, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")

print("This script 1) creates a list of ciftis to extract timeseries data from, and 2) concatenates timeseries and makes connectivity matrices ")

# get subjects and filepaths
HCPD_CIFTI_summary  <- readRDS("/cbica/projects/network_replication/input/HCPD/sample_selection/HCPD_FinalSample_withCUBIDS_restOnly_20230203.RData") 
participants <- unique(HCPD_CIFTI_summary[[1]]$sub)
#demographics <- read.csv("/cbica/projects/network_replication/input/HCPD/hcpd_demographics.csv")


# paths to all the cifti files (task carit, emotion, guessing, rest for each atlas)
HCPD_CIFTI_filepaths <- HCPD_CIFTI_summary[[2]]$path



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
  task <- gsub("task-", "", str_extract(list_filepath[index_filepath], "task-[a-z]+"))
  run <- gsub("_", "", str_extract(list_filepath[index_filepath], "_run-[0-9]{2}")) # need to add "run" since HCP-D has multiple runs per rest scan for some participants
  direction <- gsub("dir-", "", str_extract(list_filepath[index_filepath], "dir-[A-Z]{2}"))
  
  cifti_list <- lapply(list_filepath[index_filepath], read_cifti) # creates list of cifti files in atlas order for each participant
  cifti_names <- as.data.frame(cbind(subject, atlas, task, direction, run))
  cifti_names <- cifti_names %>% mutate(ciftiList_names = paste0(subject, "_", task, "_", run,  "_", direction, "_", atlas))
  names(cifti_list) <- cifti_names$ciftiList_names
  print(paste(which(participants == subject), "/", length(participants), "Cifti list done for", subject))
  saveRDS(cifti_list, paste0("/cbica/projects/network_replication/input/HCPD/connMatricesData/cifti_lists/", subject, "_timeseriesTaskandRest.RData"))
  
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
  saveRDS(matrix_output, paste0("/cbica/projects/network_replication/input/HCPD/connMatricesData/connectivity_matrices_restOnly/", subject,"_ConnMatrices.RData"))
  rm(cifti_list)
  rm(matrix_output)
  gc()
}

# make list of ciftis, extract and concatenate timeseries for task and rest scans, compute connectivity matrices for subjects, ages 5-22 
lapply(participants, listCifti, list_filepath= HCPD_CIFTI_filepaths, ciftiFiles_df = HCPD_CIFTI_summary)



