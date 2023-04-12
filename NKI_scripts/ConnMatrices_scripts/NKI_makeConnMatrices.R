
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
NKI_filepaths  <- readRDS("/cbica/projects/network_replication/input/NKI/sample_selection/NKI_FinalSample_withCUBIDS_20221219.RData") 
participants <- unique(NKI_filepaths[[1]]$sub)
demographics <- read.csv("/cbica/projects/network_replication/input/NKI/sample_selection/NKI_demographics_finalsample_20221219.csv")
participants <- participants[which(participants %in% demographics$subject)]
participants <- paste0("sub-", participants)
NKI_CIFTIfiles <- NKI_filepaths[[2]][c(which(NKI_filepaths[[2]]$subject_namesDF %in% participants)),] # number of scans x number of atlases = 993 * 4 = 3972

 
# set up functions

## This function extracts the timeseries data from cifti files and concatenates those timeseries for existing acquisitions
### @param name A character string, the name of cifti containing timeseries of interest
### @param cifti A list of ciftis containing the cifti of interest

 
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
### @param subject A character string, subject number 
### @param list_filepath A list of paths to all cifti files of interest
### @param ciftiFiles_df A dataframe with info on subject numbers, subfolders and cifti filenames

listCifti <- function(subject, list_filepath, ciftiFiles_df){
  index_filepath <- which(str_detect(list_filepath, subject))
  atlas_names <- str_extract(ciftiFiles_df$ptseries_filenamesDF[index_filepath], "atlas-(.*)_den")
  atlas_names <- gsub("atlas-", "", atlas_names)
  atlas_names <- gsub("_den", "", atlas_names)
  ses <- str_extract(ciftiFiles_df$subfolder_namesDF[index_filepath], "ses-(.*)/f")
  ses <- gsub("ses-", "", ses)
  ses <- gsub("/f", "", ses)
  acq <- str_extract(ciftiFiles_df$ptseries_filenamesDF[index_filepath], "acq-(.*)_space")
  acq <- gsub("acq-", "", acq)
  acq <- gsub("_space", "", acq)
  if(str_detect(acq[1], "VARIANT")){
    acq <- gsub("VARIANT(.*)", "", acq)
    print(paste("Variant detected for", subject))
  } else {
    acq <- acq
  }
  
  cifti_list <- lapply(list_filepath[index_filepath], read_cifti) # creates list of cifti files in atlas order for each participant
  cifti_names <- as.data.frame(cbind(subject, ses, atlas_names))
  cifti_names <- cifti_names %>% mutate(ciftiList_names = paste0(subject, "_", ses, "_", acq, "_", atlas_names))
  names(cifti_list) <- cifti_names$ciftiList_names
  file_sub <- gsub("sub-", "", subject)
  print(paste(which(participants == subject), "/", length(participants), "Cifti list done for", file_sub))
  saveRDS(cifti_list, paste0("/cbica/projects/network_replication/input/NKI/connMatricesData/cifti_lists/", file_sub,"_", ses[1], ".RData"))
  return(cifti_list)
}


 

## This function requires loadCiftiTimeseries and makeConnMatrix functions and creates list of connectivity matrices for each atlas
### @param subject A character string, subject number 
### @param cifti_list, A list of ciftis containing all acquisitions for a given subject 

computeConnMatrices <- function(subject, cifti_list){

  # if cifti_list > 0, then concatenate the timeseries for each atlas
  # if missing cifti_list, then terminate
  if(length(cifti_list[[subject]]) == 0){
    print(paste(subject, "missing cifti_list"))
    toReturn <- NA
  } else {
    print(paste(subject, ": concatenating", "timeseries and making connectivity matrices"))
    
    glasser_names <- names(cifti_list[[subject]])[grep("Glasser", names(cifti_list[[subject]]))] # get names of acquisitions with Glasser
    timeseries_glasser <- lapply(glasser_names, loadCiftiTimeseries, cifti_list[[subject]]) # concatenate all the acq's for a given atlas
    glasser_conn <- makeConnMatrix(timeseries_glasser) # make connectivity matrix 
    print(paste(subject, ":", "glasser_conn made"))
    
    
    gordon_names <- names(cifti_list[[subject]][grep("Gordon", names(cifti_list[[subject]]))])
    timeseries_gordon <- lapply(gordon_names, loadCiftiTimeseries, cifti_list[[subject]])
    gordon_conn <- makeConnMatrix(timeseries_gordon)  
    print(paste(subject, ":", "gordon_conn made"))
    
    schaefer217_names <- names(cifti_list[[subject]][grep("Schaefer217", names(cifti_list[[subject]]))])
    timeseries_schaefer217 <- lapply(schaefer217_names, loadCiftiTimeseries, cifti_list[[subject]])
    schaefer217_conn <- makeConnMatrix(timeseries_schaefer217)
    print(paste(subject, ":", "schaefer217_conn made"))
    
    schaefer417_names <- names(cifti_list[[subject]][grep("Schaefer417", names(cifti_list[[subject]]))])
    timeseries_schaefer417 <- lapply(schaefer417_names, loadCiftiTimeseries, cifti_list[[subject]])
    schaefer417_conn <- makeConnMatrix(timeseries_schaefer417)
    print(paste(subject, ":", "schaefer417_conn made"))
    
    matrix_output <- list(glasser_conn, gordon_conn, schaefer217_conn, schaefer417_conn)
    ses <- str_extract(names(cifti_list[[subject]]), "[A-Z]{3}1")[1]
    names(matrix_output) <- c(paste0(ses, "_glasser_conn"), paste0(ses, "_gordon_conn"), paste0(ses, "_schaefer217_conn"), paste0(ses, "_schaefer417_conn"))
  } 
  print(paste(which(participants == subject), "/", length(participants), subject, "matrices DONE")) 
   
  file_sub <- gsub("sub-", "", subject)
  saveRDS(matrix_output, paste0("/cbica/projects/network_replication/input/NKI/connMatricesData/connectivity_matrices/", file_sub,"_ConnMatrices.RData"))
  return(matrix_output)
}

 
# make list of ciftis (that exist), extract and concatenate timeseries, compute connectivity matrices for subjects with BAS1 or FLU1, ages 5-22
all_cifti_list <- lapply(participants, listCifti, list_filepath= NKI_CIFTIfiles$path, ciftiFiles_df = NKI_CIFTIfiles)
names(all_cifti_list) <- participants  

 
allConnMatrices <- lapply(names(all_cifti_list),computeConnMatrices, all_cifti_list) 
names(allConnMatrices) <- participants


 
