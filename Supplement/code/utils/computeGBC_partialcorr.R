library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(magrittr, lib.loc="/cbica/home/luoau/Rlibs")

atlas = "schaefer200x7"
dataset="PNC"
# only doing this for PNC, n-back only 
outdir <- sprintf("/cbica/projects/network_replication/revisions/output/%1$s/GBC/", dataset)

schaefer200x7.parcel.labels <- read.csv(sprintf("/cbica/projects/network_replication/atlases/parcellations/%1$s_regionlist_final.csv", atlas))

# Function for making empty dataframes for connectivity metric outputs (schaefer 200 only)
# @param subject_list list of participant ID's
make_output_dfs <- function(subject_list){
  n <- 201
  regionheaders <- as.character(schaefer200x7.parcel.labels$label)
  subxparcel.matrix  <- matrix(data = NA, nrow = length(subject_list), ncol = n)
  demoheaders <- c("subject")
  colheaders <- as.matrix(c(demoheaders,regionheaders))
  colnames(subxparcel.matrix) <- colheaders
  return(subxparcel.matrix)
  
}

# n-back only PNC participants 
participants <- read.csv("/cbica/projects/network_replication/revisions/input/PNC/CPAC_sample_selection/CPAC_partialcorr_PNC_demographics_finalsample_20231130.csv")
participants <- participants$sub

# Specify the path to the directory containing the files
CPAC_path <- "/cbica/projects/network_replication/revisions/input/PNC/CPAC_subfiles"

# Get a list of files ending with "partialcorr.tsv" in all subdirectories of the specified path
file_list <- list.files(path = CPAC_path, pattern = "frac2back(.*)PartialNilearn",  recursive = TRUE)
file_list_participants <- str_extract(file_list, "sub-[0-9]*")
   
file_list_participants_temp2 <- file_list_participants[which(file_list_participants %in% participants)]  

# there should be 980 file_list_participants  
final_file_list <- file_list[which(file_list_participants %in% participants)]

final_file_list <- paste0("/cbica/projects/network_replication/revisions/input/PNC/CPAC_subfiles/", final_file_list)

# Load each file into a list of data frames
data_list <- lapply(final_file_list, read.delim, header=F)
names(data_list) <- participants

# Now data_list contains your data frames, and you can access them using indices
# Function for Computing Global Brain Connectivity  
# @param subject id of subject of interest
computeGBC <- function(subject){
  connect.matrix <- data_list[[subject]]
  
  #compute average connectivity 
  GBC <- as.array(rowMeans(connect.matrix, na.rm=TRUE))
  return(GBC)
}
 
# make output df: 
##GBC.subxparcel.matrix.schaefer200x7 
GBC.subxparcel.matrix.schaefer200x7 <- make_output_dfs(participants) 

# compute GBC
for(sub in c(1:length(participants))){
  subjectID=as.character(participants[sub])
  id.data <- computeGBC(subjectID)
  GBC.subxparcel.matrix.schaefer200x7[sub,] <- cbind(subjectID, t(id.data)) #update the name of the df every iteration to GBC.subxparcel.matrix.[atlas]
  print(paste(sub, "/", length(participants), "-", subjectID))
  
}

write.csv(GBC.subxparcel.matrix.schaefer200x7, sprintf("%1$sGBC_subxparcel_matrix_partialcor.csv", outdir), row.names=F, quote=F)
