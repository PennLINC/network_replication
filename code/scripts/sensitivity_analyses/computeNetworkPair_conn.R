
library(dplyr)
library(rjson)
library(stringr)
library(tidyr)


################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]


################## 
# Set Directories 
################## 
config_data <- fromJSON(file=sprintf("/cbica/projects/network_replication/manuscript/code/config_%1$s.json", dataset))
data_root <- config_data$data_root
conn_matrices_dir <- paste0(data_root, "connMatricesData/connectivity_matrices/")
outputs_root <- config_data$outputs_root
sample_selection_dir <- config_data$sample_selection_data_dir
metric_output_dir <- paste0(outputs_root, "sensitivity_analyses/networkpair/")

if (!dir.exists(outputs_root)) {
  # If directory doesn't exist, create it
  dir.create(outputs_root, recursive = TRUE)
  print(paste("Directory", outputs_root, "created."))
} else {
  print(paste("Directory",outputs_root, "already exists."))
}

if (!dir.exists(metric_output_dir)) {
  # If directory doesn't exist, create it
  dir.create(metric_output_dir, recursive = TRUE)
  print(paste("Directory", metric_output_dir, "created."))
} else {
  print(paste("Directory", metric_output_dir, "already exists."))
}
 

################## 
# Read files 
################## 
parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
communityAffil <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Schaefer_7Network/schaefer200x7CommunityAffiliation_final.csv")
parcel.labels <- data.frame(cbind(parcel.labels$label, communityAffil$CommAffil))
names(parcel.labels) <- c("label", "network")
networks <- sort(unique(parcel.labels$network))

commNames <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Schaefer_7Network/schaefer7CommunityNames.csv", header=F)
commNames <- cbind(commNames, c(1:7)) %>% setNames(c("commNames", "network"))
commNames$network <- as.character(commNames$network)
parcel.labels <- merge(parcel.labels, commNames)

participants <- read.table(sprintf("%1$s%2$s_final_subjectlist.txt", sample_selection_dir, dataset))
if (dataset=="NKI") {
  participants <- gsub("sub-", "", participants$V1)
} else {
  participants <- participants$V1
}



################## 
# Define Functions
################## 
# Function for computing average connectivity between and within networks
# @param subject String for id of subject of interest
# @param atlas String for name of atlas 
computeMeanNetworkConn <- function(subject, atlas) {
  print(subject)
  #read in connectivity matrix
  if(dataset == "PNC" | dataset == "HCPD" | dataset == "HBN") {
    connect.matrix <- readRDS(sprintf("%1$s%2$s_ConnMatrices.RData", conn_matrices_dir, subject))
    
    if(atlas=="schaefer200" | atlas=="schaefer400") {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(atlas_name, "_conn")]]
  } else if(dataset == "NKI") {  
    connect.matrix <- readRDS(sprintf("%1$s%2$s_ConnMatrices.RData", conn_matrices_dir, subject))
    ses_name <- str_extract(names(connect.matrix), "[A-Z]{3}1")[1]
    if(atlas=="schaefer200" | atlas=="schaefer400") {
      atlas_name <- paste0(str_extract(atlas, "schaefer[0-9]"), "17")
    } else {
      atlas_name <- atlas
    } 
    connect.matrix <- connect.matrix[[paste0(ses_name, "_", atlas_name, "_conn")]]
  } else {
    print("Provide valid dataset")
  }
  connect.matrix <- data.frame(connect.matrix) %>% setNames(parcel.labels$label)
  rownames(connect.matrix) <- parcel.labels$label
  
  mean_BNC_list <- list()
  between_list <- list() # have a different between_list for each pair of networks
  within_mat <-  matrix(data=NA, nrow=length(networks), ncol=3)
  colnames(within_mat) <- c("network", "commNames", "mean_WNC")
  
  # for parcels in each network, find the connectivity (correlation) to all parcels in each of the other networks
  # then take the average connectivity of the initial network with the other networks
  for(n in c(1:length(networks))) {
    parcels <- parcel.labels$label[which(parcel.labels$network==networks[n])]
    for(o in c(1:length(networks[-1]))) {
      outside_networks <- sort(networks[-c(which(networks == networks[n]))])
      between_mat <- matrix(data=NA, nrow=length(parcels), ncol=length(which(parcel.labels$network==outside_networks[o]))) 
      #ncol = number of parcels in outside_network[o]
      colnames(between_mat) <- c(parcel.labels$label[which(parcel.labels$network==outside_networks[o])])
      rownames(between_mat) <- c(parcel.labels$label[which(parcel.labels$network==networks[n])])
      within_vec <- c()
      
      for(p in c(1:length(parcels))) { # for each parcel in a given network:
        conn_values <- connect.matrix[[parcels[p]]] # connectivity of parcel[p] to parcels in outside_network[o]
        conn_values_df <- data.frame(bind_cols(conn_values, parcel.labels$network)) %>% 
          setNames(c("conn_values", "network"))
        
        # calculate within-network connectivity (WNC) for parcel[p]
        within_conn <- conn_values_df$conn_values[which(conn_values_df$network==networks[n])] 
        within_conn <- within_conn[-c(which(within_conn==1))] # remove self-connection
        within_vec <- append(within_vec, within_conn)
        
        # calculate between-network connectivity (BNC) of parcel[p] with outside_network[o]  
        between_conn <- conn_values_df$conn_values[which(conn_values_df$network==outside_networks[o])]  
        between_mat[p,] <- c(between_conn) 
      }
      
      # compute mean WNC for networks[n]
      mean_within_conn <- mean(within_vec, na.rm=TRUE)
      within_mat[n,] <- cbind(networks[n], 
                              unique(parcel.labels$commNames[which(parcel.labels$network==networks[n])]), 
                              mean_within_conn)
      
      # compute mean BNC for networks[n] and outside_networks[o]
      between_list[[o]] <- between_mat
      mean_BNC <- mean(as.vector(between_list[[o]]))
      network_pair <- cbind(unique(parcel.labels$commNames[which(parcel.labels$network==networks[n])]),
                            unique(parcel.labels$commNames[which(parcel.labels$network %in% outside_networks[o])]))
      mean_BNC_df <- data.frame(bind_cols(mean_BNC, network_pair)) %>% 
        setNames(c("mean_BNC", "network1", "network2"))
      mean_BNC_list <- append(mean_BNC_list, list(mean_BNC_df))
      print(paste(networks[n], outside_networks[o], "finished"))
    }
    
  }
  
  # finalize BNC and WNC dataframes 
  mean_BNC_final <- dplyr::bind_rows(mean_BNC_list)
  # no longer removing duplicates here
  #mean_BNC_final$sorted_pairs <- apply(mean_BNC_final[c("network1", "network2")], 1, 
                                       #function(x) toString(sort(x))) # remove duplicated pairs from mean_BNC_final
  #mean_BNC_final <- mean_BNC_final[!duplicated(mean_BNC_final$sorted_pairs), ] 
  #mean_BNC_final <- mean_BNC_final[, -ncol(mean_BNC_final)] # remove sorted_pairs col
   
  
  mean_WNC_temp <- data.frame(within_mat)  
  mean_WNC_temp$mean_WNC <- as.numeric(mean_WNC_temp$mean_WNC)
  mean_WNC_temp <- mean_WNC_temp %>% mutate(network1 = commNames, network2=commNames)
  mean_WNC_final <- mean_WNC_temp %>% select(mean_WNC, network1, network2)
  
  # combine BNC and WNC
  mean_BNC_final <- mean_BNC_final %>% rename(mean_connectivity = mean_BNC)
  mean_WNC_final <- mean_WNC_final %>% rename(mean_connectivity = mean_WNC)
  
  mean_network_conn <- bind_rows(mean_WNC_final, mean_BNC_final)
  mean_network_conn <- mean_network_conn %>% mutate(network_pair = paste0(network1, "_", network2)) %>%
    select(mean_connectivity, network_pair)
  
  mean_network_conn <- pivot_wider(mean_network_conn, names_from="network_pair", values_from="mean_connectivity")
  mean_network_conn <- mean_network_conn %>% mutate(subject=subject) %>% relocate(subject)
  
  print(paste(which(participants %in% subject), "/", length(participants), "-", subject))
  return(mean_network_conn)
}


################## 
# Compute Metric
################## 
subjects_mean_netconn <- lapply(participants, computeMeanNetworkConn,  atlas="schaefer200")
names(subjects_mean_netconn) <- participants

saveRDS(subjects_mean_netconn, sprintf("%1$snetworkpair_subxnetpair_matrix_schaefer200x7_orig_notbound.RData", metric_output_dir))
subjects_mean_netconn <- bind_rows(subjects_mean_netconn)

# creates a subject x network_pair csv (nrow=length(subjects), ncol=length(network_pair))
write.csv(subjects_mean_netconn, sprintf("%1$snetworkpair_subxnetpair_matrix_schaefer200x7_orig.csv", metric_output_dir), row.names=F, quote=F)


