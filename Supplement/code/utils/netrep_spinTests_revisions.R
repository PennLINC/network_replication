# code for running spin-based spatial permutation tests

 
#Load Spin Test Parcel Rotation Matrix
source("/cbica/projects/network_replication/software/perm.sphere.p.R")
 

# create df.dev.spin formatted for spatial permutation testing

# load parcel labels   
schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")

# function for creating the final df.dev.spin dataframe for permutation testing
prepDfSpin <- function(dataset, conn_type) {
  parcel.labels <- schaefer200x7.parcel.labels$label
  axis.df <- get(paste0("GBC", ".axis_", dataset, "_", conn_type))
  figureDF <- axis.df
  df.dev.spin <- rbind(figureDF[1:c(length(parcel.labels)/2),], figureDF[c(length(parcel.labels)/2+1):length(parcel.labels),]) #format df as left hemisphere -> right hemisphere for spin tests
  return(df.dev.spin)
}


  
  

 



 
