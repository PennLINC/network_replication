# code for running spin-based spatial permutation tests

 
#Load Spin Test Parcel Rotation Matrix
source("/cbica/projects/network_replication/software/perm.sphere.p.R")
 

# create df.dev.spin formatted for spatial permutation testing

# load parcel labels  
glasser.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/glasser360_regionlist_final.csv")

gordon.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/gordon_regionlist_final.csv")

schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")

schaefer200x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
schaefer200.parcel.labels <- schaefer200x17.parcel.labels #gbc

schaefer400x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x7_regionlist_final.csv")

schaefer400x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x17_regionlist_final.csv")
schaefer400.parcel.labels <- schaefer400x17.parcel.labels #gbc


# function for creating the final df.dev.spin dataframe for permutation testing
prepDfSpin <- function(atlas, metric) {
  parcel.labels <- get(paste0(atlas, ".parcel.labels"))
  parcel.labels <- parcel.labels$label
  axis.df <- get(paste0(metric, ".axis_", atlas))
  figureDF <- axis.df
  if(atlas=="gordon"){
    df.dev.spin <- rbind(figureDF[1:161,], figureDF[162:333,]) #format df as left hemisphere -> right hemisphere for spin tests
  } else if(str_detect(atlas, "schaefer") | str_detect(atlas, "glasser")) {
    df.dev.spin <- rbind(figureDF[1:c(length(parcel.labels)/2),], figureDF[c(length(parcel.labels)/2+1):length(parcel.labels),]) #format df as left hemisphere -> right hemisphere for spin tests
  }  
  return(df.dev.spin)
}


  
  

 



 
