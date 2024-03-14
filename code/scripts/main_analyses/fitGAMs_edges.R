library(dplyr)
library(gratia)
library(janitor)
library(magrittr)
library(mgcv)
library(purrr)
library(rjson)
library(tidyr)
library(tidyverse)
source("/cbica/projects/network_replication/manuscript/code/scripts/main_analyses/GAM_functions.R")

################## 
# Set Variables 
################## 
args <- commandArgs(trailingOnly=TRUE)
config_file <- args[1]


################## 
# Set Directories 
################## 
config <- fromJSON(file = config_file)
dataset <- config$dataset
sample_selection_dir <- config$sample_selection_data_dir

outputs_root <- config$outputs_root
outdir_edge <- paste0(config$outputs_root, "edge/GAM/")

if (!dir.exists(outdir_edge)) {
  # If directory doesn't exist, create it
  dir.create(outdir_edge, recursive = TRUE)
  print(paste("Directory", outdir_edge, "created."))
} else {
  print(paste("Directory", outdir_edge, "already exists."))
}


################## 
# Define Functions
################## 

# Function for making edges_demographics dataframe
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_EdgeDemog <- function(atlas_network, demog_csv, subj_id, dataset){
  subxedge.matrix <- readRDS(sprintf("%1$sedge/subxedge_%2$s.RData", outputs_root, atlas_network))
  print("subxedge file loaded")
  subxedge.matrix[subxedge.matrix == 1] <- NA
  subxedge.matrix <- remove_empty(subxedge.matrix, which = "cols") # remove columns that are all ==1 (connection to self)
  print("rm columns of self-connections")
  if(dataset=="HCPD" | dataset=="HBN") {
    subxedge.matrix$subject <- rownames(subxedge.matrix) 
    subxedge.matrix <- subxedge.matrix %>% relocate(subject)
  }
  edgeDemog_df <- merge(subxedge.matrix, demog_csv, by=subj_id)
  edgeDemog_df$sex <- as.factor(edgeDemog_df$sex)
  names(edgeDemog_df) <- gsub("7Networks", "Networks", names(edgeDemog_df))
  saveRDS(edgeDemog_df, sprintf("%1$sedge/%2$s_demographics_finalsample.RData", outputs_root, atlas_network))
  print(paste(atlas_network, dataset, "done"))
  return(edgeDemog_df)
}


# Function for **Fitting GAMs** for parcel-by-parcel edges
# Fit GAM (edge ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param edge.labels A vector of edge.labels
# @param metric A character string of connectivity metric (e.g. "edge")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs_edge <- function(edge.labels, metric, atlas_name, dataset_name) {
  gam.age <- matrix(data=NA, nrow=length(edge.labels), ncol=9) #empty matrix to save gam.fit output to
  for(row in c(1:length(edge.labels))){ #for each region
    edge_name <- edge.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name), dataset = dataset_name, region = edge_name, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS
    print(paste(row, "/", length(edge.labels)))}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  saveRDS(gam.age, sprintf("%1$sGAMresults.%2$s.age.%3$s.RData", outdir_edge, metric, atlas_name))
  return(gam.age)
}


################## 
# Load files
################## 

# load edge labels 
schaefer200x7_edge <- read.csv("/cbica/projects/network_replication/atlases/edge/schaefer200x7_edge.csv")
schaefer200x7_edge <- schaefer200x7_edge[,2] 
schaefer200x7_edge <- gsub("7Networks", "Networks", schaefer200x7_edge)

# load demographics
demographics <- read.csv(paste0(sample_selection_dir, dataset, "_demographics_finalsample.csv"))
demographics <- demographics %>% dplyr::rename(subject=sub)
if (dataset=="HCPD") {
  demographics$age <- demographics$interview_age/12
} else if (dataset == "NKI") {
  demographics$subject <- gsub("sub-", "", demographics$subject)
}
demographics <- demographics %>% select(subject, sex, age, meanFD_avgSes)


######################
# Fit GAMs for Edges #
######################
# create edge.<atlas>.<dataset> dataframes and csv's
metric = "edge"
df_name <- paste0(metric, ".", "schaefer200x7", ".", dataset)
assign(df_name, make_EdgeDemog("schaefer200x7", demographics, "subject", dataset))
print("make_EdgeDemog complete")
 
  
# fit GAMs on edges
print("fitting GAMs on edges")
gam.edge.age.schaefer200x7 <- fitGAMs_edge(schaefer200x7_edge,"edge","schaefer200x7", dataset)
rm(gam.edge.age.schaefer200x7)
rm(get(df_name))
print("done fitting GAMs")
 

 