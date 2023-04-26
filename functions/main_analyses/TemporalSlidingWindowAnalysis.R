library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(purrr, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/GAM_functions.R")

args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
metric = args[2]

GBC.glasser <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBCglasser_demographics_finalsample.csv", dataset)) 
GBC.gordon <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBCgordon_demographics_finalsample.csv", dataset)) 

GBC.schaefer200 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBCschaefer200_demographics_finalsample.csv", dataset)) 
names(GBC.schaefer200) <- gsub('X17Networks', "Networks", names(GBC.schaefer200))

GBC.schaefer400 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBCschaefer400_demographics_finalsample.csv", dataset)) 
names(GBC.schaefer400) <- gsub('X17Networks', "Networks", names(GBC.schaefer400)) 


# BNC
BNC.gordon <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/BNC/BNCgordon_demographics_finalsample.csv", dataset))   
BNC.schaefer200x7 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/BNC/BNCschaefer200x7_demographics_finalsample.csv", dataset)) 
names(BNC.schaefer200x7) <- gsub("X7", "", names(BNC.schaefer200x7))

BNC.schaefer200x17 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/BNC/BNCschaefer200x17_demographics_finalsample.csv", dataset))    
names(BNC.schaefer200x17) <- gsub("X17", "", names(BNC.schaefer200x17))

BNC.schaefer400x7 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/BNC/BNCschaefer400x7_demographics_finalsample.csv", dataset))
names(BNC.schaefer400x7) <- gsub("X7", "", names(BNC.schaefer400x7))

BNC.schaefer400x17 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/BNC/BNCschaefer400x17_demographics_finalsample.csv", dataset))   
names(BNC.schaefer400x17) <- gsub("X17", "", names(BNC.schaefer400x17))

# WNC
WNC.gordon <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNCgordon_demographics_finalsample.csv", dataset))    
WNC.schaefer200x7 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNCschaefer200x7_demographics_finalsample.csv", dataset))
names(WNC.schaefer200x7) <- gsub("X7", "", names(WNC.schaefer200x7))

WNC.schaefer200x17 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNCschaefer200x17_demographics_finalsample.csv", dataset)) 
names(WNC.schaefer200x17) <- gsub("X17", "", names(WNC.schaefer200x17))
WNC.schaefer400x7 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNCschaefer400x7_demographics_finalsample.csv", dataset)) 

names(WNC.schaefer400x7) <- gsub("X7", "", names(WNC.schaefer400x7))
WNC.schaefer400x17 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/WNC/WNCschaefer400x17_demographics_finalsample.csv", dataset))  
names(WNC.schaefer400x17) <- gsub("X17", "", names(WNC.schaefer400x17))

# make sex a factor
if (dataset == "NKI"){
  GBC.glasser$sex <- as.factor(GBC.glasser$gender)
  GBC.gordon$sex <- as.factor(GBC.gordon$gender)
  GBC.schaefer200$sex <- as.factor(GBC.schaefer200$gender)
  GBC.schaefer400$sex <- as.factor(GBC.schaefer400$gender)
  
  
  BNC.gordon$sex <- as.factor(BNC.gordon$gender)
  BNC.schaefer200x7$sex <- as.factor(BNC.schaefer200x7$gender)
  BNC.schaefer200x17$sex <- as.factor(BNC.schaefer200x17$gender)
  BNC.schaefer400x7$sex <- as.factor(BNC.schaefer400x7$gender)
  BNC.schaefer400x17$sex <- as.factor(BNC.schaefer400x17$gender)
  
  
  WNC.gordon$sex <- as.factor(WNC.gordon$gender)   
  WNC.schaefer200x17$sex <- as.factor(WNC.schaefer200x17$gender) 
  WNC.schaefer400x17$sex <- as.factor(WNC.schaefer400x17$gender) 
  WNC.schaefer200x7$sex <- as.factor(WNC.schaefer200x7$gender) 
  WNC.schaefer400x7$sex <- as.factor(WNC.schaefer400x7$gender) 
} else {
  GBC.glasser$sex <- as.factor(GBC.glasser$sex)
  GBC.gordon$sex <- as.factor(GBC.gordon$sex)
  GBC.schaefer200$sex <- as.factor(GBC.schaefer200$sex)
  GBC.schaefer400$sex <- as.factor(GBC.schaefer400$sex)
  
  
  BNC.gordon$sex <- as.factor(BNC.gordon$sex)
  BNC.schaefer200x7$sex <- as.factor(BNC.schaefer200x7$sex)
  BNC.schaefer200x17$sex <- as.factor(BNC.schaefer200x17$sex)
  BNC.schaefer400x7$sex <- as.factor(BNC.schaefer400x7$sex)
  BNC.schaefer400x17$sex <- as.factor(BNC.schaefer400x17$sex)
  
  
  WNC.gordon$sex <- as.factor(WNC.gordon$sex)   
  WNC.schaefer200x17$sex <- as.factor(WNC.schaefer200x17$sex) 
  WNC.schaefer400x17$sex <- as.factor(WNC.schaefer400x17$sex) 
  WNC.schaefer200x7$sex <- as.factor(WNC.schaefer200x7$sex) 
  WNC.schaefer400x7$sex <- as.factor(WNC.schaefer400x7$sex) 
}


# renaming dataframes to be metric.atlas.dataset for gam functions
assign(paste0("GBC.glasser.", dataset), GBC.glasser) 
assign(paste0("GBC.gordon.", dataset), GBC.gordon) 
assign(paste0("GBC.schaefer200.", dataset), GBC.schaefer200) 
assign(paste0("GBC.schaefer400.", dataset), GBC.schaefer400) 

assign(paste0("BNC.gordon.", dataset), BNC.gordon) 
assign(paste0("BNC.schaefer200x7.", dataset), BNC.schaefer200x7) 
assign(paste0("BNC.schaefer200x17.", dataset), BNC.schaefer200x17) 
assign(paste0("BNC.schaefer400x7.", dataset), BNC.schaefer400x7) 
assign(paste0("BNC.schaefer400x17.", dataset), BNC.schaefer400x17) 


assign(paste0("WNC.gordon.", dataset), WNC.gordon) 
assign(paste0("WNC.schaefer200x7.", dataset), WNC.schaefer200x7) 
assign(paste0("WNC.schaefer200x17.", dataset), WNC.schaefer200x17) 
assign(paste0("WNC.schaefer400x7.", dataset), WNC.schaefer400x7) 
assign(paste0("WNC.schaefer400x17.", dataset), WNC.schaefer400x17) 

# load parcellated S-A axis
glasser_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/glasser_SAaxis.csv")
gordon_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/gordon_SAaxis.csv")

schaefer200x7_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x7_SAaxis.csv")
schaefer200x7_SAaxis$label <- gsub("7Network", "Network", schaefer200x7_SAaxis$label) 

schaefer200x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x17_SAaxis.csv")
schaefer200x17_SAaxis$label <- gsub("17Network", "Network", schaefer200x17_SAaxis$label) 
schaefer200_SAaxis <- schaefer200x17_SAaxis 

schaefer400x7_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer400x7_SAaxis.csv")
schaefer400x7_SAaxis$label <- gsub("7Network", "Network", schaefer400x7_SAaxis$label)

schaefer400x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer400x17_SAaxis.csv")
schaefer400x17_SAaxis$label <- gsub("17Network", "Network", schaefer400x17_SAaxis$label)
schaefer400_SAaxis <- schaefer400x17_SAaxis

# load parcel labels 

glasser.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/glasser360_regionlist_final.csv")

gordon.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/gordon_regionlist_final.csv")

schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")

schaefer200x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
schaefer200.parcel.labels <- schaefer200x17.parcel.labels #gbc

schaefer400x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x7_regionlist_final.csv")

schaefer400x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x17_regionlist_final.csv")
schaefer400.parcel.labels <- schaefer400x17.parcel.labels #gbc


## Function to calculate Region-wise Posterior Smooth Derivatives Correlation with SA Axis 
#function to estimate npd posterior draw derivatives for each region 
#and compute the correlation between regional derivative and regional S-A axis rank at each age, for each draw
#executed via rerun below

np <- 200 #number of ages to get the derivative at
npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)


compute_axis_correlation <- function(metric, atlas_name, dataset_name){ 
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  parcel.labels <- get(paste0(atlas_name, ".parcel.labels"))
  parcel.labels <- parcel.labels$label
  print("Get posterior derivatives")
  gam.derivatives.atlas <- map_dfr(parcel.labels, function(x){
    gam.derivatives(measure = metric, atlas = atlas_name, dataset = dataset_name, region = as.character(x), smooth_var = "age", covariates = "sex + meanFD_avgSes", knots = 3, set_fx = TRUE, 
                    draws = npd, increments = np, return_posterior_derivatives = TRUE)}) #run gam.derivatives to get simulated derivatives
  # Compute and save correlations with S-A axis
  gam.derivatives.atlas <- left_join(SAaxis, gam.derivatives.atlas, by="label", sort=F) #assign axis rank to each label
  corr_values <- gam.derivatives.atlas %>%
    group_by(draw,age) %>%  
    do(SAcorrelation = cor(as.numeric(.$SA.axis_rank), as.numeric(.$posterior.derivative), method=c("spearman"))) %>% #correlation between parcel rank and derivative at each age for each draw 
    # SAcorrelation is computed  by taking the correlation between parcel rank and parcel posterior derivative at a given draw and age. 
    # So for draw 1 of age 8, take the SAcorrelation for all the parcels (a correlation of length(atlas) parcels). 
    # Then for draw 1 of age 8.067, again take the correlation of all the parcels. Etc.
    unnest(cols = c(SAcorrelation))
  corr_values.wide <- corr_values %>% pivot_wider(names_from = "draw", values_from = "SAcorrelation", names_sort = FALSE)
  corr_values.wide <- corr_values.wide %>% select(contains("draw"))
  print(paste(metric, atlas_name, "finished"))
  rm(gam.derivatives.atlas)
  gc()
  return(corr_values.wide)
}

make_deriv.SAaxis.posteriorcor <- function(metric, atlas_name, dataset_name){
  deriv.SAaxis.posteriorcor <- rerun(10, compute_axis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; 
  # each column is a draw and has the correlation between parcel derivative and parcel axis ranking at each age (row)
  # basically is a df of SA correlations across all parcels, for a given draw and age
  saveRDS(deriv.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorderivative_correlation_byage_%3$s.RData", dataset_name, metric, atlas_name))
  rm(deriv.SAaxis.posteriorcor)
  gc()
}

make_deriv.SAaxis.posteriorcor_covbat <- function(metric, atlas_name, dataset_name){
  deriv.SAaxis.posteriorcor <- rerun(10, compute_axis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; each column is a draw and has the correlation between parcel derivative and parcel axis ranking at each age (row)
  write.table(deriv.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorderivative_correlation_byage_%3$s_covbat.csv", dataset_name, metric, atlas_name), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")
  rm(deriv.SAaxis.posteriorcor)
  gc()
}

if (dataset="HCPD" | dataset="HBN"){
  if (metric=="GBC") {
    # GBC - Region-wise Posterior Smooth Derivatives Correlation with SA Axis
    deriv.SAaxis.posteriorcor_glasser.GBC <- make_deriv.SAaxis.posteriorcor_covbat("GBC", "glasser", dataset)
    deriv.SAaxis.posteriorcor_gordon.GBC <- make_deriv.SAaxis.posteriorcor_covbat("GBC", "gordon", dataset)
    deriv.SAaxis.posteriorcor_schaefer200.GBC <- make_deriv.SAaxis.posteriorcor_covbat("GBC", "schaefer200", dataset)
    deriv.SAaxis.posteriorcor_schaefer400.GBC <- make_deriv.SAaxis.posteriorcor_covbat("GBC", "schaefer400", dataset)
  } else if (metric=="BNC") {
    # BNC - Region-wise Posterior Smooth Derivatives Correlation with SA Axis
    deriv.SAaxis.posteriorcor_gordon.BNC <- make_deriv.SAaxis.posteriorcor_covbat("BNC", "gordon", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x7.BNC <- make_deriv.SAaxis.posteriorcor_covbat("BNC", "schaefer200x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x7.BNC <- make_deriv.SAaxis.posteriorcor_covbat("BNC", "schaefer400x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x17.BNC <- make_deriv.SAaxis.posteriorcor_covbat("BNC", "schaefer200x17", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x17.BNC <- make_deriv.SAaxis.posteriorcor_covbat("BNC", "schaefer400x17", dataset)
    
  } else if (metric=="WNC") {
    # WNC - Region-wise Posterior Smooth Derivatives Correlation with SA Axis
    deriv.SAaxis.posteriorcor_gordon.WNC <- make_deriv.SAaxis.posteriorcor_covbat("WNC", "gordon", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x7.WNC <- make_deriv.SAaxis.posteriorcor_covbat("WNC", "schaefer200x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x7.WNC <- make_deriv.SAaxis.posteriorcor_covbat("WNC", "schaefer400x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x17.WNC <- make_deriv.SAaxis.posteriorcor_covbat("WNC", "schaefer200x17", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x17.WNC <- make_deriv.SAaxis.posteriorcor_covbat("WNC", "schaefer400x17", dataset)
    
  } else {
    print("Supply metric as arg[2]")
  }
} else {
  if (metric=="GBC") {
    # GBC - Region-wise Posterior Smooth Derivatives Correlation with SA Axis
    deriv.SAaxis.posteriorcor_glasser.GBC <- make_deriv.SAaxis.posteriorcor("GBC", "glasser", dataset)
    deriv.SAaxis.posteriorcor_gordon.GBC <- make_deriv.SAaxis.posteriorcor("GBC", "gordon", dataset)
    deriv.SAaxis.posteriorcor_schaefer200.GBC <- make_deriv.SAaxis.posteriorcor("GBC", "schaefer200", dataset)
    deriv.SAaxis.posteriorcor_schaefer400.GBC <- make_deriv.SAaxis.posteriorcor("GBC", "schaefer400", dataset)
  } else if (metric=="BNC") {
    # BNC - Region-wise Posterior Smooth Derivatives Correlation with SA Axis
    deriv.SAaxis.posteriorcor_gordon.BNC <- make_deriv.SAaxis.posteriorcor("BNC", "gordon", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x7.BNC <- make_deriv.SAaxis.posteriorcor("BNC", "schaefer200x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x7.BNC <- make_deriv.SAaxis.posteriorcor("BNC", "schaefer400x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x17.BNC <- make_deriv.SAaxis.posteriorcor("BNC", "schaefer200x17", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x17.BNC <- make_deriv.SAaxis.posteriorcor("BNC", "schaefer400x17", dataset)
    
  } else if (metric=="WNC") {
    # WNC - Region-wise Posterior Smooth Derivatives Correlation with SA Axis
    deriv.SAaxis.posteriorcor_gordon.WNC <- make_deriv.SAaxis.posteriorcor("WNC", "gordon", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x7.WNC <- make_deriv.SAaxis.posteriorcor("WNC", "schaefer200x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x7.WNC <- make_deriv.SAaxis.posteriorcor("WNC", "schaefer400x7", dataset)
    deriv.SAaxis.posteriorcor_schaefer200x17.WNC <- make_deriv.SAaxis.posteriorcor("WNC", "schaefer200x17", dataset)
    deriv.SAaxis.posteriorcor_schaefer400x17.WNC <- make_deriv.SAaxis.posteriorcor("WNC", "schaefer400x17", dataset)
    
  } else {
    print("Supply metric as arg[2]")
  }
}

 


