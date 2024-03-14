library(dplyr)
library(purrr)
library(rjson)
library(stringr)
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
outputs_root <- config$outputs_root
gam_dir <- paste0(outputs_root, "GBC/GAM/")

if (!dir.exists(gam_dir)) {
  # If directory doesn't exist, create it
  dir.create(outdir_GBC, recursive = TRUE)
  print(paste("Directory", gam_dir, "created."))
} else {
  print(paste("Directory", gam_dir, "already exists."))
}

################## 
# Define Functions
################## 

## Function to calculate Region-wise Fitted GBC values  
# derived from https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/developmental_effects/hierarchical_development.Rmd#L699-L709)
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param dataset_name A character string of dataset (e.g. "PNC")

dev_fitted <- function(parcel.labels, metric, atlas_name, dataset_name){
  np <- 200 #number of ages to get the fitted at
  npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)
  
  gam.fitted.atlas <- matrix(data=NA, ncol=9)
  colnames(gam.fitted.atlas) <- c("age","fitted","se","lower","upper","significant","significant.fit","index","label")
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  
  if(metric == "GBC" & str_detect(atlas_name, "x")){
    atlas_name <- gsub("x(.*)", "", atlas_name)
  } else {
    atlas_name <- atlas_name
  }
  for(row in c(1:nrow(SAaxis))){ #for each region
    region <- SAaxis$label[row] 
    GAM.fitted.output <- gam.smooth.predict_posterior(measure = metric, atlas = atlas_name, dataset = dataset_name, 
                                                      region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes", 
                                                      knots = 3, set_fx = TRUE, draws = npd, increments = np, return_posterior_fits = FALSE) #run the gam.fitted.atlas function to get true model fitted
    GAM.fitted.output$index <- as.numeric(rep(x=row, np)) #region index
    GAM.fitted.output$label <- as.character(rep(x=region, np)) #region label
    gam.fitted.atlas <- rbind(gam.fitted.atlas, GAM.fitted.output)
  }
  
  gam.fitted.atlas <- gam.fitted.atlas[-1,] #remove empty initialization row
  gam.fitted.atlas <- left_join(gam.fitted.atlas, SAaxis, by="label", sort = F)
  filename <- paste0("gam.", metric, ".fitted.", atlas_name)
  write.csv(gam.fitted.atlas, sprintf("%1$s%2$s_%3$s.csv", gam_dir, dataset_name, filename), row.names = F, quote = F)
  return(gam.fitted.atlas)
}


## Function to calculate GAM Fitted Value Predictions alignment with S-A Axis
#function to estimate npd posterior draw derivs for each region 
#and compute the correlation between regional GBC fitted value and regional S-A axis rank at each age, for each draw
#executed via rerun below

np <- 200 #number of ages to get the deriv at
npd <- 1000 #number of posterior draws for each run of rerun (repeated 10 times)

compute_GBCaxis_correlation <- function(metric, atlas_name, dataset_name){ 
  SAaxis <- get(paste0(atlas_name, "_SAaxis"))
  parcel.labels <- get(paste0(atlas_name, ".parcel.labels"))
  parcel.labels <- parcel.labels$label
  print("Get posterior predicted fitted values")
  gam.fits.atlas <- map_dfr(parcel.labels, function(x){
    gam.smooth.predict_posterior(measure = metric, atlas = atlas_name, dataset = dataset_name, region = as.character(x), smooth_var = "age", covariates = "sex + meanFD_avgSes", knots = 3, set_fx = TRUE, 
                                 draws = npd, increments = np, return_posterior_fits = TRUE)}) #run gam.fits to get simulated fits
  # Compute and save correlations with S-A axis
  gam.fits.atlas <- left_join(SAaxis, gam.fits.atlas, by="label", sort=F) #assign axis rank to each label
  corr_values <- gam.fits.atlas %>%
    group_by(draw,age) %>%  
    do(SAcorrelation = cor(as.numeric(.$SA.axis_rank), as.numeric(.$posterior.fits), method=c("spearman"))) %>% #correlation between parcel rank and fit at each age for each draw 
    # SAcorrelation is computed  by taking the correlation between parcel rank and parcel posterior fit at a given draw and age. 
    # So for draw 1 of age 8, take the SAcorrelation for all the parcels (a correlation of length(atlas) parcels). 
    # Then for draw 1 of age 8.067, again take the correlation of all the parcels. Etc.
    unnest(cols = c(SAcorrelation))
  corr_values.wide <- corr_values %>% pivot_wider(names_from = "draw", values_from = "SAcorrelation", names_sort = FALSE)
  corr_values.wide <- corr_values.wide %>% select(contains("draw"))
  print(paste(metric, atlas_name, "finished"))
  rm(gam.fits.atlas)
  gc()
  return(corr_values.wide)
}


make_fits.SAaxis.posteriorcor <- function(metric, atlas_name, dataset_name){
  fits.SAaxis.posteriorcor <- rerun(10, compute_GBCaxis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; 
  # each column is a draw and has the correlation between parcel fits and parcel axis ranking at each age (row)
  # basically is a df of SA correlations across all parcels, for a given draw and age
  saveRDS(fits.SAaxis.posteriorcor, sprintf("%1$s%2$s_SAaxis_posteriorfits_correlation_byage_%3$s.RData", gam_dir, dataset, atlas_name))
  rm(fits.SAaxis.posteriorcor)
  gc()
}


################## 
# Load files
################## 

# load subxparcel + demographics df
GBC.schaefer200 <- read.csv(sprintf("%1$sGBC/GBCschaefer200_demographics_finalsample.csv", outputs_root)) 
GBC.schaefer200$sex <- as.factor(GBC.schaefer200$sex)

# rename dataframe to be metric.atlas.dataset for gam functions
assign(paste0("GBC.schaefer200.", dataset), GBC.schaefer200) 
 
# load parcellated S-A axis
schaefer200_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x7_SAaxis.csv")
schaefer200_SAaxis$label <- gsub("7Network", "Network", schaefer200_SAaxis$label) 
 
# load parcel labels 
schaefer200.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
 
 

######################################################
# Calculate Region-wise Predicted Fitted GBC values  #
######################################################
gam.GBC.fitted.schaefer200 <- dev_fitted(schaefer200.parcel.labels, "GBC", "schaefer200", dataset)

###########################################################################################
# Calculate GAM Fitted Value Predictions alignment with S-A Axis in Age-Resolved analysis #
###########################################################################################
fits.SAaxis.posteriorcor_schaefer200.GBC <- make_fits.SAaxis.posteriorcor("GBC", "schaefer200", dataset)
 


