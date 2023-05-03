library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")
library(purrr, lib.loc="/cbica/home/luoau/Rlibs")
library(stringr, lib.loc="/cbica/home/luoau/Rlibs")
source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/gam.smooths.predict_posterior.R")

args <- commandArgs(trailingOnly = TRUE) 
dataset = args[1]
metric = args[2]

if (dataset=="HBN" | dataset=="HCPD") {
  GBC.schaefer200 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBCschaefer200_demographics_finalsample_covbat.csv", dataset)) 
  names(GBC.schaefer200) <- gsub('X17Networks', "Networks", names(GBC.schaefer200))
} else {
  
  GBC.schaefer200 <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/GBC/GBCschaefer200_demographics_finalsample.csv", dataset)) 
  names(GBC.schaefer200) <- gsub('X17Networks', "Networks", names(GBC.schaefer200))
}
  
# make sex a factor
if (dataset == "NKI"){
  GBC.schaefer200$sex <- as.factor(GBC.schaefer200$gender)
} else {
  GBC.schaefer200$sex <- as.factor(GBC.schaefer200$sex)
}


# renaming dataframes to be metric.atlas.dataset for gam functions
 
assign(paste0("GBC.schaefer200.", dataset), GBC.schaefer200) 
 

# load parcellated S-A axis
 
schaefer200x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x17_SAaxis.csv")
schaefer200x17_SAaxis$label <- gsub("17Network", "Network", schaefer200x17_SAaxis$label) 
schaefer200_SAaxis <- schaefer200x17_SAaxis 
 

# load parcel labels 
 
schaefer200x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
schaefer200.parcel.labels <- schaefer200x17.parcel.labels #gbc
 
## Function to calculate GAM Fitted Value Predictions alignment with S-A Axis
#function to estimate npd posterior draw fitted GBC values for each region 
#and compute the correlation between regional GBC fitted value and regional S-A axis rank at each age, for each draw
#executed via rerun below

np <- 200 #number of ages to get the fitted GBC values at
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
  saveRDS(fits.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorfits_correlation_byage_%3$s.RData", dataset_name, metric, atlas_name))
  rm(fits.SAaxis.posteriorcor)
  gc()
}

make_fits.SAaxis.posteriorcor_covbat <- function(metric, atlas_name, dataset_name){
  fits.SAaxis.posteriorcor <- rerun(10, compute_GBCaxis_correlation(metric, atlas_name, dataset_name)) %>% bind_cols() #np rows by npd*10 columns; each column is a draw and has the correlation between parcel fits and parcel axis ranking at each age (row)
  saveRDS(fits.SAaxis.posteriorcor, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/SAaxis_posteriorfits_correlation_byage_%3$s_covbat.RData", dataset_name, metric, atlas_name))
  rm(fits.SAaxis.posteriorcor)
  gc()
}

if (dataset=="HCPD" | dataset=="HBN"){
  if (metric=="GBC") {
    # GBC - Region-wise Posterior Smooth fits Correlation with SA Axis
    fits.SAaxis.posteriorcor_schaefer200.GBC <- make_fits.SAaxis.posteriorcor_covbat("GBC", "schaefer200", dataset)
  } else {
    print("Supply metric as arg[2]")
  }
} else {
  if (metric=="GBC") {
    # GBC - Region-wise Posterior Smooth fits Correlation with SA Axis
    fits.SAaxis.posteriorcor_schaefer200.GBC <- make_fits.SAaxis.posteriorcor("GBC", "schaefer200", dataset)
     
  } else {
    print("Supply metric as arg[2]")
  }
}




