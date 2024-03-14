
library(dplyr)
library(magrittr)
library(rjson)
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
outdir_GBC <- paste0(config$outputs_root, "GBC/GAM/")
outdir_BNC <- paste0(config$outputs_root, "BNC/GAM/")
outdir_WNC <- paste0(config$outputs_root, "WNC/GAM/")


if (!dir.exists(outdir_GBC)) {
  # If directory doesn't exist, create it
  dir.create(outdir_GBC, recursive = TRUE)
  print(paste("Directory", outdir_GBC, "created."))
} else {
  print(paste("Directory", outdir_GBC, "already exists."))
}

if (!dir.exists(outdir_BNC)) {
  # If directory doesn't exist, create it
  dir.create(outdir_BNC, recursive = TRUE)
  print(paste("Directory", outdir_BNC, "created."))
} else {
  print(paste("Directory", outdir_BNC, "already exists."))
}

if (!dir.exists(outdir_WNC)) {
  # If directory doesn't exist, create it
  dir.create(outdir_WNC, recursive = TRUE)
  print(paste("Directory", outdir_WNC, "created."))
} else {
  print(paste("Directory", outdir_WNC, "already exists."))
}


################## 
# Define Functions
################## 

# Function for making connMetrics_demographics dataframe
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs (e.g. "rbcid" for PNC, "subject" for NKI)
# @param dataset A character string, name of dataset
make_connMetricsDemog <- function(atlas_network, metric, demog_csv, subj_id, dataset){
  subxparcel.matrix <- read.csv(sprintf("%1$s%2$s/%2$s_subxparcel_matrix_%3$s.csv", outputs_root, metric, atlas_network))
  print(paste(atlas_network, metric, "file loaded"))
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  connMetricsDemog_df$sex <- as.factor(connMetricsDemog_df$sex)
  write.csv(connMetricsDemog_df, sprintf("%1$s%2$s/%2$s%3$s_demographics_finalsample.csv", outputs_root, metric, atlas_network), row.names = F)
  return(connMetricsDemog_df)
}



# Function for **Fitting GAMs** for GBC, BNC, and WNC
# Fit GAM (func_conn_metric ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "BNC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
fitGAMs <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name) {
  parcel.labels <- parcel.labels$label
  gam.age <- matrix(data=NA, nrow=length(parcel.labels), ncol=9) #empty matrix to save gam.fit output to
  
  for(row in c(1:length(parcel.labels))){ #for each region
    region <- parcel.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name, network_parcellation), dataset = dataset_name, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  write.csv(gam.age, sprintf("%1$s%2$s/GAM/GAMresults.%2$s.age.%3$s%4$s.csv", outputs_root, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  return(gam.age)
}


# Function for estimating GAM smooths based on model-predicted data and save out predicted y data
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric (e.g. "GBC")
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
estimate_GAMsmooths <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name) {
  parcel.labels <- parcel.labels$label
  gam.smooths <- matrix(data=NA, ncol=7) #empty matrix to save gam.predsmooth fits to
  colnames(gam.smooths) <- c("age","fit","se.fit","selo","sehi","index","label")
  
  for(row in c(1:length(parcel.labels))){ #for each region
    
    region <- parcel.labels[row] #get the region name
    GAM.SMOOTH <- gam.predsmooth(measure = metric, atlas = paste0(atlas_name, network_parcellation), dataset = dataset_name, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.predsmooth function
    preddata <- as.data.frame(GAM.SMOOTH[3]) #get predicted.smooth df from function output
    preddata$index <- rep(x=row, 1000) #region index
    preddata$label <- rep(x=GAM.SMOOTH[1], 1000) #label
    gam.smooths <- rbind(gam.smooths, preddata)
  }
  gam.smooths <- gam.smooths[-1,] #remove empty initilization row
  gam.smooths$label <- as.character(gam.smooths$label)
  
  write.csv(gam.smooths, sprintf("%1$s%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s%4$s.csv", outputs_root, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  return(gam.smooths)
  
}




## Makes smooths df for developmental trajectories, centered on zero
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
# @param dataset A character string of dataset (i.e. "NKI")

make_gam.smooths_centered <- function(metric, atlas, dataset){
  subxparcel_dem <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample.csv", dataset, metric, atlas))
  subxparcel_dem$sex <- as.factor(subxparcel_dem$sex)
  
  #make wide_df to include parcels and age, sex, meanFD_avgSes 
  start_parcels <- which(names(subxparcel_dem) == "subject") +1
  end_parcels <- start_parcels + nrow(get(paste0(atlas, ".parcel.labels"))) - 1
  ind_age <- which(names(subxparcel_dem) == "age")
  ind_sex <- which(names(subxparcel_dem) == "sex")
  ind_meanFD_avgSes <- which(names(subxparcel_dem) == "meanFD_avgSes")
  wide_df <- subxparcel_dem[,c(start_parcels:end_parcels, ind_age, ind_sex, ind_meanFD_avgSes)]
  wide_df$age <- as.numeric(wide_df$age)
  long_df<- wide_df %>% pivot_longer(cols = !contains(c("age","sex", "meanFD_avgSes")),names_to = "parcel",values_to = metric) 
  if(metric == "GBC") { 
    smooth_fits <- long_df %>% group_by(parcel) %>% do(fit = smooth_estimates(gam(GBC ~ s(age,k=3,fx=F),data=.))) %>% unnest(fit)
  } else if(metric=="BNC") {
    smooth_fits <- long_df %>% group_by(parcel) %>% do(fit = smooth_estimates(gam(BNC ~ s(age,k=3,fx=F),data=.))) %>% unnest(fit)
  } else if (metric=="WNC") {
    smooth_fits <- long_df %>% group_by(parcel) %>% do(fit = smooth_estimates(gam(WNC ~ s(age,k=3,fx=F),data=.))) %>% unnest(fit)
  }
  smooth_fits$label <- smooth_fits$parcel
  SAaxis <- get(paste0(atlas, "_SAaxis"))
  if (str_detect(atlas, "00x7")) {
    smooth_fits$label <- gsub("X7Network", "Network", smooth_fits$label)
  } else if (str_detect(atlas, "00x17")) {
    smooth_fits$label <- gsub("X17Network", "Network", smooth_fits$label)
  } else if (atlas=="schaefer200" | atlas=="schaefer400") {
    smooth_fits$label <- gsub("X17Network", "Network", smooth_fits$label)
  }
  if (atlas=="glasser") {
    to_replace_in_glasser <- setdiff(SAaxis$label, unique(smooth_fits$label))
    glasser_labels_to_change <- gsub("-", '.', to_replace_in_glasser)
    index <- which(SAaxis$label %in% to_replace_in_glasser)
    SAaxis$label[index] <- glasser_labels_to_change
  }
  smooth_fits <- left_join(smooth_fits, SAaxis, by = "label")
  smooth_fits$SA.axis_rank <- as.integer(smooth_fits$SA.axis_rank)
  write.csv(smooth_fits, sprintf("%1$s%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_centered.csv", outputs_root, metric, atlas))
  return(smooth_fits)
} 


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
  write.csv(gam.fitted.atlas, sprintf("%1$s%2$s_%3$s.csv", outdir_GBC, dataset_name, filename), row.names = F, quote = F)
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
  saveRDS(fits.SAaxis.posteriorcor, sprintf("%1$s%2$s_SAaxis_posteriorfits_correlation_byage_%3$s.RData", outdir_GBC, dataset, atlas_name))
  rm(fits.SAaxis.posteriorcor)
  gc()
}

################## 
# Load files
################## 

# load parcel labels
glasser.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/glasser360_regionlist_final.csv")
gordon.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/gordon_regionlist_final.csv")
schaefer200x7.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
schaefer200.parcel.labels <- schaefer200x7.parcel.labels
schaefer400.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer400x7_regionlist_final.csv")


# load parcellated S-A axis
glasser_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/glasser_SAaxis.csv")
gordon_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/gordon_SAaxis.csv")
schaefer200x7_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x7_SAaxis.csv")
schaefer200x7_SAaxis$label <- gsub("7Network", "Network", schaefer200x7_SAaxis$label) 
schaefer200_SAaxis <- schaefer200x7_SAaxis
schaefer400_SAaxis  <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer400x7_SAaxis.csv")
schaefer400_SAaxis$label <-  gsub("7Networks", "Networks", schaefer400_SAaxis$label)


# load demographics
demographics <- read.csv(paste0(sample_selection_dir, dataset, "_demographics_finalsample.csv"))
demographics <- demographics %>% dplyr::rename(subject=sub)
if (dataset=="HCPD") {
  demographics$age <- demographics$interview_age/12
} else if (dataset == "NKI") {
  demographics$subject <- gsub("sub-", "", demographics$subject)
}
demographics <- demographics %>% select(subject, sex, age, meanFD_avgSes)

####################
# Fit GAMs for GBC #
####################
# create GBC.<atlas>.<dataset> dataframes and csv's
metric = "GBC"

if (dataset == "PNC") {
  atlases <- c("glasser", "gordon", "schaefer200", "schaefer400")
  for(i in c(1:length(atlases))){
    df_name <- paste0(metric, ".", atlases[i], ".", dataset)
    assign(df_name, make_connMetricsDemog(atlases[i], metric, demographics, "subject", dataset))
  }
} else {
  df_name <- paste0(metric, ".", "schaefer200", ".", dataset)
  assign(df_name, make_connMetricsDemog("schaefer200", metric, demographics, "subject", dataset))
  
}

# fit GAMs
if (dataset=="PNC") {
  gam.GBC.age.glasser <- fitGAMs(glasser.parcel.labels, "GBC", "glasser", "",dataset)
  gam.GBC.age.gordon <- fitGAMs(gordon.parcel.labels, "GBC", "gordon", "", dataset)
  gam.GBC.age.schaefer200 <- fitGAMs(schaefer200.parcel.labels, "GBC", "schaefer200", "", dataset)
  gam.GBC.age.schaefer400 <- fitGAMs(schaefer400.parcel.labels, "GBC", "schaefer400", "", dataset)
} else {
  gam.GBC.age.schaefer200 <- fitGAMs(schaefer200.parcel.labels, "GBC", "schaefer200", "", dataset)
}
print("GAMs fit for GBC")

# estimate centered GAM smooths 
smooth_fits_centered <- make_gam.smooths_centered(metric, "schaefer200", dataset)
print("GBC smooth_fits_centered made")
if (dataset=="PNC") {
  smooth_fits_centered <- make_gam.smooths_centered(metric, "glasser", dataset)
  print("GBC smooth_fits_centered made for glasser")
  smooth_fits_centered <- make_gam.smooths_centered(metric, "gordon", dataset)
  print("GBC smooth_fits_centered made gordon")
  smooth_fits_centered <- make_gam.smooths_centered(metric, "schaefer400", dataset)
  print("GBC smooth_fits_centered made for schaefer400")
  
}


# estimate GAM smooths 
smooth_fits_noncentered <- estimate_GAMsmooths(schaefer200x7.parcel.labels, "GBC", "schaefer200", "", dataset)
print("GBC smooth_fits_noncentered made")


####################
# Fit GAMs for BNC #
####################
# create BNC.schaefer200.<dataset> dataframes and csv's
metric = "BNC"
df_name <- paste0(metric, ".", "schaefer200x7", ".", dataset)
assign(df_name, make_connMetricsDemog("schaefer200x7", metric, demographics, "subject", dataset))

# fit GAM
gam.BNC.age.schaefer200x7 <- fitGAMs(schaefer200x7.parcel.labels, "BNC", "schaefer200", "x7", dataset)
print("GAMs fit for BNC")

####################
# Fit GAMs for WNC #
####################
# create WNC.<atlas>.<dataset> dataframes and csv's
metric = "WNC"
df_name <- paste0(metric, ".", "schaefer200x7", ".", dataset)
assign(df_name, make_connMetricsDemog("schaefer200x7", metric, demographics, "subject", dataset))

# fit GAM
gam.WNC.age.schaefer200x7 <- fitGAMs(schaefer200x7.parcel.labels, "WNC", "schaefer200", "x7", dataset)
print("GAMs fit for WNC")





######################################################
# Calculate Region-wise Predicted Fitted GBC values  #
######################################################

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

gam.GBC.fitted.schaefer200 <- dev_fitted(schaefer200.parcel.labels, "GBC", "schaefer200", dataset)

###########################################################################################
# Calculate GAM Fitted Value Predictions alignment with S-A Axis in Age-Resolved analysis #
###########################################################################################
fits.SAaxis.posteriorcor_schaefer200.GBC <- make_fits.SAaxis.posteriorcor("GBC", "schaefer200", dataset)


 
