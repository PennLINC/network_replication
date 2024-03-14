
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
metric_output_dir <- paste0(outputs_root, "sensitivity_analyses/GBC/")
gam_dir <- paste0(metric_output_dir, "GAM/")

if (!dir.exists(gam_dir)) {
  # If directory doesn't exist, create it
  dir.create(gam_dir, recursive = TRUE)
  print(paste("Directory", gam_dir, "created."))
} else {
  print(paste("Directory", gam_dir, "already exists."))
}

################## 
# Define Functions
################## 

# Function for making connMetrics_demographics dataframe
# @param atlas_network A character string, name of atlas of interest (e.g. "gordon", "schaefer200x7", "schaefer400x17" - specify network parcellation for schaefer atlases if metric = BNC or WNC)
# @param metric A character string, name of connectivity metric ("GBC", "BNC" or "WNC")
# @param demog_csv A df of demographic info
# @param subj_id A character string of column name of subject IDs 
# @param dataset A character string, name of dataset
make_connMetricsDemog <- function(atlas_network, metric, demog_csv, subj_id, dataset){
  subxparcel.matrix <- read.csv(sprintf("%1$ssensitivity_analyses/%2$s/%2$s_subxparcel_matrix_%3$s.csv", outputs_root, metric, atlas_network))
  print(paste(atlas_network, metric, "file loaded"))
  connMetricsDemog_df <- merge(subxparcel.matrix, demog_csv, by=subj_id)
  connMetricsDemog_df$sex <- as.factor(connMetricsDemog_df$sex)
  write.csv(connMetricsDemog_df, sprintf("%1$s%2$s%3$s_demographics_finalsample_restOnly.csv", metric_output_dir, metric, atlas_network), row.names = F)
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
  write.csv(gam.age, sprintf("%1$sGAMresults.%2$s.age.%3$s%4$s_restOnly.csv", gam_dir, metric, atlas_name, network_parcellation), row.names = F, quote = F)
  return(gam.age)
}

  

## Makes smooths df for developmental trajectories, centered on zero
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
# @param dataset A character string of dataset (i.e. "NKI")

make_gam.smooths_centered <- function(metric, atlas, dataset){
  subxparcel_dem <- read.csv(sprintf("%1$ssensitivity_analyses/%2$s/%2$s%3$s_demographics_finalsample_restOnly.csv", outputs_root, metric, atlas))
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
  write.csv(smooth_fits, sprintf("%1$ssensitivity_analyses/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_centered.csv", outputs_root, metric, atlas))
  return(smooth_fits)
} 


################## 
# Load files
################## 

# load parcel labels
schaefer200.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")

# load parcellated S-A axis

schaefer200_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x7_SAaxis.csv")
schaefer200_SAaxis$label <- gsub("7Network", "Network", schaefer200_SAaxis$label) 
 
# load demographics
demographics <- read.csv(paste0(sample_selection_dir, dataset, "_demographics_finalsample_restOnly.csv"))
demographics <- demographics %>% dplyr::rename(subject=sub)
if (dataset=="HCPD") {
  demographics$age <- demographics$interview_age/12
} else if (dataset == "NKI") {
  demographics$subject <- gsub("sub-", "", demographics$subject)
}
demographics <- demographics %>% select(subject, sex, age, meanFD_avgSes)
demographics$sex <- as.factor(demographics$sex)

####################
# Fit GAMs for GBC #
####################
# create GBC.<atlas>.<dataset> dataframes and csv's
metric = "GBC"
df_name <- paste0(metric, ".", "schaefer200", ".", dataset)
assign(df_name, make_connMetricsDemog("schaefer200", metric, demographics, "subject", dataset))


# fit GAMs
gam.GBC.age.schaefer200 <- fitGAMs(schaefer200.parcel.labels, "GBC", "schaefer200", "", dataset)


# estimate centered GAM smooths 
smooth_fits_centered <- make_gam.smooths_centered(metric, "schaefer200", dataset)
print("GBC smooth_fits_centered made")
  