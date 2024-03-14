
library(dplyr)
library(gratia)
library(magrittr)
library(mgcv)
library(rjson)
library(stringr)

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
gam_dir <- paste0(config$outputs_root, "sensitivity_analyses/GBC/GAM/")

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

#FIT GAM FUNCTION - for nonGSR, absvalue, and thresholded
##Function to fit a GAM (measure ~ s(smooth_var) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit <- function(measure, atlas, conn_type, region, smooth_var, covariates, stats_only = FALSE){
  
  #Fit the gam
  dataname <- sprintf("%s.%s_%s", measure, atlas, conn_type) 
  gam.data <- get(dataname) # Return the Value of a Named Object
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k=3, fx=T) + %s", region, smooth_var, covariates)) # 3 knots, fx = should term be unpenalized (so here, it is unpenalized)
  gam.model <- gam(modelformula, data=gam.data)
  gam.results <- summary(gam.model) 
  
  #Get the F value for the smooth term and the GAM-based significance of the term
  # (F-tests on smooth terms (rather than for the full model) are 
  # joint tests for equality to zero for all of the coefficients making up a single spline term)
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  
  #Calculate the magnitude and significance of the **smooth term** effect based on delta adjusted R^2
  ##Compare a full model GAM (with the smooth term) to a nested, reduced model (with covariates only)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, covariates))
  gam.nullmodel <- gam(nullmodel, data=gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  anova.smooth.pvalue <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2] #anova between reduced and full models
  
  gam.adjRsq <- abs(gam.results$r.sq - gam.nullmodel.results$r.sq) #delta R.sq (adj)
  
  #to get the sign of the effect, fit a linear model and extract the sign
  linearmodel <- as.formula(sprintf("%s ~ %s + %s", region, smooth_var, covariates))
  lm.model.t <- summary(lm(linearmodel, data=gam.data))$coefficients[2,3] #t-value for smooth_var
  if(lm.model.t < 0){ #if the linear model t-value for smooth_var is less than 0, make the delta adj R.sq negative
    gam.adjRsq <- gam.adjRsq*-1}
  
  #Get derivatives of the smooth function using the gratia derivatives function; gratia estimates derivatives from GAM smooths via finite differences
  derv <- derivatives(gam.model,term=sprintf('s(%s)',smooth_var))
  
  #Identify derivative significance window (basically finding where there is significant change in a region)
  derv <- derv %>%
    mutate(sig = !(0 >lower & 0 < upper)) #add "sig" column (TRUE/FALSE) to derv. Derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., in the CI does not include 0).
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv column that has all the significant derivative values, but all non-significant derivatives are set to 0. Aka, mask out non-significant deriv points by the derv$sig TRUE/FALSE column
  
  #Age of developmental change onset
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    age.dev.onset <- min(derv$data[derv$sig==T])} #find minimum age (stored in derv$data) where sig = T
  if(sum(derv$sig) == 0){ #if gam derivate is not significant, assign NA
    age.dev.onset <- NA}  
  
  #Age of peak development
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value the sig_deriv column 
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    ages.peak.change <- derv$data[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is largest
    age.peak.change <- mean(ages.peak.change)} #identify the age of peak developmental change. If multiple ages have the same change rate, identify the average age
  if(sum(derv$sig) == 0){ #if gam is not significant, assign NA
    age.peak.change <- NA}  
  
  #Oldest age of significant increase
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    increasing.agerange <- derv$data[derv$sig_deriv > 0] #identify ages where measure is significantly increasing (positive derivative)
    if(length(increasing.agerange) > 0)
      max.age.increase <- max(increasing.agerange) #identify max age
    if(length(increasing.agerange) == 0)
      max.age.increase <- NA}
  if(sum(derv$sig) == 0){ #if gam is not significant, assign NA
    max.age.increase <- NA}  
  
  #Age of maturation
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    age.maturation <- max(derv$data[derv$sig==T])} #maximum age where sig = T
  if(sum(derv$sig) == 0){ #if gam is not significant, assign NA
    age.maturation <- NA}  
  
  full.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, gam.adjRsq, anova.smooth.pvalue, age.dev.onset, age.peak.change, max.age.increase, age.maturation)
  stats.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, gam.adjRsq, anova.smooth.pvalue)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}
 

# Function for **Fitting GAMs** for GBC (absolute value and thresholded)
# Fit GAM (func_conn_metric ~ s(age) + sex + meanFD_avgSes)) per each region in atlas and save out statistics and derivative-based characteristics
# @param parcel.labels A vector of parcel labels
# @param metric A character string of connectivity metric  
# @param atlas_name A character string of atlas  
# @param network_parcellation A character string for schaefer parcellation ("" for non-schaefer atlas; "x7" or "x17" for schaefer atlas)
# @param dataset_name A character string of dataset (e.g. "PNC")
# @param conn_type A character string, which connectivity type (e.g. "absvalue" or "thresholded)
fitGAMs <- function(parcel.labels, metric, atlas_name, network_parcellation, dataset_name, conn_type) {
  parcel.labels <- parcel.labels$label
  gam.age <- matrix(data=NA, nrow=length(parcel.labels), ncol=9) #empty matrix to save gam.fit output to
  
  for(row in c(1:length(parcel.labels))){ #for each region
    region <- parcel.labels[row] 
    GAM.RESULTS <- gam.fit(measure = metric, atlas = paste0(atlas_name, network_parcellation), conn_type = conn_type, region = region, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
    gam.age[row,] <- GAM.RESULTS}
  gam.age <- as.data.frame(gam.age)
  colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
  cols = c(2:9)    
  gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
  gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))
  write.csv(gam.age, sprintf("%1$sGAMresults.%2$s.age.%3$s%4$s_%5$s.csv", gam_dir, metric, atlas_name, network_parcellation, conn_type), row.names = F, quote = F)
  return(gam.age)
}

################## 
# Load files
################## 

# load parcel labels
schaefer200.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x7_regionlist_final.csv")
 
# load demographics
demographics <- read.csv("/cbica/projects/network_replication/revisions/input/PNC/CPAC_sample_selection/CPAC_nonGSR_PNC_demographics_finalsample_20231130.csv")
demographics <- demographics %>% dplyr::rename(subject=sub)
subxparcel.matrix <- read.csv("/cbica/projects/network_replication/revisions/output/PNC/GBC/GBC_subxparcel_matrix_nonGSR.csv")
GBC.schaefer200_nonGSR <- merge(subxparcel.matrix, demographics, by="subject")
GBC.schaefer200_nonGSR$sex <- as.factor(GBC.schaefer200_nonGSR$sex) 

####################
# Fit GAMs for GBC #
####################

# create GBC.<atlas>.<dataset>.<conn_type> dataframes and csv's
metric = "GBC"
 
# fit GAMs 
gam.GBC.age.schaefer200.nonGSR <- fitGAMs(schaefer200.parcel.labels, "GBC", "schaefer200", "", dataset, "nonGSR")
print("GAMs fit for nonGSR")

 


