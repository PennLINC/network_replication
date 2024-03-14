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
metric_output_dir <- paste0(outputs_root, "sensitivity_analyses/networkpair/")
gam_dir <- paste0(config$outputs_root, "sensitivity_analyses/networkpair/GAM/")

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

#FIT GAM FUNCTION - for nonGSR, absvalue, and thresholded (and network pairs)
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



#PREDICT GAM SMOOTHS FUNCTION
##Function to fit GAM smooths based on model-predicted data
gam.predsmooth <- function(measure, atlas, conn_type, region, smooth_var, covariates){
  
  #Fit the gam
  dataname <- sprintf("%s.%s_%s", measure, atlas, conn_type) 
  gam.data <- get(dataname)  
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k=3, fx=T) + %s", region, smooth_var, covariates))
  model <- gam(modelformula, data=gam.data)
  
  #Extract model summary and input data
  s <- summary(model) #model summary
  df <- model$model #extract df of y + terms, for all participants (i.e., the data used to build the model)
  
  #Create a dataframe for predictions
  np <- 1000 #number of predictions to make
  thisPred <- data.frame(init = rep(0,np)) #initiate a np matrix
  
  theseVars <- attr(model$terms,"term.labels") #the GAM terms
  varClasses <- attr(model$terms,"dataClasses") #classes of the GAM terms
  thisResp <- as.character(model$terms[[2]]) #the y 
  for (v in c(1:length(theseVars))) { #create "thisPred" df with data for predictions. This will be used to predict the y variable based on different values of the smooth_var, holding all other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #seq generates a sequence of np data points, starting from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
  #Generate model-based predictions for y based on smooth_var values 
  p <- data.frame(predict(model,pred,se.fit = T)) #using the GAM model, predict the response variable y (and the SE of the prediction) based on prediction df generated above
  pred <- cbind(pred,p)
  pred$selo <- pred$fit - 2*pred$se.fit #lower bound of prediction interval
  pred$sehi <- pred$fit + 2*pred$se.fit #upper bound of prediction interval
  
  predicted.smooth <- pred %>% select(smooth_var, fit, se.fit, selo, sehi)
  
  #Get the point at which y is maximal, based on the predicted smooth
  maxval <- max(predicted.smooth$fit)
  peak <- predicted.smooth[,smooth_var][predicted.smooth$fit == maxval]
  
  smooth.fit <- list(parcel, peak, predicted.smooth)
  return(smooth.fit)
}

################## 
# Load files
################## 

# load demographics
demographics <- read.csv(paste0(sample_selection_dir, dataset, "_demographics_finalsample.csv"))
demographics <- demographics %>% dplyr::rename(subject=sub)
if (dataset=="HCPD") {
  demographics$age <- demographics$interview_age/12
} else if (dataset == "NKI") {
  demographics$subject <- gsub("sub-", "", demographics$subject)
}
demographics <- demographics %>% select(subject, sex, age, meanFD_avgSes)


# load mean connectivity - combine with demographics
subxnetpair.matrix <- read.csv(sprintf("%1$snetworkpair_subxnetpair_matrix_schaefer200x7_orig.csv", metric_output_dir))
netpair.schaefer200_orig <- merge(subxnetpair.matrix, demographics, by="subject")
netpair.schaefer200_orig$sex <- as.factor(netpair.schaefer200_orig$sex)


# extract network-network labels 
netpair_labels <- names(subxnetpair.matrix)[c(2:length(names(subxnetpair.matrix)))]


##############################
# Fit GAMs for Network Pairs #
##############################

# network connectivity GAMs 
gam.age <- matrix(data=NA, nrow=length(netpair_labels), ncol=9) #empty matrix to save gam.fit output to

for(row in c(1:length(netpair_labels))){ #for each region
  netpair <- netpair_labels[row] 
  GAM.RESULTS <- gam.fit(measure = "netpair", atlas = "schaefer200", conn_type = "orig", region = netpair, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.fit function
  gam.age[row,] <- GAM.RESULTS}
gam.age <- as.data.frame(gam.age)
colnames(gam.age) <- c("label","GAM.age.Fvalue","GAM.age.pvalue","GAM.age.AdjRsq","Anova.age.pvalue","age.onsetchange","age.peakchange","maxage.BNC.increase","age.maturation" )
cols = c(2:9)    
gam.age[,cols] = apply(gam.age[,cols], 2, function(x) as.numeric(as.character(x)))
gam.age <- gam.age %>% mutate(significant = (Anova.age.pvalue < 0.05))

write.csv(gam.age, sprintf("%1$sGAMresults.networkpair.age.csv", gam_dir), row.names = F, quote = F)




# Estimate GAM smooths based on model-predicted data and save out predicted y data**   
gam.smooths <- matrix(data=NA, ncol=7) #empty matrix to save gam.predsmooth fits to
colnames(gam.smooths) <- c("age","fit","se.fit","selo","sehi","index","label")

for(row in c(1:length(netpair_labels))){ #for each region
  netpair <- netpair_labels[row] 
  GAM.SMOOTH <- gam.predsmooth(measure = "netpair", atlas = "schaefer200", conn_type = "orig", region = netpair, smooth_var = "age", covariates = "sex + meanFD_avgSes") #run the gam.predsmooth function
  
  preddata <- as.data.frame(GAM.SMOOTH[3]) #get predicted.smooth df from function output
  preddata$index <- rep(x=row, 1000) #region index
  preddata$label <- rep(x=GAM.SMOOTH[1], 1000) #label
  gam.smooths <- rbind(gam.smooths, preddata)
  
  
}
gam.smooths <- gam.smooths[-1,] #remove empty initialization row
gam.smooths$label <- as.character(gam.smooths$label)


write.csv(gam.smooths, sprintf("%1$sGAMsmoothfits.networkpair.csv", gam_dir), row.names = F, quote = F)

 