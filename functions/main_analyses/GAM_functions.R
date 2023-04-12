#library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)

#FIT GAM FUNCTION
##Function to fit a GAM (measure ~ s(smooth_var) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit <- function(measure, atlas, dataset, region, smooth_var, covariates, stats_only = FALSE){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
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


#FIT GAM FUNCTION with intercept
##Function to fit a GAM (measure ~ s(smooth_var) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit_intercept <- function(measure, atlas, dataset, region, smooth_var, covariates, stats_only = FALSE){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname) # Return the Value of a Named Object
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k=3, fx=T) + %s", region, smooth_var, covariates)) # 3 knots, fx = should term be unpenalized (so here, it is unpenalized)
  gam.model <- gam(modelformula, data=gam.data)
  gam.results <- summary(gam.model) # get intercept ( -> intercept)
  
  #Get the intercept
  gam.intercept <- gam.results$p.coeff["(Intercept)"]
  
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
  
  full.results <- cbind(parcel, gam.intercept, gam.smooth.F, gam.smooth.pvalue, gam.adjRsq, anova.smooth.pvalue, age.dev.onset, age.peak.change, max.age.increase, age.maturation)
  stats.results <- cbind(parcel, gam.intercept, gam.smooth.F, gam.smooth.pvalue, gam.adjRsq, anova.smooth.pvalue)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}

 

#### FIT GAM SMOOTH WITH A COVARIATE OF INTEREST ####
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariate of interest + control covariates)) and save out statistics for the first covariate
gam.fit.covariate <- function(measure, atlas, dataset, region, smooth_var, covariate.interest, covariates.noninterest, knots, set_fx = FALSE){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s + %s", region, smooth_var, knots, set_fx, covariate.interest, covariates.noninterest))
  gam.model <- gam(modelformula, method = "REML", data=gam.data)
  gam.results <- summary(gam.model)
  
  #GAM statistics
  #t-value for the covariate of interest term and GAM-based significance of this term
  gam.cov.tvalue <- gam.results$p.table[2,3]
  #GAM based significance of the term
  gam.cov.pvalue <- gam.results$p.table[2,4]
  
  #Calculate the magnitude and significance of the covariate of interest effect by comparing full and reduced models
  ##Compare a full model GAM (with the covariate of interst) to a nested, reduced model (without covariate of interst)
  nullmodel <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates.noninterest))
  gam.nullmodel <- gam(nullmodel, method = "REML", data=gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ##Full versus reduced model anova p-value
  anova.cov.pvalue <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  if(is.na(anova.cov.pvalue)){ #if residual deviance is exactly equal between full and reduced models and p=value = NA, set p = 1
    anova.cov.pvalue <- 1}
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  ### effect direction
  if(gam.cov.tvalue < 0){ #if the gam t-value for covariate of interest is less than 0, make the partialRsq negative
    partialRsq <- partialRsq*-1}
  
  results <- cbind(parcel, gam.cov.tvalue, gam.cov.pvalue, anova.cov.pvalue, partialRsq)
  return(results)
}




#PREDICT GAM SMOOTHS FUNCTION
##Function to fit GAM smooths based on model-predicted data
gam.predsmooth <- function(measure, atlas, dataset, region, smooth_var, covariates){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset)
  print(dataname)
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


#PREDICT GAM SMOOTHS FUNCTION
##Function to fit GAM smooths based on model-predicted data
gam.predsmooth_env <- function(gam.data, measure, atlas, region, smooth_var, covariates){
  
  #Fit the gam
  #print(dataname)
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



#### FIT GAM SMOOTH FUNCTION ####
##Function to fit a GAM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) per each region in atlas and save out statistics and derivative-based characteristics
gam.fit.smooth <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, stats_only = FALSE){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #GAM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to 0
  
  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  
  #Calculate the magnitude and significance of the smooth term effect by comparing full and reduced models
  ##Compare a full model GAM (with the smooth term) to a nested, reduced model (with covariates only)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, covariates)) #no smooth term
  gam.nullmodel <- gam(nullmodel, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  ##Full versus reduced model anova p-value
  anova.smooth.pvalue <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  ### effect direction
  mean.derivative <- mean(derv$derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    partialRsq <- partialRsq*-1}
  
  #Derivative-based temporal characteristics
  #Age of developmental change onset
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    change.onset <- min(derv$data[derv$sig==T])} #find first age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ #if gam derivative is never significant
    change.onset <- NA} #assign NA
  
  #Age of maximal developmental change
  if(sum(derv$sig) > 0){ 
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value significant derivatives
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    window.peak.change <- derv$data[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is greatest in absolute magnitude
    peak.change <- mean(window.peak.change)} #identify the age of peak developmental change
  if(sum(derv$sig) == 0){ 
    peak.change <- NA}  
  
  #Age of decrease onset
  if(sum(derv$sig) > 0){ 
    decreasing.range <- derv$data[derv$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
    if(length(decreasing.range) > 0)
      decrease.onset <- min(decreasing.range) #find youngest age with a significant negative derivative
    if(length(decreasing.range) == 0)
      decrease.onset <- NA}
  if(sum(derv$sig) == 0){
    decrease.onset <- NA}  
  
  #Age of increase offset
  if(sum(derv$sig) > 0){ 
    increasing.range <- derv$data[derv$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
    if(length(increasing.range) > 0)
      increase.offset <- max(increasing.range) #find oldest age with a significant positive derivative
    if(length(increasing.range) == 0)
      increase.offset <- NA}
  if(sum(derv$sig) == 0){ 
    increase.offset <- NA}  
  
  #Age of maturation
  if(sum(derv$sig) > 0){ 
    change.offset <- max(derv$data[derv$sig==T])} #find last age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ 
    change.offset <- NA}  
  
  full.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq, anova.smooth.pvalue, change.onset, peak.change, decrease.onset, increase.offset, change.offset)
  stats.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq, anova.smooth.pvalue)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}

#### PREDICT GAM SMOOTH FITTED VALUES FUNCTION ####
##Function to predict fitted values of a measure based on a fitted GAM smooth (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and a prediction df
gam.smooth.predict <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, increments){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset)  # e.g. GBC.schaefer200.PNC
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
  #Generate predictions based on the gam model and predication data frame
  predicted.smooth <- fitted_values(object = gam.model, data = pred)
  predicted.smooth <- predicted.smooth %>% select(all_of(smooth_var), fitted, se, lower, upper)
   
  smooth.fit <- list(parcel, predicted.smooth)
  return(smooth.fit)
}



#### PREDICT GAM SMOOTH FITTED VALUES FUNCTION with posterior ####
##Function to predict fitted values of a measure based on a fitted GAM smooth (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and a prediction df
# and for individual draws from the simulated posterior distribution
gam.smooth.predict_posterior <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_posterior_fits = TRUE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset)  # e.g. GBC.schaefer200.PNC
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
   
  #Generate predictions based on the gam model and predication data frame
  predicted.smooth <- fitted_values(object = gam.model, data = pred)
  predicted.smooth.fulldf <- predicted.smooth %>% select(all_of(smooth_var), fitted, se, lower, upper)
  predicted.smooth.fulldf <- predicted.smooth.fulldf %>% mutate(significant = !(0 > lower & 0 < upper))
  predicted.smooth.fulldf$significant.fit = predicted.smooth.fulldf$fitted*predicted.smooth.fulldf$significant
  colnames(predicted.smooth.fulldf) <- c(sprintf("%s", smooth_var), "fitted", "se", "lower", "upper", "significant", "significant.fit")
    
   
  #Estimate posterior fitted values from simulated GAM posterior distribution
  if(return_posterior_fits == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors for pred
       
    predicted.smooth.values <- X0 %*% t(sims) #Xp * simulated model coefficients = simulated posterior smooths Each column of predicted.smooth.values contains fitted values for a different draw from the simulated posterior distribution
    predicted.smooth.values <- as.data.frame(predicted.smooth.values)
    colnames(predicted.smooth.values) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    predicted.smooth.values <- cbind(as.numeric(pred[,smooth_var]), predicted.smooth.values) #add smooth_var increments from pred df to first column
    
    colnames(predicted.smooth.values)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    predicted.smooth.values <- as.data.frame(predicted.smooth.values) 
    predicted.smooth.values <- cbind(as.character(parcel), predicted.smooth.values) #add parcel label to first column
    colnames(predicted.smooth.values)[1] <- "label" #label the column
    predicted.smooth.values.long <- predicted.smooth.values %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.fits")
  } #np*npd rows, 3 columns (smooth_var, draw, posterior.fits)
  
  if(return_posterior_fits == FALSE)
    return(predicted.smooth.fulldf)
  if(return_posterior_fits == TRUE)
    return(predicted.smooth.values.long)
  
   
}

#### PREDICT GAM SMOOTH FITTED VALUES FOR A SPECIFIED VALUE OF AN INTERACTING COVARIATE ####
##Function to predict fitted values of a measure for a given value of a covariate, using a varying coefficients smooth-by-linear covariate interaction
gam.smooth.predict.covariateinteraction <- function(measure, atlas, dataset, region, smooth_var, int_var, int_var.predict, covariates, knots, set_fx = FALSE, increments){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s) + s(%2$s, by=%5$s, k=%3$s, fx=%4$s) + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  gam.model <- gam(modelformula, method = "REML", data=gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  pred[,int_var] <- as.numeric(int_var.predict)
  
  #Generate fitted (predicted) values based on the gam model and predication data frame
  predicted.smooth <- fitted_values(object = gam.model, data = pred)
  predicted.smooth$fitted.centered <- (predicted.smooth$fitted-gam.results$p.table[1,1]) #subtract the intercept from fitted values
  predicted.smooth <- predicted.smooth %>% select(all_of(smooth_var), fitted, se, lower, upper, fitted.centered)
  
  return(predicted.smooth)
}

#### PREDICT GAM SMOOTH FITTED VALUES FUNCTION with posterior - for envSES analyses ####
##Function to predict fitted values of a measure based on a fitted GAM smooth (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates)) and a prediction df
# and for individual draws from the simulated posterior distribution
gam.smooth.predict_posterior_env <- function(gam.data, measure, atlas, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_posterior_fits = TRUE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  
  #Fit the gam

  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
  #Generate predictions based on the gam model and predication data frame
  predicted.smooth <- fitted_values(object = gam.model, data = pred)
  predicted.smooth.fulldf <- predicted.smooth %>% select(all_of(smooth_var), fitted, se, lower, upper)
  predicted.smooth.fulldf <- predicted.smooth.fulldf %>% mutate(significant = !(0 > lower & 0 < upper))
  predicted.smooth.fulldf$significant.fit = predicted.smooth.fulldf$fitted*predicted.smooth.fulldf$significant
  colnames(predicted.smooth.fulldf) <- c(sprintf("%s", smooth_var), "fitted", "se", "lower", "upper", "significant", "significant.fit")
  
  
  #Estimate posterior fitted values from simulated GAM posterior distribution
  if(return_posterior_fits == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors for pred
    
    predicted.smooth.values <- X0 %*% t(sims) #Xp * simulated model coefficients = simulated posterior smooths Each column of predicted.smooth.values contains fitted values for a different draw from the simulated posterior distribution
    predicted.smooth.values <- as.data.frame(predicted.smooth.values)
    colnames(predicted.smooth.values) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    predicted.smooth.values <- cbind(as.numeric(pred[,smooth_var]), predicted.smooth.values) #add smooth_var increments from pred df to first column
    
    colnames(predicted.smooth.values)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    predicted.smooth.values <- as.data.frame(predicted.smooth.values) 
    predicted.smooth.values <- cbind(as.character(parcel), predicted.smooth.values) #add parcel label to first column
    colnames(predicted.smooth.values)[1] <- "label" #label the column
    predicted.smooth.values.long <- predicted.smooth.values %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.fits")
  } #np*npd rows, 3 columns (smooth_var, draw, posterior.fits)
  
  if(return_posterior_fits == FALSE)
    return(predicted.smooth.fulldf)
  if(return_posterior_fits == TRUE)
    return(predicted.smooth.values.long)
  
  
}
#### CALCULATE SMOOTH ESTIMATES FUNCTION ####
##Function to estimate the zero-averaged gam smooth function 
gam.estimate.smooth <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, increments){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
  #Estimate the smooth trajectory 
  estimated.smooth <- smooth_estimates(object = gam.model, data = pred)
  estimated.smooth <- estimated.smooth %>% select(age, est)
  
  return(estimated.smooth)
}

#### POSTERIOR DISTRIBUTION SMOOTHS FUNCTION ####
##Function to simulate the posterior distribution from a fitted GAM, calculate smooths for individual posterior draws, and return smooth max and min values + 95% credible intervals
gam.posterior.smooths <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_draws = TRUE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior smooths estimated
  np <- as.numeric(increments) #number of smooth_var increments to predict fit at from posterior model coefficients
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame, used to estimate posterior model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to fit for each posterior draw
  for (v in c(1:length(theseVars))) { #fill the prediction df with data 
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
  #Estimate posterior smooth functions (fitted values) from simulated GAM posterior distribution  
  ##Each of the posterior draws has a fitted spline (+ intercept + covariate coefficients) that includes the uncertainty in the estimated model coefficients  
  ##Akin to fitted_samples from gratia  
  Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
  sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance. 
  X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors that maps model parameters to the smooth fit (outcome measure scale)
  predicted.smooth.values <- X0 %*% t(sims) #generate posterior smooths (fitted y for each set of posterior draw model parameters)
  colnames(predicted.smooth.values) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
  predicted.smooth.values <- cbind(as.numeric(pred[,smooth_var]), predicted.smooth.values) #add smooth_var increments from pred df to first column
  colnames(predicted.smooth.values)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
  predicted.smooth.values <- as.data.frame(predicted.smooth.values) 
  
  #Smooth minimum/maximum values and credible intervals
  #Smooth max + 95% credible interval
  max.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is largest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x){
                       round(predicted.smooth.values[,smooth_var][which.max(x)],2)
                     }))
  max.y.range <- t(max.y.range)
  max.y <- median(max.y.range) #median value 
  max.y.CI <- quantile(max.y.range, probs = c(0.025, 0.975)) #credible interval
  max.y.CI.lower <- max.y.CI[[1]]
  max.y.CI.upper <- max.y.CI[[2]]
  #Smooth min + 95% credible interval
  min.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is lowest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x){
                       round(predicted.smooth.values[,smooth_var][which.min(x)],2)
                     }))  
  min.y.range <- t(min.y.range)
  min.y <- median(min.y.range) #median value
  min.y.CI <- quantile(min.y.range, probs = c(0.025, 0.975)) #credible interval
  min.y.CI.lower <- min.y.CI[[1]]
  min.y.CI.upper <- min.y.CI[[2]]
  
  if(return_draws == TRUE)
    return(predicted.smooth.values)
  if(return_draws == FALSE)
    smooth.features <- list(parcel, max.y, max.y.CI.lower, max.y.CI.upper, min.y, min.y.CI.lower, min.y.CI.upper)
  names(smooth.features) <- c("parcel", sprintf("%s at max y", smooth_var), "max y credible interval lower", "max y credible interval upper", sprintf("%s at min y", smooth_var), "min y credible interval lower", "min y credible interval upper")
  return(smooth.features)
}

#### DERIVATIVES FUNCTION ####
##Function to compute smooth derivatives for a main GAM model and for individual draws from the simulated posterior distribution
gam.derivatives <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_posterior_derivatives = TRUE){
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame, used to estimate posterior model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to fit for each posterior draw
  for (v in c(1:length(theseVars))) { #fill the prediction df with data 
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% select(-init) #prediction df
  pred2 <- pred #second prediction df
  pred2[,smooth_var] <- pred[,smooth_var] + EPS #finite differences
  
  #Estimate smooth derivatives
  derivs <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = UNCONDITIONAL, newdata = pred) #derivative at 200 indices of smooth_var with a simultaneous CI
  derivs.fulldf <- derivs %>% select(data, derivative, se, lower, upper)
  derivs.fulldf <- derivs.fulldf %>% mutate(significant = !(0 > lower & 0 < upper))
  derivs.fulldf$significant.derivative = derivs.fulldf$derivative*derivs.fulldf$significant
  colnames(derivs.fulldf) <- c(sprintf("%s", smooth_var), "derivative", "se", "lower", "upper", "significant", "significant.derivative")
  
  #Estimate posterior smooth derivatives from simulated GAM posterior distribution
  if(return_posterior_derivatives == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
    sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors for pred
    X1 <- predict(gam.model, newdata = pred2, type = "lpmatrix") #get matrix of linear predictors for pred2
    Xp <- (X1 - X0) / EPS 
    posterior.derivs <- Xp %*% t(sims) #Xp * simulated model coefficients = simulated derivatives. Each column of posterior.derivs contains derivatives for a different draw from the simulated posterior distribution
    posterior.derivs <- as.data.frame(posterior.derivs)
    colnames(posterior.derivs) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    posterior.derivs <- cbind(as.numeric(pred[,smooth_var]), posterior.derivs) #add smooth_var increments from pred df to first column
    colnames(posterior.derivs)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    posterior.derivs <- cbind(as.character(parcel), posterior.derivs) #add parcel label to first column
    colnames(posterior.derivs)[1] <- "label" #label the column
    posterior.derivs.long <- posterior.derivs %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.derivative")
  } #np*npd rows, 3 columns (smooth_var, draw, posterior.derivative)
  
  if(return_posterior_derivatives == FALSE)
    return(derivs.fulldf)
  if(return_posterior_derivatives == TRUE)
    return(posterior.derivs.long)
}
 



