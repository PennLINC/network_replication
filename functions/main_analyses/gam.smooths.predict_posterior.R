library(nlme, lib.loc="/cbica/home/luoau/Rlibs")
library(tidyr, lib.loc="/cbica/home/luoau/Rlibs")
library(mgcv, lib.loc="/cbica/home/luoau/Rlibs")
library(gratia, lib.loc="/cbica/home/luoau/Rlibs")
library(dplyr, lib.loc="/cbica/home/luoau/Rlibs")



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