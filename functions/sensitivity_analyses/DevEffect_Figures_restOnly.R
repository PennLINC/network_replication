# code for making figures in 3_<dataset>_devEffectFigures_restOnly.Rmd
source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/GAM_functions.R")
library(dplyr)
library(magrittr)
library(ggplot2)
library(paletteer)

  
   

# load parcellated S-A axis 
schaefer200x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x17_SAaxis.csv")
schaefer200x17_SAaxis$label <- gsub("17Network", "Network", schaefer200x17_SAaxis$label) 
schaefer200_SAaxis <- schaefer200x17_SAaxis 
 


# load parcel labels  
schaefer200x17.parcel.labels <- read.csv("/cbica/projects/network_replication/atlases/parcellations/schaefer200x17_regionlist_final.csv")
schaefer200.parcel.labels <- schaefer200x17.parcel.labels #gbc
 



# Function for combining final GAM output dfs for figures
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name
combineGAMdfs <- function(metric, atlas){
  # Combine into Final Dfs
  gam_output <- get(paste0("gam.", metric, ".age.", atlas))
  SAaxis <- get(paste0(atlas, "_SAaxis"))
  df.list <- list(SAaxis, gam_output)
  axis <- Reduce(function(x,y) merge(x,y, all=TRUE, sort=F), df.list) 
  axis$SA.axis_rank <- as.numeric(axis$SA.axis_rank)
  axis <- axis[,-c(which(names(axis) =="X"))]

  return(axis)
}
  
# Function for calculating number of significant parcels (FDR corrected)
# @param axis A dataframe with SA-axis rank and GAM results
sig_parcels <- function(axis) {
  axis$Anova.age.pvalue.fdr <- p.adjust(axis$Anova.age.pvalue, method=c("fdr"))
  cat(sprintf("There are %s/%s significant parcels", sum(axis$Anova.age.pvalue.fdr < 0.05, na.rm=TRUE), nrow(axis)))
  
  axis$significant.fdr <- axis$Anova.age.pvalue.fdr < 0.05
  axis$significant.fdr[axis$significant.fdr == TRUE] <- 1
  axis$significant.fdr[axis$significant.fdr == FALSE] <- 0
  return(axis)
}
 

# Function for creating gam.smooths dfs
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name
make_gam.smooths <- function(metric, atlas, dataset){
  gam.smooths <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s.csv", dataset, metric, atlas))
  print("Smooths loaded")
  parcel.labels <- get(paste0(atlas, ".parcel.labels"))
  parcel.labels <- data.frame(parcel.labels$label)
  names(parcel.labels) <- "label"
  SAaxis <- get(paste0(atlas, "_SAaxis"))
  gam.smooths <- left_join(gam.smooths, parcel.labels, by = "label")

  
  a <- merge(parcel.labels, SAaxis, by="label", sort=F)
  
  gam.smooths <- left_join(gam.smooths, a, by = "label")
  gam.smooths$SA.axis_rank <-as.numeric(gam.smooths$SA.axis_rank) 
  return(gam.smooths)
}


# Function for creating gam.smooths dfs - covbat
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name
make_gam.smooths_covbat <- function(metric, atlas, dataset){
  gam.smooths <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_covbat.csv", dataset, metric, atlas))
  print("Smooths loaded")
  parcel.labels <- get(paste0(atlas, ".parcel.labels"))
  parcel.labels <- data.frame(parcel.labels$label)
  names(parcel.labels) <- "label"
  SAaxis <- get(paste0(atlas, "_SAaxis"))
  gam.smooths <- left_join(gam.smooths, parcel.labels, by = "label")
  
  
  a <- merge(parcel.labels, SAaxis, by="label", sort=F)
  
  gam.smooths <- left_join(gam.smooths, a, by = "label")
  gam.smooths$SA.axis_rank <-as.numeric(gam.smooths$SA.axis_rank) 
  return(gam.smooths)
}

 
 

## FIGURE: SA axis average rank vs. age effect - with spin test
# included all datapoints
# @param axis A dataframe with SA-axis rank and GAM results
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name
# @param r_text An expression including info on p_spin and r
# @param x_pos A numeric for horizontal position of r_text
# @param y_pos A numeric for vertical position of r_text

gam_figure_p.spin <- function(axis, metric, atlas, r_text, x_pos, y_pos){
  axis <- axis %>% mutate(SA.axis_rank_signif = ifelse(significant.fdr == 1, SA.axis_rank, NA))
  AgeEffect_figure <- ggplot(axis, aes(x=SA.axis_rank, y=GAM.age.AdjRsq, fill = SA.axis_rank_signif)) + 
    geom_point(color = "gray", shape = 21, size=3.5) + 
    paletteer::scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(axis$SA.axis_rank), max(axis$SA.axis_rank)), oob = squish) +
    paletteer::scale_color_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(axis$SA.axis_rank), max(axis$SA.axis_rank)), oob = squish) +
    labs(fill = "SA Axis Rank", x="\nSensorimotor-Association Axis Rank\n", y=expression(paste("Age Effect (Delta Adj", " R"^2, ")"))) +
     ylim(-0.07, 0.10) +  # ylim(-0.065, 0.15) for GBC; ylim(-0.07, 0.10) for BNC;  
    geom_smooth(data = axis, method='lm', se=TRUE, fill=alpha(c("gray70"),.9), col="black") +
    theme(
      axis.title.x=element_text(size=20, color = "black", margin = margin(t = -10, r = 0, b = 0, l = 0)),
      axis.title.y=element_text(size=20, color = "black"),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=20, color = "black"),
      panel.background=element_blank(),  
      legend.position = "none") +
    annotate(geom="text", x=x_pos, y=y_pos, label=r_text, color="black", size=7)
  return(AgeEffect_figure)
}

 
## Makes smooths df for developmental trajectories, centered on zero
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
# @param dataset A character string of dataset  

make_gam.smooths_centered <- function(metric, atlas, dataset){
  demog <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/%2$s%3$s_demographics_finalsample.csv", dataset, metric, atlas))
  demog$sex <- as.factor(demog$sex)
  #make wide_df to include parcels and age, sex, meanFD_avgSes 
  start_parcels <- which(names(demog) == "subject") +1
  end_parcels <- start_parcels + nrow(get(paste0(atlas, ".parcel.labels"))) - 1
  
  ind_age <- which(names(demog) == "age")
  ind_sex <- which(names(demog) == "sex")
  ind_meanFD_avgSes <- which(names(demog) == "meanFD_avgSes")
  wide_df <- demog[,c(start_parcels:end_parcels, ind_age, ind_sex, ind_meanFD_avgSes)]
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
  smooth_fits <- smooth_fits %>% mutate(SA.axis_percentile_bin = ntile(smooth_fits$SA.axis_rank, 10))
  write.csv(smooth_fits, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_centered.csv", dataset, metric, atlas))
  return(smooth_fits)
} 




## Makes smooths df for developmental trajectories, centered on zero - covbat
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
# @param dataset A character string of dataset  

make_gam.smooths_centered_covbat <- function(metric, atlas, dataset){
  demog <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample_covbat.csv", dataset, metric, atlas))
  demog$sex <- as.factor(demog$sex)
  demog <- demog[-c(which(names(demog) == "X"))]
  #make wide_df to include parcels and age, sex, meanFD_avgSes 
  start_parcels <- which(names(demog) == "subject") +1
  end_parcels <- start_parcels + nrow(get(paste0(atlas, ".parcel.labels"))) - 1

  ind_age <- which(names(demog) == "age")
  ind_sex <- which(names(demog) == "sex")
  ind_meanFD_avgSes <- which(names(demog) == "meanFD_avgSes")
  wide_df <- demog[,c(start_parcels:end_parcels, ind_age, ind_sex, ind_meanFD_avgSes)]
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
  smooth_fits <- smooth_fits %>% mutate(SA.axis_percentile_bin = ntile(smooth_fits$SA.axis_rank, 10))
  write.csv(smooth_fits, sprintf("/cbica/projects/network_replication/output/%1$s/sensitivity_analyses/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_centered_covbat.csv", dataset, metric, atlas))
  return(smooth_fits)
} 


## FIGURE: developmental trajectories, centered on zero
# @param atlas A character string of atlas name 
# @param metric A character string of connectivity metric (i.e. "GBC")
  
make_smooths_fig_centered <- function(atlas, metric){
  smooth_fits <- get(paste0("devTraj.centered.", atlas, ".", metric, "_df"))
  smooths_fig_centered <- ggplot(smooth_fits,aes(age,est,group=SA.axis_rank)) + 
    geom_line(data = smooth_fits, size=.7, alpha = .6, aes(color=SA.axis_rank)) + 
    ylim(-0.042, 0.055) + xlim(5, 23) +  
    paletteer::scale_color_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(smooth_fits$SA.axis_rank), max(smooth_fits$SA.axis_rank)), oob = squish) + 
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=20, color = "black"),
      panel.background=element_blank(),  
      legend.position = "none")  
  return(smooths_fig_centered)
}
  



