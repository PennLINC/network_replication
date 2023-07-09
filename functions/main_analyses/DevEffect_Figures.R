# code for making figures in 3_<dataset>_devEffectFigures.Rmd
source("/cbica/projects/network_replication/Rscripts/functions/main_analyses/GAM_functions.R")
library(dplyr)
library(magrittr)
library(ggplot2)
library(paletteer)
 

# load parcellated S-A axis
glasser_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/glasser_SAaxis.csv")
gordon_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/gordon_SAaxis.csv")

schaefer200x7_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x7_SAaxis.csv")
schaefer200x7_SAaxis$label <- gsub("7Network", "Network", schaefer200x7_SAaxis$label) 

schaefer200x17_SAaxis <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer200x17_SAaxis.csv")
schaefer200x17_SAaxis$label <- gsub("17Network", "Network", schaefer200x17_SAaxis$label) 
schaefer200_SAaxis <- schaefer200x17_SAaxis 

schaefer400x7_SAaxis  <- read.csv("/cbica/projects/network_replication/SAaxis/schaefer400x7_SAaxis.csv")
schaefer400x7_SAaxis$label <-  gsub("7Networks", "Networks", schaefer400x7_SAaxis$label)

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
  gam.smooths <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s.csv", dataset, metric, atlas))
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
  gam.smooths <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_covbat.csv", dataset, metric, atlas))
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
 
gam_figure_p.spin <- function(axis, metric, atlas, r_text, x_pos, y_pos, ylim_lower, ylim_upper){
  axis <- axis %>% mutate(SA.axis_rank_signif = ifelse(significant.fdr == 1, SA.axis_rank, NA))
  AgeEffect_figure <- ggplot(axis, aes(x=SA.axis_rank, y=GAM.age.AdjRsq, fill = SA.axis_rank_signif)) + 
    geom_point(color = "gray", shape = 21, size=3.5) + 
    paletteer::scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(axis$SA.axis_rank), max(axis$SA.axis_rank)), oob = squish) +
    paletteer::scale_color_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(axis$SA.axis_rank), max(axis$SA.axis_rank)), oob = squish) +
      labs(fill = "SA Axis Rank", x="\nSensorimotor-Association Axis Rank\n", y=expression(paste("Age Effect (Delta Adj", " R"^2, ")"))) +
      #ggtitle(paste(metric, "-", atlas)) + 
    ylim(ylim_lower, ylim_upper) +  # ylim(-0.065, 0.15) for GBC;  ylim(-0.07, 0.10) for BNC;  ylim(-0.055, 0.17) for WNC; 
      geom_smooth(data = axis, method='lm', se=TRUE, fill=alpha(c("gray70"),.9), col="black") +
      theme(
        axis.title.x=element_text(size=24, color = "black", margin = margin(t = -15, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=24, color = "black"),
        axis.line = element_line(color = "black"),
        axis.text=element_text(size=24, color = "black"),
        panel.background=element_blank(),  
        legend.position = "none") +
     annotate(geom="text", x=x_pos, y=y_pos, label=r_text, color="black", size=7)
  return(AgeEffect_figure)
}
 


## Makes smooths df for developmental trajectories, centered on zero
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
# @param dataset A character string of dataset (i.e. "NKI")

make_gam.smooths_centered <- function(metric, atlas, dataset){
  demog <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample.csv", dataset, metric, atlas))
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
  write.csv(smooth_fits, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_centered.csv", dataset, metric, atlas))
  return(smooth_fits)
} 



## Makes smooths df for developmental trajectories, centered on zero - covbat
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
# @param dataset A character string of dataset (i.e. "NKI")

make_gam.smooths_centered_covbat <- function(metric, atlas, dataset){
  demog <- read.csv(sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/%2$s%3$s_demographics_finalsample_covbat.csv", dataset, metric, atlas))
  if (dataset=="NKI"){
    demog$sex <- as.factor(demog$gender)
  } else {
    demog$sex <- as.factor(demog$sex)
  }
  if (dataset=="HBN") {
    demog <- demog[-c(which(names(demog) == "X"))]
  }
  #make wide_df to include parcels and age, sex, meanFD_avgSes 
  start_parcels <- which(names(demog) == "subject") +1
  if (dataset=="HCPD") {
    demog <- demog[,-2]
    end_parcels <- start_parcels + nrow(get(paste0(atlas, ".parcel.labels"))) - 1
  } else if (dataset=="HBN") {
    end_parcels <- start_parcels + nrow(get(paste0(atlas, ".parcel.labels"))) - 1
  } else {
    end_parcels <- which(names(demog) == "X")-1
  }
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
  write.csv(smooth_fits, sprintf("/cbica/projects/network_replication/output/%1$s/%2$s/GAM/GAMsmoothfits.%2$s.age.%3$s_centered_covbat.csv", dataset, metric, atlas))
  return(smooth_fits)
} 


## FIGURE: developmental trajectories, centered on zero
# @param atlas A character string of atlas name 
# @param metric A character string of connectivity metric (i.e. "GBC")
  
make_smooths_fig_centered <- function(atlas, metric){
  smooth_fits <- get(paste0("devTraj.centered.", atlas, ".", metric, "_df"))
  smooths_fig_centered <- ggplot(smooth_fits,aes(age,est,group=SA.axis_rank)) + 
    geom_line(data = smooth_fits, size=.7, alpha = .6, aes(color=SA.axis_rank)) + 
    ylim(-0.035, 0.055) + 
    #ylab("Global Brain Connectivity \n(Zero-Centered)") + xlab("Age") + 
    paletteer::scale_color_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(smooth_fits$SA.axis_rank), max(smooth_fits$SA.axis_rank)), oob = squish) + 
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      #axis.title.x=element_text(size=24, color = "black"),
      #axis.title.y=element_text(size=24, color = "black"),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=32, color = "black"),
      panel.background=element_blank(),  
      legend.position = "none")  
  return(smooths_fig_centered)
}
 
 

## FIGURE: developmental trajectories for representative parcels
# @param atlas A character string of atlas name 
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param parcels A character string of parcel(s) whose trajectory you want to plot (i.e. "SomMot", "visual")
# @param color_hexcode A character string for color of geom_line
# @param SA_rank A numeric for mean SA rank of parcels
# @param ggseg_atlas A character for 
make_repParcel.fig <- function(metric, atlas, parcels, color_hexcode) {
  df <- get(paste0("repParcels_", atlas))
  smooth_fits <- df[[parcels]]
  main.plot <- ggplot(smooth_fits,aes(age,est,group=SA.axis_rank)) + 
    geom_line(data = smooth_fits, size=1.5, alpha = .6,  colour=color_hexcode) + ylim(-0.035, 0.055) + xlim(min(smooth_fits$age), max(smooth_fits$age)) + 
    #ylab("Global Brain Connectivity \n(Zero-Centered)") + xlab("Age") + 
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      #axis.title.x=element_text(size=24, color = "black"),
      #axis.title.y=element_text(size=24, color = "black"),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=32, color = "black"),
      panel.background=element_blank(),  
      legend.position = "none")  
  #plot.with.inset <- main.plot + annotation_custom(ggplotGrob(fig_inset), xmin = 16, xmax = 20, ymin = ggplot_coords[3], ymax = ggplot_coords[4])
  return(main.plot)
}



## FIGURE: inset for developmental trajectories for representative parcels
# @param atlas A character string of atlas name 
# @param parcels A character string of parcel(s) whose trajectory you want to plot (i.e. "SomMot", "visual")
# @param color_hexcode A character string for color of geom_line
# @param SA_rank A numeric for mean SA rank of parcels
# @param ggseg_atlas A character for 
make_repParcel.inset <- function(atlas, parcels, color_hexcode, SA_rank, ggseg_atlas) { 
  SAaxis_fig.df <- get(paste0(atlas, "_SAaxis_fig.df")) 
  names(SAaxis_fig.df)[which(names(SAaxis_fig.df) == parcels)] <- "repParcel"
  fig_inset <- ggplot() + geom_brain(data=SAaxis_fig.df, 
                                     atlas=get(ggseg_atlas), 
                                     hemi='left', side = 'lateral',
                                     mapping=aes(fill=repParcel, colour=I(color_hexcode),  size=I(0.6)), 
                                     show.legend=FALSE) + scale_fill_gradient(low="white", high = color_hexcode) + theme_void() +
    labs(title=paste("Mean S-A Rank =", SA_rank), color="black", fontface="bold")
  return(fig_inset)
  }


 
## FIGURE: Fitted GBC at each age -
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas_name A character string of atlas name 
plot_cred.intervals <- function(atlas_name, metric) {
  cor.credible.intervals <- get(paste0("cor.credible.intervals_", metric)) 
  cor.credible.intervals <- cor.credible.intervals[[atlas_name]]
  cor.credible.intervals$max.cor.window
  fig <- ggplot(cor.credible.intervals, aes(x = age, y = SA.correlation, ymin = lower, ymax = upper)) + 
    geom_line(size = .7) +
    geom_ribbon(alpha = .2, fill = c("grey60")) +
    labs(x="Age", y="Spatial Alignment of \nFitted GBC to S-A Axis (r)\n") + ggtitle(paste0(metric, "-", atlas_name, ": ", dataset)) +
    theme_classic() + ylim(-0.61, 0) +  
    theme(axis.text = element_text(size = 20, 
                                   color = c("black")), 
          axis.title = element_text(size = 20, color = c("black")), 
          axis.line = element_line(size =.22), 
          axis.ticks = element_line(size = .22)) +
    scale_x_continuous(breaks=c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22), expand = c(0, .1)) 
  return(fig)
}


##########################
######### Extras ######### 
##########################

# Functions for making figures
# load gordon parcel xcp mapping
gordon_mapping <- read.csv("/cbica/projects/network_replication/atlases/parcellations/Gordon333Atlas.Parcels.LabelKey.csv")
## FIGURE: Surface View of Age Effect
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name
# NOTE: glasser and schaefer17_400 need to have separate hemispheres; gordon needs labels to be re-mapped

ageEffect_surf_figure <- function(atlas, metric){
  if (metric=="GBC" & str_detect(atlas,"schaefer200")) {
    ggseg_atlas <- paste0("schaefer17_200")
  } else if (metric=="GBC" & str_detect(atlas,"schaefer400")){
    ggseg_atlas <- paste0("schaefer17_400")  
  } else if (str_detect(atlas, "x")) {
    network <- str_extract(atlas, "x[0-9]+")[1]
    network <- gsub("x", "", network)
    parcellation <- str_extract(atlas, "r[0-9]{3}")[1]
    parcellation <- gsub("r", "", parcellation)
    ggseg_atlas <- paste0("schaefer", network, "_", parcellation)
  } else if (atlas == "gordon" | atlas=="glasser") {
    ggseg_atlas <- atlas
  } else (
    print("Please provide valid atlas.")
  )
  df.dev <- paste0(metric, ".axis_", atlas)
  df.dev <- get(df.dev)
  df.dev <- df.dev %>% mutate(SA.axis_rank_signif = ifelse(significant.fdr == 1, SA.axis_rank, NA))
  df.dev$region <- df.dev$label
  df.dev <- data.frame(cbind(df.dev$region, df.dev$GAM.age.AdjRsq))
  names(df.dev) <- c("region", "GAM.age.AdjRsq")
  df.dev$GAM.age.AdjRsq <- as.numeric(df.dev$GAM.age.AdjRsq)
  # adjust region names
  if(atlas =="schaefer200" | atlas =="schaefer400" | str_detect(atlas, "17")) {
    df.dev$region <- gsub("Network", "17Network", df.dev$region)
  } else if(str_detect(atlas, "7")) {
    df.dev$region <- gsub("Network", "7Network", df.dev$region)
  } else if(atlas=="gordon"){
    df.dev$region <- gordon_mapping$parcel.balsa
    df.dev$region <- gsub("L_", "", df.dev$region)
    df.dev$region <- gsub("R_", "", df.dev$region)
  } else if(atlas=="glasser") {
    df.dev$region <- df.dev$region  # nothing needs to be done for glasser
  } else{
    print("Please provide valid atlas.")
  }
  
  # split df's into right and left hemisphere for schaefer400x17 and glasser 
  # make ggseg figure
  if(atlas=="glasser") { 
    r_df.dev <-  df.dev
    r_df.dev$region<- gsub(x = r_df.dev$region, pattern = "R_", replacement = "")
    r_df.dev$region <-  gsub(x = r_df.dev$region, pattern = "_ROI", replacement = "")
    r_df.dev <- r_df.dev[-which(grepl("L_", r_df.dev$region)),]
    
    l_df.dev <-  df.dev
    l_df.dev$region<- gsub(x = l_df.dev$region, pattern = "L_", replacement = "")
    l_df.dev$region <-  gsub(x = l_df.dev$region, pattern = "_ROI", replacement = "")
    l_df.dev <- l_df.dev[-which(grepl("R_", l_df.dev$region)),]
    
    
    r_glasser_ageEffect <- ggplot() + geom_brain(data=r_df.dev, 
                                                 atlas=get(ggseg_atlas), 
                                                 mapping=aes(fill=GAM.age.AdjRsq), 
                                                 show.legend=TRUE, 
                                                 hemi = "right") + 
      scale_fill_gradientn(colors= c("#00A3A7FF", "#FFFFFF","#C75DAAFF"), 
                           limits = range(df.dev$GAM.age.AdjRsq), 
                           values=rescale(c(min(df.dev$GAM.age.AdjRsq),0,max(df.dev$GAM.age.AdjRsq)))) +
      theme_void()
    
    l_glasser_ageEffect <- ggplot() + geom_brain(data=l_df.dev, 
                                                 atlas=glasser, 
                                                 mapping=aes(fill=GAM.age.AdjRsq), 
                                                 show.legend=TRUE, 
                                                 hemi = "left") + 
      labs(title = paste(atlas, metric)) +
      scale_fill_gradientn(colors= c("#00A3A7FF", "#FFFFFF","#C75DAAFF"), 
                           limits = range(df.dev$GAM.age.AdjRsq), 
                           values=rescale(c(min(df.dev$GAM.age.AdjRsq),0,max(df.dev$GAM.age.AdjRsq)))) +
      theme_void()
    
    
    ageEffect_surf <- plot_grid(l_glasser_ageEffect, r_glasser_ageEffect, ncol = 1)
    
  } else if (atlas=="schaefer400x17" | atlas=="schaefer400"){
    r_df.dev <-  df.dev
    r_df.dev$region<- gsub(x = df.dev$region, pattern = "17Networks_RH_", replacement = "")
    r_df.dev <- r_df.dev[-which(grepl("17Networks_LH_", r_df.dev$region)),]
    
    l_df.dev <-  df.dev
    l_df.dev$region<- gsub(x = df.dev$region, pattern = "17Networks_LH_", replacement = "")
    l_df.dev <- l_df.dev[-which(grepl("17Networks_RH_", l_df.dev$region)),]
    
    r_schaefer17_400_ageEffect <- ggplot() + geom_brain(data=r_df.dev, 
                                                        atlas=get(ggseg_atlas), 
                                                        mapping=aes(fill=GAM.age.AdjRsq), 
                                                        show.legend=TRUE, 
                                                        hemi = "right") + 
      scale_fill_gradientn(colors= c("#00A3A7FF", "#FFFFFF","#C75DAAFF"), 
                           limits = range(df.dev$GAM.age.AdjRsq), 
                           values=rescale(c(min(df.dev$GAM.age.AdjRsq),0,max(df.dev$GAM.age.AdjRsq)))) +
      theme_void()
    
    l_schaefer17_400_ageEffect <- ggplot() + geom_brain(data=l_df.dev, 
                                                        atlas=schaefer17_400, 
                                                        mapping=aes(fill=GAM.age.AdjRsq), 
                                                        show.legend=TRUE, 
                                                        hemi = "left") + 
      labs(title = paste(atlas, metric)) +
      scale_fill_gradientn(colors= c("#00A3A7FF", "#FFFFFF","#C75DAAFF"), 
                           limits = range(df.dev$GAM.age.AdjRsq), 
                           values=rescale(c(min(df.dev$GAM.age.AdjRsq),0,max(df.dev$GAM.age.AdjRsq)))) +
      theme_void()
    
    
    ageEffect_surf <- plot_grid(l_schaefer17_400_ageEffect, r_schaefer17_400_ageEffect, ncol = 1)
  } else {
    ageEffect_surf <- ggplot() + geom_brain(data=df.dev, 
                                            atlas=get(ggseg_atlas), 
                                            mapping=aes(fill=GAM.age.AdjRsq), 
                                            show.legend=TRUE, 
                                            position = position_brain(hemi ~ side)) + 
      labs(title = paste(atlas, metric)) +
      scale_fill_gradientn(colors= c("#00A3A7FF", "#FFFFFF","#C75DAAFF"), 
                           limits = range(df.dev$GAM.age.AdjRsq), 
                           values=rescale(c(min(df.dev$GAM.age.AdjRsq),0,max(df.dev$GAM.age.AdjRsq)))) +
      theme_void()
    
  }
  return(ageEffect_surf)
}

## FIGURE: SA axis average rank vs. age effect - before spin test
# @param axis A dataframe with SA-axis rank and GAM results
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name
gam_figure <- function(axis, metric, atlas){
  axis <- axis %>% mutate(SA.axis_rank_signif = ifelse(significant.fdr == 1, SA.axis_rank, NA))
  AgeEffect_figure <- ggplot(axis, aes(x=SA.axis_rank, y=GAM.age.AdjRsq, fill = SA.axis_rank_signif)) + 
    geom_point(color = "white", shape = 21, size=2) +
    scale_fill_gradient2(low="#24bde0",high = "#dc665e", mid = "#F8E97D", midpoint = median(axis$SA.axis_rank), na.value="gray") +
    labs(fill = "SA Axis Rank", x="\nSensorimotor-Association Axis Rank\n", y=expression(paste("Age Effect (Delta Adj", " R"^2, ")"))) +
    ggtitle(paste(metric, "-", atlas)) +
    geom_smooth(data = axis, method='lm', se=TRUE, fill=alpha(c("gray70"),.9), col="black") +
    theme(plot.title = element_text(size=12, hjust=0.5),
          axis.title.x=element_text(size=12, color = "black"),
          axis.title.y=element_text(size=12, color = "black"),
          axis.line = element_line(color = "black"),
          axis.text=element_text(size=12, color = "black"),
          panel.background=element_blank(),
          legend.position = "right") +
    stat_cor(method = "spearman",label.x.npc = 0.6, label.y.npc = 0.8,size = 4, r.digits = 3)
  return(AgeEffect_figure)
}



## FIGURE: SA axis average rank vs. age effect - with spin test, with Val's color palette
# included all datapoints
# @param axis A dataframe with SA-axis rank and GAM results
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name
# @param r_text An expression including info on p_spin and r
# @param x_pos A numeric for horizontal position of r_text
# @param y_pos A numeric for vertical position of r_text

gam_figure_p.spin_Val <- function(axis, metric, atlas, r_text, x_pos, y_pos){
  axis <- axis %>% mutate(SA.axis_rank_signif = ifelse(significant.fdr == 1, SA.axis_rank, NA))
  AgeEffect_figure <- ggplot(axis, aes(x=SA.axis_rank, y=GAM.age.AdjRsq, fill = SA.axis_rank_signif)) + 
    geom_point(color = "gray", shape = 21, size=3.5) + 
    scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", guide = "colourbar", aesthetics = "fill", name = NULL, midpoint = 100) +
    scale_colour_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", guide = "colourbar", aesthetics = "color", name = NULL, midpoint = 100) +
    labs(fill = "SA Axis Rank", x="\nSensorimotor-Association Axis Rank\n", y=expression(paste("Age Effect (Delta Adj", " R"^2, ")"))) +
    #ggtitle(paste(metric, "-", atlas)) + 
    ylim(-0.055, 0.17) +  # ylim(-0.065, 0.15) for GBC; ylim(-0.055, 0.17) for WNC; ylim(-0.07, 0.10) for BNC;  
    geom_smooth(data = axis, method='lm', se=TRUE, fill=alpha(c("gray70"),.9), col="black") +
    theme(
      axis.title.x=element_text(size=18, color = "black", margin = margin(t = -15, r = 0, b = 0, l = 0)),
      axis.title.y=element_text(size=18, color = "black"),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=18, color = "black"),
      panel.background=element_blank(),  
      legend.position = "none") +
    annotate(geom="text", x=x_pos, y=y_pos, label=r_text, color="black", size=7)
  return(AgeEffect_figure)
}



## FIGURE: Developmental Trajectories 
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
devTraj_figure <- function(metric, atlas) {
  gam.smooths <- get(paste0("gam.smooths.", atlas, '.', metric))
  axis <- get(paste0(metric, ".axis_", atlas))  
  smooths_figure <- ggplot(gam.smooths,aes(age,fit,group=index)) + 
    geom_line(data = gam.smooths, size=1, alpha = .6, aes(color=SA.axis_rank)) + 
    paletteer::scale_color_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(axis$SA.axis_rank), max(axis$SA.axis_rank)), oob = squish) + 
    labs(color = "SA Axis Rank",
         x="Age", 
         y="Predicted Global Brain Connectivity") +
    ggtitle(paste(metric, "-", atlas)) + 
    theme(
      axis.title.x=element_text(size=12, color = "black"),
      axis.title.y=element_text(size=12, color = "black"),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=12, color = "black"),
      panel.background=element_blank(),
      legend.position = "bottom") + coord_cartesian(expand = FALSE, xlim = c(8, 23))
  return(smooths_figure)
}


## FIGURE: developmental trajectories, centered on zero - separated by position on SA axis
# @param atlas A character string of atlas name 
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param SA_position A character string of SM, middle axis, or association pole (i.e. "SM", "middle.axis", "assoc.pole")

make_SM.middle.assoc.fig_centered <- function(atlas, metric, SA_position){
  smooth_fits <- get(paste0(metric, "_SM.middle.assoc"))
  smooth_fits <- smooth_fits[[atlas]][[SA_position]]
  if(SA_position == "SM"){
    SM.middle.assoc_centered <- ggplot(smooth_fits,aes(age,est,group=SA.axis_rank)) + 
      geom_line(data = smooth_fits, size=.7, alpha = .5, aes(color=SA.axis_rank)) + 
      ggtitle(paste(metric, "-", atlas, "Sensorimotor Pole")) + 
      scale_color_gradient2(low="#006CA5FF",high = "#A0CEB6FF", mid = "#009BA8FF", midpoint = median(smooth_fits$SA.axis_rank)) +
      theme_classic()
  } else if (SA_position =="middle.axis") {
    SM.middle.assoc_centered <- ggplot(smooth_fits,aes(age,est,group=SA.axis_rank)) + 
      geom_line(data = smooth_fits, size=.7, alpha = .5, aes(color=SA.axis_rank)) + 
      ggtitle(paste(metric, "-", atlas, "Middle Axis")) +   
      scale_color_gradient2(low="#F8E69DFF",high = "#DB4500FF", mid = "#F0BB55FF", midpoint = median(smooth_fits$SA.axis_rank)) + 
      theme_classic()
  } else {
    SM.middle.assoc_centered <- ggplot(smooth_fits,aes(age,est,group=SA.axis_rank)) + 
      geom_line(data = smooth_fits, size=.7, alpha = .5, aes(color=SA.axis_rank)) + 
      ggtitle(paste(metric, "-", atlas, "Association Pole")) +   
      scale_color_gradient2(low="#F0BB55FF",high = "#A51122FF", mid = "#DB4500FF", midpoint = median(smooth_fits$SA.axis_rank)) +
      theme_classic()
  }
  return(SM.middle.assoc_centered)
}



## FIGURE: Temporal Sliding Window -- Magnitude and direction of developmental change in connectivity for each cortical region
# @param metric A character string of connectivity metric (i.e. "GBC")
# @param atlas A character string of atlas name 
# @param dataset A character string of dataset (i.e. "NKI")

tempSliding_fig <- function(metric, atlas) {
  gam.derivatives.atlas <- get(paste0("gam.", metric, ".derivatives.", atlas))
  derivative.colorbar <- c("#009593FF","#269D91FF","#53A592FF","#71AD95FF","#8AB49BFF",
                           "#FEFFFFFF",
                           "#D79D86FF","#D68F7FFF","#D5807AFF","#D46F78FF","#D35C79FF")
  derivative.colorbar <- rev(derivative.colorbar)
  tempSliding_ageDerivatives <- ggplot(data=gam.derivatives.atlas,aes(x = age, y = SA.axis_rank, group = index)) +
    geom_line(aes(color=derivative), size=1) +
    scale_color_gradientn(colours = derivative.colorbar, oob = squish) +
    theme_classic() +
    ylab("Sensorimotor-Association Axis Rank\n") +
    xlab("\nAge") +
    labs(color="") +
    theme(axis.text = element_text(size = 7, family = "Arial", color = c("black")), axis.title = element_text(size = 7, family = "Arial", color = c("black")), axis.line = element_line(size =.22), axis.ticks = element_line(size = .22), legend.text = element_text(size = 6, family = "Arial", color = c("black"))) +
    scale_x_continuous(breaks=c(8, 10, 12, 14, 16, 18, 20, 22), expand = c(0,.45))  
  return(tempSliding_ageDerivatives)
}


## FIGURE: Surface Plot
# @param atlas A character string of atlas name  
make_surface_plot <- function(atlas) {
  double_age_SA.diff <- get(paste0("double_age_SA.diff_", atlas))
  g2 <- gam(GAM.age.AdjRsq ~ te(parcel2.SA.vec,parcel1.SA.vec,k=3),data = double_age_SA.diff)
  a <- quantile(double_age_SA.diff$GAM.age.AdjRsq, prob=c(0.01, 0.5, 0.99))[1]
  b <- quantile(double_age_SA.diff$GAM.age.AdjRsq, prob=c(0.01, 0.5, 0.99))[3]
  
  surface.plot.fig <- gg_tensor(g2)+theme_classic() +
    theme(
      axis.title.x=element_text(size=24, color = "black"),
      axis.title.y=element_text(size=24, color = "black"),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=32, color = "black")) + 
    scale_fill_gradient2(  
      low = "blue",high = "red",mid = "white",midpoint = 0,
      name=expression(paste('Age Effect (',Delta,R^2[adj],')'))) +
    xlab("Parcel SA Rank") + 
    ylab("Parcel SA Rank") + 
    geom_vline(xintercept = mean(range(double_age_SA.diff$parcel1.SA.vec)),
               linetype='dashed',size=1) + 
    geom_hline(yintercept = mean(range(double_age_SA.diff$parcel1.SA.vec)),
               linetype='dashed',size=1) +
    theme(legend.position='none')
  return(surface.plot.fig)
  
}