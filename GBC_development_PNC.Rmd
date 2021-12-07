---
title: "GBC Developmental Trajectories: PNC"
author: "Valerie Jill Sydnor"
output:
  html_document:
    code_folding: show
    highlight: haddock
    theme: lumen
    toc: yes
    toc_depth: 4
    toc_float: yes
  pdf_document:
    toc: yes
    toc_depth: '4'
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
require(ggplot2)
library(cowplot)
library(cifti)
library(ggseg)
library(ggsegExtra)
library(ggsegGlasser)
library(ggcorrplot)
library(viridis)
library(dplyr)
library(stringr)
library(tidyr)
```

```{r warning=FALSE, include=FALSE}
## Sensorimotor-Association Axis

SAaxis.glasser <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/Glasser360_MMP_WholeBrain/Sensorimotor_Association_Axis_AverageRanks.csv")
glasser.parcel.labels <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/glasser360_regionlist.csv", header = T) 
```

```{r include=FALSE}
## Brain Organization

brains.glasser <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/S-A_ArchetypalAxis/Glasser360_MMP_WholeBrain/brainmaps_glasser.csv")
brains.glasser$label <- glasser.parcel.labels$label
```

```{r include=FALSE}
## GAM Results

gam.gbc.age.glasser <- read.csv("GAMresults.GBC.age.glasser.csv")
gam.gbc.age.glasser <- gam.gbc.age.glasser %>% select(-label)
gam.gbc.age.glasser$label <- glasser.parcel.labels$label

gam.smooths.glasser <- read.csv("GAMsmoothfits.GBC.age.glasser.csv")
gam.smooths.glasser$orig_parcelname <- gam.smooths.glasser$label
gam.smooths.glasser <- gam.smooths.glasser %>% select(-label)

gam.agepeaks.glasser <- read.csv("GAMpeaks.GBC.age.glasser.csv")
gam.agepeaks.glasser <- gam.agepeaks.glasser %>% select(-label)
gam.agepeaks.glasser$label <- glasser.parcel.labels$label
```


```{r include=FALSE}
## SNR Parcel Mask

SNR.mask <- read.csv("/cbica/projects/spatiotemp_dev_plasticity/Maps/parcellations/surface/SNRmask_glasser360.csv")
```


```{r include=FALSE}
# Combine into Final Dfs

df.list <- list(SAaxis.glasser,interneurons,brains.glasser,gam.gbc.age.glasser,SNR.mask)
gbc.axis.all <- Reduce(function(x,y) merge(x,y, all=TRUE, sort=F), df.list) 
gbc.axis <- gbc.axis.all %>% filter(SNR.mask != 0)
```
```{r include=FALSE}
SAaxis <- merge(glasser.parcel.labels, SAaxis.glasser, by="label", sort=F)
gam.smooths.glasser <- left_join(gam.smooths.glasser, SAaxis, "orig_parcelname")
gam.smooths.glasser <- left_join(gam.smooths.glasser, SNR.mask, by="label")
gam.smooths.glasser <- gam.smooths.glasser %>% filter(SNR.mask != 0)
```
```{r include=FALSE}
gam.agepeaks.glasser <- left_join(gam.agepeaks.glasser, SAaxis.glasser, by="label")
gam.agepeaks.glasser <- left_join(gam.agepeaks.glasser, SNR.mask, by="label")
gam.agepeaks.glasser <- gam.agepeaks.glasser %>% filter(SNR.mask != 0)
```


# GBC Developmental Effects

## Number of significant parcels (FDR corrected)

```{r echo=FALSE}
gbc.axis$Anova.age.pvalue.fdr <- p.adjust(gbc.axis$Anova.age.pvalue, method=c("fdr"))
cat(sprintf("There are %s/%s significant parcels", sum(gbc.axis$Anova.age.pvalue.fdr < 0.05), nrow(gbc.axis)))

gbc.axis$significant.fdr <- gbc.axis$Anova.age.pvalue.fdr < 0.05
gbc.axis$significant.fdr[gbc.axis$significant.fdr == TRUE] <- 1
gbc.axis$significant.fdr[gbc.axis$significant.fdr == FALSE] <- 0
```

```{r, echo=F, include=F}
ggseg(.data = gbc.axis, atlas = "glasser", mapping=aes(fill=as.factor(significant.fdr)), position = c("stacked")) + theme_void() + scale_fill_discrete(name = "significant.fdr", labels = c("Not Sig","Sig"))
```

## Age Effect by S-A Axis position

```{r}
cor.test(gbc.axis$GAM.age.AdjRsq, gbc.axis$finalrank.wholebrain, method=c("spearman"), exact=F)
```

```{r, echo=F}
ggplot(gbc.axis, aes(x=finalrank.wholebrain, y=GAM.age.AdjRsq, fill = finalrank.wholebrain)) + 
geom_point(color = "white",shape=21, size=2.5) +
scale_fill_gradient2(low="goldenrod1",high = "#6f1282", midpoint = median(gbc.axis$finalrank.wholebrain)) +
labs(x="\nSensorimotor-Association Axis Rank\n", y="\nAge Effect (Delta Adj R^2)\n") +
geom_smooth(method='lm', se=TRUE, fill=alpha(c("gray70"),.7), col="black") +
theme(
axis.title.x=element_text(size=20, color = "black"),
axis.title.y=element_text(size=20, color = "black"),
axis.line = element_line(color = "black"),
axis.text=element_text(size=15, color = "black"),
panel.background=element_blank(),
legend.position = "none")
```

## GBC Developmental Trajectories

```{r}
ggplot(gam.smooths.glasser,aes(age,fit,group=index)) + geom_line(data = gam.smooths.glasser, size=.75, alpha = .3, aes(color=finalrank.wholebrain)) + scale_color_gradient2(low="goldenrod1",high = "#6f1282", midpoint = median(gam.smooths.glasser$finalrank.wholebrain)) + theme(legend.position = "none")
```

## Age of Maximal GBC 

```{r echo=FALSE}
ggseg(.data = gam.agepeaks.glasser, atlas = "glasser", mapping=aes(fill=age.peak), position = c("stacked")) + theme_void() + scale_fill_gradientn(colors = viridis_pal(option="B")(10))
```
