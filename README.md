# MiCloud
Title: A unified web platform for comprehensive microbiome data analysis

Maintainer: Won Gu won.gu@stonybrook.edu

Depends: R (≥ 3.5)

## Description
MiCloud is a web application for statistical analysis of microbiome data. It provides user-friendly environment with extended options for data processsing, analysis and graphical procedures. MiCloud can conduct ecological (alpha- and beta- diversity) and taxonomic analyses for binary or continuous trait, cross-sectional or longitudinal study, and with or without covariate adjustments.

The application is available at:

URLs
*Web Platform: http://223.194.200.160:3838/app/micloud
*Local Implementation: https://github.com/wg99526/MiCloudGit

## Launch App 
Users can launch the app in their local computer by typing the following command to R terminal:

```
install.packages("shiny")

library(shiny)

runGitHub("MiCloudGit", "wg99526", ref = "main")
```


### Prerequisites

phyloseq
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
```

biomformat
```
source("https://bioconductor.org/biocLite.R")
biocLite("biomformat")
```

GLMM-MiRKAT
```
library(devtools)
install_github("hk1785/GLMM-MiRKAT", force=T)
```

NBZIMM
```
library(remotes)
install_github("nyiuab/NBZIMM", force=T, build_vignettes=F)
```


CRAN Packages

```
list.of.packages <- c('seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable', 'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante', 'entropart', 'lme4', 'lmerTest', 'broom.mixed', 'gee', 'geepack', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn', 'CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr', 'devtools', 'betareg', 'reticulate',   'nlme', 'glmmTMB', 'glmm', 'remotes', 'gridGraphics', 'compositions')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```


