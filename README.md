# MiCloud
Title: A web-based microbiome data analysis tool for binary or continuous trait, 
cross-sectional or longitudinal study, and with or without covariate adjustment

Maintainer: Won Gu won.gu@stonybrook.edu

Depends: R (≥ 3.0.2)

## Description
MiCloud is a web application for statistical analysis of microbiome data. It provides user-friendly environment with extended options for data processsing, analysis and graphical procedures. MiCloud can conduct ecological (alpha- and beta- diversity) and taxonomic analyses for binary or continuous trait, cross-sectional or longitudinal study, and with or without covariate adjustments.

URL: https://github.com/wg99526/micloud

The application is available at:

https://

## Launch app

In R terminal:
```
library(shiny)

runGithub("micloud", "wg99526", ref = "main")
```

Using Rstudio is recommended, but not required

## Prerequisites

shiny
```
install.packages("shiny")
```

```
list.of.packages <- c('seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', <br> 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable', 'DT', <br> 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', <br> 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante', <br> 'entropart', 'lme4', 'lmerTest', 'broom.mixed', 'gee', <br> 'geepack', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn', <br> 'CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', <br> 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', <br> 'DiagrammeR', 'stringr', 'devtools', 'betareg', 'reticulate', <br> 'nlme', 'glmmTMB', 'glmm', 'remotes', 'gridGraphics', 'compositions')
```

```
if(!require('phyloseq')) remotes::install_github('joey711/phyloseq')
if(!require('biomformat')) remotes::install_github('joey711/biomformat')
if(!require('GLMM-MiRKAT')) remotes::install_github('hk1785/GLMM-MiRKAT')
if(!require('NBZIMM')) remotes::install_github('nyiuab/NBZIMM')
```
