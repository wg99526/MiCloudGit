# MiCloud
Title: A unified web platform for comprehensive microbiome data analysis

Maintainer: Won Gu won.gu@stonybrook.edu

Depends: R (≥ 3.5)

## Description
MiCloud is a unified web platform for comprehensive microbiome data analysis. MiCloud provides step-by-step user-friendly web environments for a breadth of data processing, analytic and graphical procedures. MiCloud performs both ecological (alpha- and beta-diversity) and taxonomical (phylum, class, order, family, genus, species) analyses for various types of host phenotypes (or disease status) and study designs with or without covariate adjustment(s). More details are as follows.

- Interactive procedures for various data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), quality controls (kingdom, library size, mean proportion, taxonomic name) and data transformations (alpha- and beta-diversity, rarefaction, proportion, centered log-ratio).
- Comparative analysis for both ecological (alpha- and beta-diversity) indices and microbial taxa in relative abundance (phylum, class, order, family, genus, species).
- Comparative analysis for both binary and continuous traits of host phenotypes (or medical interventions, disease status or environmental/behavioral factors).
- Comparative analysis with or without covariate (e.g., age, gender) adjustment(s) for either cross-sectional or longitudinal/family-based microbiome study design.
- Adjustable/downloadable/publishable data, tables and graphs.

<br/>

URLs
- Web Platform: http://223.194.200.160:3838/app/micloud
- Local Implementation: https://github.com/wg99526/MiCloudGit

<br/>

Reference: Gu, W., Moon, J., Chisina, C., Kang, B., Park, T., Koh, H. MiCloud: A unified web platform for comprehensive microbiome data analysis. (Under review)

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


