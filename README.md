# MiCloud
Title: A unified web platform for comprehensive microbiome data analysis

Maintainer: Won Gu won.gu@stonybrook.edu

Depends: R (≥ 3.5)

## Description
MiCloud is a web application for statistical analysis of microbiome data. It provides user-friendly environment with extended options for data processsing, analysis and graphical procedures. MiCloud can conduct ecological (alpha- and beta- diversity) and taxonomic analyses for binary or continuous trait, cross-sectional or longitudinal study, and with or without covariate adjustments.

The application is available at:

URL: http://223.194.200.160:3838/

## Launch App
The web application supports up to five concurrent users. If the application is too crowded, you can launch the app from your computer by typing the following command after installing prerequisite packages.

In R terminal:
```
library(shiny)

runGitHub("MiCloudGit", "wg99526", ref = "main")
```

Rstudio is not required to launch MiCloud.

## Data Input
Four components are required to get started: feature table, taxonomic table, metadata, and phylogenetic tree. Users can upload them in phyloseq (McMurdie and Holmes, 2013) or individually. More information on compatible data format for required data components can be found below

### Phyloseq
Phyloseq is a data format that integrates all four components (feature table, taxonomic table, metadata, and phylogenetic tree) in a single R object. Users can upload it using .Rdata and .rds files.

### Individual Data
1. Feature table (.txt, .csv, .biom) should contain counts
2. Taxonomic table (.txt, .csv) should contain taxonomic names
3. Metadata/Sample (.txt, .csv) shuold contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors
4. Phylogenetic tree (.tre, .nwk) should be rooted tree


### Example Data
Two example sets (16s rRNA and Shotgun sequencing data) are available. The datasets will be downloaded as 'biom.Rdata' for phyloseq and 'biom.zip' for individual data. 'biom.zip' file contains feature table (otu.tab.txt), taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre). The example data is UK twin study data (Goodrich et al, 2014) and publicly available in the European Bioinformatics Institute (EMBL-EBI) database. 

## Quality Control


### Prerequisites

shiny
```
install.packages("shiny")
```

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


Following codes includes rest of the packages required to launch the app. You can simply copy and paste these codes into R terminal, then it will examine which packages do not exist in your workspace and install them automatically.

```
list.of.packages <- c('seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable', 'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante',  'entropart', 'lme4', 'lmerTest', 'broom.mixed', 'gee', 'geepack', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn', 'CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr', 'devtools', 'betareg', 'reticulate',   'nlme', 'glmmTMB', 'glmm', 'remotes', 'gridGraphics', 'compositions')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```


