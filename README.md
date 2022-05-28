# MiCloud

Title: MiCloud: A Unified Web Platform for Comprehensive Microbiome Data Analysis

Version: 1.21

Date: 2022-05-28

Maintainer: Won Gu <wpg5129@psu.edu>

Description: MiCloud is a unified web platform for comprehensive microbiome data analysis. MiCloud provides step-by-step user-friendly web environments for a breadth of data processing, analytic and graphical procedures. MiCloud performs both ecological (alpha- and beta-diversity) and taxonomical (phylum, class, order, family, genus, species) analyses for various types of host phenotypes (or disease status) and study designs with or without covariate adjustment(s). More details are as follows.

* Interactive procedures for various data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), quality controls (kingdom, library size, mean proportion, taxonomic name) and data transformations (alpha- and beta-diversity, rarefaction, proportion, centered log-ratio).
* Comparative analysis for both ecological (alpha- and beta-diversity) indices and microbial taxa in relative abundance (phylum, class, order, family, genus, species).
* Comparative analysis for both binary and continuous traits of host phenotypes (or medical interventions, disease status or environmental/behavioral factors).
* Comparative analysis with or without covariate (e.g., age, gender) adjustment(s) for either cross-sectional or longitudinal/family-based microbiome study design.
* Adjustable/downloadable/publishable data, tables and graphs.

Depends: R(>= 3.5), 'phyloseq', 'biomformat', 'GLMM-MiRKAT', 'NBZIMM', 'seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable', 'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante', 'entropart', 'lme4', 'lmerTest', 'broom.mixed', 'gee', 'geepack', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn', 'CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr', 'devtools', 'betareg', 'reticulate',   'nlme', 'glmmTMB', 'glmm', 'remotes', 'gridGraphics', 'compositions'

## URLs

* Web application (online implementation): http://micloud.kr   
* GitHub repository (local implementation): https://github.com/wg99526/MiCloudGit
 
## References

* Gu, W., Moon, J., Chisina, C., Kang, B., Park, T., Koh, H. MiCloud: A unified web platform for comprehensive microbiome data analysis. (*_Under revision_*)

# Prerequites

shiny
```
install.packages("shiny")
```

# Launch App

```
library(shiny)

runGitHub("MiCloudGit", "wg99526", ref = "main")
```

# Troubleshooting Tips

If you have any problems for using MiCloud, please report in Issues (https://github.com/wg99526/MiCloudGit/issues) or email Ms. Won Gu (won.gu@stonybrook.edu).


# External Resources

MiCloud does not take raw sequence data. For the raw sequence data processing and microbiome profiling, we recommend other popular and well-established bioinformatic pipelines, such as 
* Nephele (https://nephele.niaid.nih.gov), Qiita (https://qiita.ucsd.edu), QIIME2 (q2studio) (https://qiime2.org) and PUMAA (https://sites.google.com/g.ucla.edu/pumaa) for web platforms, or
* QIIME (http://qiime.org), QIIME2 (q2cli) (https://qiime2.org), MG-RAST (https://www.mg-rast.org), Mothur (https://mothur.org), MEGAN (http://ab.inf.uni-tuebingen.de/software/megan6) and MetaPhlAn (https://huttenhower.sph.harvard.edu/metaphlan) for command line interfaces.

