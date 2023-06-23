rm(list = ls())

list.of.packages <- c('seqinr', 'shinydashboard', 'dashboardthemes', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable', 
                      'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'zCompositions', 'dplyr', 'forestplot', 'quantreg', 'fossil', 'picante',
                      'entropart', 'lme4', 'lmerTest', 'broom.mixed', 'gee', 'geepack', 'dirmult', 'robustbase', 'robCompositions', 'BiasedUrn',
                      'CompQuadForm', 'GUniFrac', 'ecodist', 'MiRKAT', 'gridExtra', 'ggplot2', 'patchwork', 'ggthemes', 'erer', 'DiagrammeR', 'stringr',
                      'devtools', 'betareg', 'reticulate', 'nlme', 'glmmTMB', 'glmm', 'remotes', 'gridGraphics', 'compositions', 'kableExtra', 'modelr')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!require('phyloseq')) remotes::install_github('joey711/phyloseq')
if(!require('biomformat')) remotes::install_github('joey711/biomformat')
if(!require('GLMM-MiRKAT')) remotes::install_github('hk1785/GLMM-MiRKAT', force = T)
if(!require('NBZIMM')) remotes::install_github('nyiuab/NBZIMM')


library(seqinr)
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(tidyverse)
library(phyloseq)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(biomformat)
library(phangorn)
library(bios2mds)
library(kableExtra)
library(modelr)
library(zip)

source("Source/MiDataProc.Data.Upload.R")
source("Source/MiDataProc.Alpha.Cross.Sectional.R")
source("Source/MiDataProc.Alpha.Longitudinal.R")
source("Source/MiDataProc.Beta.Cross.Sectional.R")
source("Source/MiDataProc.Beta.Longitudinal.R")
source("Source/MiDataProc.Taxa.Cross.Sectional.R")
source("Source/MiDataProc.Taxa.Longitudinal.R")

# COMMENTS
{
  TITLE = p("MiCloud: A Unified Web Platform for Comprehensive Microbiome Data Analysis", style = "font-size:18pt")
  HOME_COMMENT = p(strong("MiCloud", style = "font-size:15pt"), "is a unified web platform for comprehensive microbiome data analysis. 
                   MiCloud provides step-by-step user-friendly web environments for a breadth of data processing, analytic and graphical procedures. 
                   MiCloud performs both ecological (alpha- and beta-diversity) and taxonomical (phylum, class, order, family, genus, 
                   species) analyses for various types of host phenotypes (or disease status) and study designs with or without covariate 
                   adjustment(s). More details are as follows.", style = "font-size:13pt")
  HOME_COMMENT1 = ("Interactive procedures for various data inputs (.rdata, .rds, .biom, .txt, .csv, .tsv, .tre), quality controls (kingdom, 
                   library size, mean proportion, taxonomic name) and data transformations (alpha- and beta-diversity, rarefaction, 
                   proportion, centered log-ratio).")
  HOME_COMMENT2 = ("Comparative analysis for both ecological (alpha- and beta-diversity) indices and microbial taxa in relative abundance 
                   (phylum, class, order, family, genus, species).")
  HOME_COMMENT3 = ("Comparative analysis for both binary and continuous traits of host phenotypes 
                   (or medical interventions, disease status or environmental/behavioral factors).")
  HOME_COMMENT4 = ("Comparative analysis with or without covariate (e.g., age, gender) adjustment(s) for either 
                   cross-sectional or longitudinal/family-based microbiome study design.")
  HOME_COMMENT5 = ("Adjustable/downloadable/publishable data, tables and graphs.")
  HOME_COMMENT6 = p(strong("URLs:", style = "font-size:15pt"), a("Web Server (http://micloud.kr)", href = "http://micloud.kr", target = "blank"), 
                                                                 a("Github (https://github.com/wg99526/MiCloudGit)", href = "https://github.com/wg99526/MiCloudGit", target = "blank"), style = "font-size:13pt")
  HOME_COMMENT7 = p(strong("Maintainer:", style = "font-size:13pt"), "Won Gu (wpg5129@psu.edu)", style = "font-size:13pt")
  HOME_COMMENT8 = p(strong("Reference:", style = "font-size:13pt"), 
                    " Gu W, Moon J, Chisina C, Kang B, Park T, Koh H (2022) MiCloud: A unified web platform for comprehensive microbiome data analysis.",  
                    "PLoS ONE 17(8): e0272354. https://doi.org/10.1371/journal.pone.0272354", style = "font-size:13pt")
  
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), "This should be an '.Rdata' or '.rds' file, and the data should be in 'phyloseq' format (see ", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain all the four necessary data, feature (OTU or ASV) table, taxonomic table, 
                              metadata/sample information, and phylogenetic tree.", 
                              br(), br(), "Details:", br(), br(), 
                              "1) The feature table should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                              (row names are feature IDs and column names are subject IDs).", br(),
                              "2) The taxonomic table should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                              (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                              'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(),
                              "3) The metadata/sample information should contain variables for the subjects about host phenotypes, medical interventions, 
                              disease status or environmental/behavioral factors, where rows are subjects and columns are variables 
                              (row names are subject IDs, and column names are variable names).", br(), 
                              "4) The phylogenetic tree should be a rooted tree. Otherwise, MiCloud automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiCloud will analyze only the matched features and subjects."
                              , style = "font-size:11pt")
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in 'phyloseq' format. The name of the 
                              phyloseq object should be 'biom'. For more details about 'phyloseq', see ", 
                              a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(), 
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'biom.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)", br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > tree <- phy_tree(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), br(), 
                              "You can check if the features are matched and identical across feature table, taxonomic table and 
                              phylogenetic tree, and the subjects are matched and identical between feature table and metadata/sample information 
                              using following code.", br(), br(), 
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                              " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The feature table (.txt or .csv) should contain counts, where rows are features (OTUs or ASVs) and columns are subjects 
                                    (row names are feature IDs and column names are subject IDs). Alternatively, you can upload .biom file processed by QIIME.", br(), 
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                                    (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                                    'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). Alternatively, you can upload .tsv file processed by QIIME.", br(), 
                                    "3) The metadata/sample information (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, 
                                    disease status or environmental/behavioral factors, where rows are subjects and columns are variables (row names are subject IDs, and 
                                    column names are variable names).", br(), 
                                    "4) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiCloud automatically roots the tree through midpoint 
                                    rooting (phangorn::midpoint). The tip labels of the phylogenetic tree are feature IDs.", br(), br(), 
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. 
                                    MiCloud will analyze only the matched features and subjects.", style = "font-size:11pt")
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'biom.zip'. This zip file contains four necessary data, feature table (otu.tab.txt), 
                                     taxonomic table (tax.tab.txt), metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), br(),
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE) ", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     " > sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE) ", br(),
                                     "> tree <- read.tree(file = 'tree.tre')", br(), br(), 
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                                     and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt")
  
  EXTERNAL_RESOURCE_COMMENT = p("MiCloud does not take raw sequence data. For the raw sequence data processing and microbiome profiling, we recommend following popular and well-established bioinformatic pipelines.", br(), br(),
                              "For web platforms:", p(" ", style = "margin-bottom: 10px;"),"Nephele (https://nephele.niaid.nih.gov), Qiita (https://qiita.ucsd.edu), QIIME2 (q2studio) (https://qiime2.org) and PUMAA (https://sites.google.com/g.ucla.edu/pumaa)", br(), br(),
                              "For command line interfaces:", p(" ", style = "margin-bottom: 10px;"), "QIIME (http://qiime.org), QIIME2 (q2cli) (https://qiime2.org), MG-RAST (https://www.mg-rast.org), Mothur (https://mothur.org), MEGAN (http://ab.inf.uni-tuebingen.de/software/megan6) and MetaPhlAn (https://huttenhower.sph.harvard.edu/metaphlan)", style = "font-size:11pt")
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. Alternatively, you can type 'Fungi' for ITS data 
                         or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove subjects that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per subject.", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.002%.",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p('Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. Default is "", "metagenome", "gut metagenome", "mouse gut metagenome".',
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p('Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (i.e., taxonomic names that contain 
                            the specified character strings). Multiple character strings should be separated by a comma. Default is "uncultured", "incertae", "Incertae",
                            "unidentified", "unclassified", "unknown".',style = "font-size:11pt")
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative β-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  DATA_TRANSFORM_COMMENT = p("Transform the data into four different formats 1) count, 2) count (rarefied), 3) proportion, 
                             4) CLR (centered log ratio) (Aitchison, 1982) for each taxonomic rank (phylum, class, order, familiy, genus, species).")
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44:2:139-77")
}

# UI
{
  ui = dashboardPage(
    title = "MiCloud",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
    dashboardSidebar(
      tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
      sidebarMenu(id = "side_menu",
                  menuItem("Home", tabName = "home", icon = icon("home")),
                  menuItem("Data Processing",  icon = icon("file-text-o"),
                           menuSubItem("Data Input", tabName = "step1", icon = icon("mouse")),
                           menuSubItem("Quality Control", tabName = "step2", icon = icon("chart-bar"))),
                  menuItem("Ecological Analysis",  icon = icon("chart-pie"),
                           menuSubItem("Diversity Calculation", tabName = "divCalculation", icon = icon("calculator")),
                           menuSubItem("Alpha Diversity", tabName = "alphaDivanalysis", icon = icon("font")),
                           menuSubItem("Beta Diversity", tabName = "betaDivanalysis", icon = icon("bold"))),
                  menuItem("Taxonomical Analysis",  icon = icon("disease"),
                           menuSubItem("Data Transformation", tabName = "dataTransform", icon = icon("th-large")),
                           menuSubItem("Comparison/Association", tabName = "taxaAnalysis", icon = icon("align-left"))))),
    dashboardBody(
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      shinyDashboardThemes(theme = "onenote"),
      uiOutput("themes"),
      tabItems(
        
        ##### HOME ####
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT,
                    tags$ol(
                      tags$li(HOME_COMMENT1), tags$li(HOME_COMMENT2), tags$li(HOME_COMMENT3), tags$li(HOME_COMMENT4), tags$li(HOME_COMMENT5),
                      style = "font-size:13pt"),
                    HOME_COMMENT6, br(),p(" ", style = "margin-bottom: -20px;"), HOME_COMMENT7, br(),p(" ", style = "margin-bottom: -20px;"), HOME_COMMENT8)),
        
        ##### DATA INPUT ####
        tabItem(tabName = "step1", br(),
                column(width = 6, style='padding-left:0px',
                       box(
                         width = NULL, status = "primary", solidHeader = TRUE,
                         title = strong("Data Input", style = "color:black"),
                         selectInput("inputOption", h4(strong("Data Type?")), c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                         div(id = "optionsInfo", tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), style = "margin-top: -15px"),
                         uiOutput("moreOptions"))),
                column(width = 6, style='padding-left:0px', uiOutput("addDownloadinfo"))),
        
        ##### QC ####
        tabItem(tabName = "step2", br(), 
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               textInput("kingdom", h4(strong("Bacteria?")), value = "Bacteria"),
                               QC_KINGDOM_COMMENT,
                               tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                               tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'), 
                               
                               sliderInput("slider1", h4(strong("Library Size?")), min=0, max=10000, value = 3000, step = 1000),
                               QC_LIBRARY_SIZE_COMMENT1,
                               QC_LIBRARY_SIZE_COMMENT2,
                               
                               sliderInput("slider2", h4(strong("Mean Proportion?")), min = 0, max = 0.1, value = 0.002, step = 0.001,  post  = " %"),
                               QC_MEAN_PROP_COMMENT1,
                               QC_MEAN_PROP_COMMENT2,
                               
                               br(),
                               p(" ", style = "margin-bottom: -20px;"),
                               
                               h4(strong("Erroneous Taxonomic Names?")),
                               textInput("rem.str", label = "Complete Match", value = ""),
                               QC_TAXA_NAME_COMMENT1,
                               
                               textInput("part.rem.str", label = "Partial Match", value = ""),
                               QC_TAXA_NAME_COMMENT2,
                               
                               actionButton("run", (strong("Run!")), class = "btn-info"), br(), br(),
                               uiOutput("moreControls")),
                  mainPanel(width = 9,
                            fluidRow(width = 12,
                                     status = "primary", solidHeader = TRUE, 
                                     valueBoxOutput("sample_Size", width = 3),
                                     valueBoxOutput("OTUs_Size", width = 3),
                                     valueBoxOutput("phyla", width = 3),
                                     valueBoxOutput("classes", width = 3)),
                            fluidRow(width = 12, 
                                     status = "primary", solidHeader = TRUE,
                                     valueBoxOutput("orders", width = 3),
                                     valueBoxOutput("families", width = 3),
                                     valueBoxOutput("genera", width = 3),
                                     valueBoxOutput("species", width = 3)),
                            fluidRow(style = "position:relative",
                                     tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                            tabPanel("Histogram",
                                                     plotlyOutput("hist"),
                                                     sliderInput("binwidth", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                     chooseSliderSkin("Round", color = "#112446")),
                                            tabPanel("Box Plot", 
                                                     plotlyOutput("boxplot"))),
                                     tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                            tabPanel("Histogram",
                                                     plotlyOutput("hist2"),
                                                     sliderInput("binwidth2", "# of Bins:",min = 0, max = 100, value = 50, width = "100%"),
                                                     chooseSliderSkin("Round", color = "#112446")),
                                            tabPanel("Box Plot",
                                                     plotlyOutput("boxplot2"))))))),
        
        ##### DIVERSITY Calculation ####
        tabItem(tabName = "divCalculation", br(),
                column(width = 6, style = 'padding-left:0px',
                       box(title = strong("Diversity Calculation", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           ALPHA_COMMENT, 
                           BETA_COMMENT, 
                           actionButton("divCalcRun", (strong("Run!")), class = "btn-info")),
                       uiOutput("divCalcDownload")),
                column(width = 6, style='padding-left:0px',
                       box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           p("Alpha Diversity", style = "font-size:12pt"),
                           ALPHA_REFERENCES,
                           p("Beta Diversity", style = "font-size:12pt"),
                           BETA_REFERENCES))),
        
        ##### ALPHA DIVERSITY ####
        tabItem(tabName = "alphaDivanalysis", br(),
                fluidRow(
                  tabBox(width = 12,
                         tabPanel(
                           title ="Cross-Sectional",
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("primvars"),
                                          uiOutput("prim_vars_types"),
                                          uiOutput("covariates"), br(), 
                                          uiOutput("alpha_downloadTable"),
                                          uiOutput("alpha_references")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("alpha_display_results"))))),
                         tabPanel(
                           title ="Longitudinal", 
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3, 
                                          uiOutput("primvars_long"),
                                          uiOutput("prim_vars_types_long"),
                                          uiOutput("covariates_long"), br(), 
                                          uiOutput("alpha_downloadTablelong"),
                                          uiOutput("alpha_references_long")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("alpha_display_resultslong")))))))),
        
        ##### BETA DIVERSITY ####
        tabItem(tabName = "betaDivanalysis", br(),
                fluidRow(
                  tabBox(width = 12,
                         tabPanel(
                           title ="Cross-Sectional",
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("beta_primvar_cross"),
                                          uiOutput("beta_prim_vars_types_cross"),
                                          uiOutput("beta_covariates_cross"), br(), 
                                          uiOutput("beta_downloadTable"),
                                          uiOutput("beta_references")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("beta_display_results_cross"))))),
                         tabPanel(
                           title ="Longitudinal",
                           sidebarLayout(
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("beta_primvars_long"),
                                          uiOutput("beta_prim_vars_types_long"),
                                          uiOutput("beta_covariates_long"), br(), 
                                          uiOutput("beta_downloadTablelong"),
                                          uiOutput("beta_references_long")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                uiOutput("beta_display_resultslong")))))))),
        
        ##### Data Transformation ####
        tabItem(tabName = "dataTransform", br(),
                column(width = 6, style='padding-left:0px',
                       box(title = strong("Data Transformation", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           DATA_TRANSFORM_COMMENT,
                           actionButton("datTransRun", (strong("Run!")), class = "btn-info") ),
                       uiOutput("datTransDownload")),
                column(width = 6, style='padding-left:0px', 
                       box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           DATA_TRANSFORM_REFERENCE))),
        
        ##### Taxa Analysis ####
        tabItem(tabName = "taxaAnalysis", br(),
                fluidRow(
                  tabBox(width = 12,
                         tabPanel(
                           title = "Cross-Sectional",
                           sidebarLayout( 
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("primvars_taxa"),
                                          uiOutput("morePrimvar_opt_taxa"),
                                          uiOutput("covariates_taxa"), br(),
                                          uiOutput("downloadTable_taxa"),
                                          uiOutput("taxa_references")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12, 
                                                div(style='height:800px;overflow-y: scroll;', uiOutput("taxa_display")), br(),br(),
                                                uiOutput("Ctaxa_display_dend"))))),
                         tabPanel(
                           title = "Longitudinal",
                           sidebarLayout( 
                             position = "left",
                             sidebarPanel(width = 3,
                                          uiOutput("primvars_taxa.long"),
                                          uiOutput("morePrimvar_opt_taxa.long"),
                                          uiOutput("covariates_taxa.long"), br(),
                                          uiOutput("downloadTable_taxalong"),
                                          uiOutput("taxa_references_long")),
                             mainPanel(width = 9,
                                       fluidRow(width = 12,
                                                div(style='height:800px;overflow-y: scroll;', uiOutput("taxa_displaylong")), br(),br(),
                                                uiOutput("Ltaxa_display_dend"))))))))
      )
    )
  )
}

# Server
server = function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  ## load example data ####
  sub.biom <- readRDS("Data/val_physeq.rds")
  biom <- sub.biom
  
  env <- new.env()
  nm <- load(file = "Data/BMI/biom.MZ.BMI.Rdata", env)[1]
  BMI <- env[[nm]]
  
  ori.biom <- biom
  otu.tab <- otu_table(ori.biom)
  tax.tab <- tax_table(ori.biom)
  tree <- phy_tree(ori.biom)
  sam.dat <- sample_data(ori.biom)
  
  B.otu.tab <- otu_table(BMI)
  B.tax.tab <- tax_table(BMI)
  B.tree <- phy_tree(BMI)
  B.sam.dat <- sample_data(BMI)
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("biom",".Rdata", sep = "")
    },
    content = function(file1) {
      save(biom, file = file1)
    })
  output$downloadData.BMI <- downloadHandler(
    filename = function() {
      paste("biom.MZ.BMI.Rdata", sep = "")
    },
    content = function(file1) {
      save(BMI, file = file1)
    })
  
  output$downloadZip <- downloadHandler(
    filename = function() {
      paste("biom",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt" ,"tree.tre")
      write.table(otu.tab, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(sam.dat, "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  output$downloadZip.BMI <- downloadHandler(
    filename = function() {
      paste("BMI",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("BMI.otu.tab.txt", "BMI.tax.tab.txt", "BMI.sam.dat.txt" ,"BMI.tree.tre")
      write.table(B.otu.tab, "BMI.otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(B.tax.tab, "BMI.tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(B.sam.dat, "BMI.sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(B.tree, "BMI.tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  ## variable define ####
  infile = reactiveValues(biom = NULL, qc_biom = NULL, rare_biom = NULL, qc_biomNA = NULL, rare_biomNA = NULL)
  ds.Ks <- reactiveValues(res = NULL)
  chooseData = reactiveValues(sam.dat = NULL, mon.sin.rev.bin.con = NULL, prim_vars = NULL, alpha.div = NULL, NAadded = NULL,
                              alpha.div.rare = NULL, alpha.div.qc = NULL, taxa.out = NULL, tax.tab = NULL, taxa.outNA = NULL, tax.tabNA = NULL)
  is.results = reactiveValues(result = NULL)
  is.results.long = reactiveValues(result = NULL)
  multi.test = reactiveValues(boolval = FALSE)
  multi.test.long = reactiveValues(boolval = FALSE)
  alpha.categos <- reactiveValues(cat1 = NULL, cat2 = NULL)
  alpha.categos.long <- reactiveValues(cat1 = NULL, cat2 = NULL)
  alpha.data.results = reactiveValues(table.out = NULL, data.q.out = NULL, table.p.out = NULL)
  alpha.results = reactiveValues(bin.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.results.cont = reactiveValues(alpha.con.out = NULL, alpha.table.out = NULL)
  alpha.reg.results = reactiveValues(bin.var = NULL, cov.var = NULL, alpha.div = NULL)
  alpha.resultslong = reactiveValues(bin.var = NULL, id.var = NULL, alpha_div = NULL, alpha.bin.sum.out = NULL)
  alpha.noncovs_res = reactiveValues(con.var = NULL, id.var = NULL, alpha.div = NULL)
  data.alphaBin_res = reactiveValues(table.output = NULL, data.output = NULL, alpha.bin.sum.out = NULL, table.p_outbin = NULL)
  data.results.cont_long = reactiveValues(table.out = NULL, data.q.out = NULL, table_p.out = NULL, alpha.table.out = NULL)
  
  beta.data.results = reactiveValues(data.q.out = NULL)
  beta.results = reactiveValues(result = NULL)
  beta.resultscont = reactiveValues(beta.cont.out = NULL)
  beta.data.results_long = reactiveValues(beta.bin.out = NULL)
  beta.resultscon_long = reactiveValues(beta.con.out = NULL)
  beta.categos <- reactiveValues(cat1 = NULL, cat2 = NULL)
  beta.categos.long <- reactiveValues(cat1 = NULL, cat2 = NULL)
  beta.down.results <- reactiveValues(CS = NULL, LONG = NULL)
  
  taxa.categos <- reactiveValues(cat1 = NULL, cat2 = NULL)
  taxa.data.results = reactiveValues(data.q.out = NULL)
  taxa.results = reactiveValues(bin.var = NULL, cov.var = NULL, id.var = NULL, taxa = NULL, taxa.bin.sum.out = NULL, con.var = NULL, taxa.con.sum.out = NULL, lib.size = NULL)
  taxa.types = reactiveValues(dataType = NULL, regression = NULL)
  taxa.outputs = reactiveValues(DAoutput = NULL, DAoutput_or = NULL, DAoutputlong = NULL)
  
  rcol = reactiveValues(selected = "lightblue") 
  
  ## options to input ####
  observeEvent(input$inputOption,{
    observe({
      if (input$inputOption == "Phyloseq") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("phyloseqData", strong("Please upload your 'phyloseq' data (.Rdata, .rds)", style = "color:black"), 
                      accept = c(".Rdata", ".rds"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Phyloseq_Data', 'Upload', class = "btn-info"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "phyloseqUpload_error",
                         shiny::tags$p("Please upload a Rdata file!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_PHYLOSEQ_COMMENT1
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                downloadButton("downloadData.BMI", "16S", width = '30%', style = "color:black; background-color: red2"),
                downloadButton("downloadData", "Shotgun", width = '30%', style = "color:black; background-color: red2"),
                br(),br(),
                INPUT_PHYLOSEQ_COMMENT2
            ),
            box(title = strong("External Resource", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                EXTERNAL_RESOURCE_COMMENT)
          )
        })
      } else if (input$inputOption == "Individual Data") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("otuTable", strong("Please upload your feature (OTU or ASV) table (.txt, .csv, .biom)", style = "color:black"), 
                      accept = c(".txt", ".csv", ".biom"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("taxTable", strong("Please upload your taxonomic table (.txt, .tsv)", style = "color:black"), 
                      accept = c(".txt", ".tsv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("samData", strong("Please upload your metadata/sample information (.txt, .csv)", style = "color:black"), 
                      accept = c(".txt", ".csv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("tree", strong("Please upload your phylogenetic tree (.tre, .nwk)", style = "color:black"), 
                      accept = c(".tre", ".nwk"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Individual_Data', 'Upload', class = "btn-info"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "textfilesUpload_error",
                         shiny::tags$p("Please upload txt and tre files!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_INDIVIDUAL_DATA_COMMENT
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                downloadButton("downloadZip.BMI", "16S", width = '30%', style = "color:black; background-color: red2"),
                downloadButton("downloadZip", "Shotgun", width = '30%', style = "color:black; background-color: red2"),
                br(),br(),
                INPUT_INDIVIDUAL_DATA_COMMENT2
            ),
            box(title = strong("External Resource", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                EXTERNAL_RESOURCE_COMMENT)
          )
        })
      }
    })
    
  }, ignoreInit = TRUE, once = TRUE, ignoreNULL = TRUE)
  
  observe({
    toggleState("Load_Phyloseq_Data", !is.null(input$phyloseqData))
    toggleState("Load_Individual_Data", 
                !(is.null(input$otuTable) | is.null(input$taxTable) | is.null(input$samData) | is.null(input$tree)))
    toggleState("run", !is.null(infile$biom))
    toggleState("skip", !is.null(infile$biom))
    toggleState("slider1", !is.null(infile$biom))
    toggleState("slider2", !is.null(infile$biom))
    toggleState("kingdom", !is.null(infile$biom))
    toggleState("binwidth", !is.null(infile$biom))
    toggleState("binwidth2", !is.null(infile$biom))
    
    toggleState("divCalcRun", !is.null(infile$rare_biom))
    toggleState("datTransRun", !is.null(infile$rare_biom))
  })
  
  observeEvent(input$Load_Phyloseq_Data, {
    
    if (!is.null(input$phyloseqData)) {
      dataInfile  = reactive({
        phyloseq.data = input$phyloseqData
        ext <- tools::file_ext(phyloseq.data$datapath)
        
        req(phyloseq.data)
        if (ext == "Rdata") {
          phyloseq.dataPath = phyloseq.data$datapath
          e = new.env()
          name <- load(phyloseq.dataPath, envir = e)
          data <- e[[name]]
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else if (ext == "rds") {
          phyloseq.dataPath = phyloseq.data$datapath
          data <- readRDS(phyloseq.dataPath)
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } else {
          shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade")
          shinyjs::delay(5000, shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade"))
          return(NULL)
        }
      })
    } else {
      return(NULL)
    }
    
    if (is.null(dataInfile)) {
      infile$biom <- NULL
      infile$qc_biom <- NULL
      infile$rare_biom = NULL
    } else {
      infile$biom <- dataInfile()
      infile$qc_biom <- dataInfile()
      infile$rare_biom = NULL
    }
    
    updateTabsetPanel(session, "side_menu",
                      selected = "step2")
    rcol$selected = "lightblue"
    
    if (!is.null(infile$biom)) QC$resume()
  })
  observeEvent(input$Load_Individual_Data, {
    shinyjs::disable("Load_Individual_Data")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(3/10, message = "File Check")
        if (!is.null(input$otuTable) & !is.null(input$taxTable) & !is.null(input$samData) & !is.null(input$tree)) {
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            tax.table = input$taxTable
            ext2 <- tools::file_ext(tax.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            tree.data = input$tree
            ext4 <- tools::file_ext(tree.data$datapath)
            
            req(otu.table, tax.table, sam.data, tree.data)
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") & (ext2 == "txt" | ext2 == "tsv") &
                (ext3 == "txt" | ext3 == "csv") & (ext4 == "tre" | ext4 == "nwk")) {
              otu.table.path = otu.table$datapath
              tax.table.path = tax.table$datapath
              sam.data.path = sam.data$datapath
              tree.data.path = tree.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              } else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext2 == "txt") {
                tax.tab <- read.table(tax.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext2 == "tsv") {
                tax.tab <- read.table(tax.table.path, header=TRUE, sep="\t")
                tax.tab = preprocess.tax.tab(tax.tab)
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              if (ext4 == "tre") {
                tree <- read.tree(file = tree.data.path)
              } else if (ext4 == "nwk") {
                tree <- read.tree(file = tree.data.path) 
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              tax.tab <- tax_table(as.matrix(tax.tab))
              sam.dat <- sample_data(sam.dat)
              tree <- phy_tree(tree)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu(otu.tab, tax.tab, tree)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                     type = "error")
                  } else {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } else if (biom.check.otu(otu.tab, tax.tab, tree)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature table, taxonomic table and tree tip labels"),
                                   type = "error")
                } else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
              return(biomData)
            } else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
        } else {
          return(NULL)
        }
        
        if (is.null(dataInfile)) {
          infile$biom <- NULL
          infile$qc_biom <- NULL
          infile$rare_biom = NULL
        } else {
          infile$biom <- dataInfile()
          infile$qc_biom <- dataInfile()
          infile$rare_biom = NULL
        }
        
        updateTabsetPanel(session, "side_menu",
                          selected = "step2")
        rcol$selected = "lightblue"
        
        if (!is.null(infile$biom)) QC$resume()
      })
    shinyjs::enable("Load_Individual_Data")
  })
  
  ######################################
  # Quality control and transformation #
  ######################################
  # This reactive expression stores the input data from either the individual data or phyloseq data
  QC = observe(suspended = T,{
    taxa.results$lib.size <- lib.size.func(infile$biom)$lib.size
    
    # Plots graphs using example infile data
    output$hist <- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      plot_ly(x = ~lib_size, nbinsx = input$binwidth,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Library Size", zeroline = FALSE))
    })
    
    output$hist2 <- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      plot_ly(x = ~mean_prop, nbinsx = input$binwidth2,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Mean Proportion", zeroline = FALSE))
    })
    
    output$boxplot<- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      
      plot_ly(x = ~lib_size, type = "box", notched=TRUE, name = "Library Size",
              color = ~"lib_size", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    output$boxplot2<- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      
      plot_ly(x = ~mean_prop, type = "box", notched=TRUE, name = "Mean Proportion",
              color = ~"mean_prop", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    ## Number of Taxonomic Rank for biom before QC
    num_tax.rank = reactive({
      tax.tab = tax_table(infile$qc_biom)
      num.tax.rank(tax.tab)
    })
    
    ## Fills value boxes using example biom data
    output$sample_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.sams), style = "font-size: 75%;"),
        "Sample Size", icon = icon("user-circle"), color = "fuchsia")
    })
    
    output$OTUs_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.otus), style = "font-size: 75%;"),
        "Number of Features", icon = icon("dna"), color = "aqua")
    })
    
    output$phyla <- renderValueBox({
      num.phyla = num_tax.rank()[1]
      valueBox(
        value = tags$p(paste0(num.phyla), style = "font-size: 75%;"),
        "Number of Phyla", icon = icon("sitemap"), color = "orange")
    })
    
    output$classes <- renderValueBox({
      num.classes = num_tax.rank()[2]
      valueBox(
        value = tags$p(paste0(num.classes), style = "font-size: 75%;"),
        "Number of Classes", icon = icon("sitemap"), color = "purple")
    })
    
    output$orders <- renderValueBox({
      num.orders = num_tax.rank()[3]
      valueBox(
        value = tags$p(paste0(num.orders), style = "font-size: 75%;"),
        "Number of Orders", icon = icon("sitemap"), color = "blue")
    })
    
    output$families <- renderValueBox({
      num.families = num_tax.rank()[4]
      valueBox(
        value = tags$p(paste0(num.families), style = "font-size: 75%;"),
        "Number of Families", icon = icon("sitemap"), color = "red")
    })
    
    output$genera <- renderValueBox({
      num.genera = num_tax.rank()[5]
      valueBox(
        value = tags$p(paste0(num.genera), style = "font-size: 75%;"),
        "Number of Genera", icon = icon("sitemap"), color = "lime")
    })
    
    output$species <- renderValueBox({
      num.species = num_tax.rank()[6]
      valueBox(
        value = tags$p(paste0(num.species), style = "font-size: 75%;"),
        "Number of Species", icon = icon("sitemap"), color = "teal" )
    })
    
    ## This event handler checks whether there is an input file and updates the slider options accordingly
    maxi.slider1 = as.numeric(lib.size.func(infile$qc_biom)$lib.size.sum["3rd quartile"])
    max.mean.prop = as.numeric(mean.prop.func(infile$qc_biom)$mean.prop.sum["3rd quartile"])
    maxi.slider2 = round(max.mean.prop, digits = 6)
    
    if (maxi.slider2 < 2e-05) {
      maxi.slider2 = 2e-05
    }
    
    updateSliderInput(session, "slider1", min = 0, max = round(maxi.slider1,-3))
    updateSliderInput(session, "slider2", min = 0, max = maxi.slider2*100)
  })
  
  ######################################
  ######################################
  #########  Data Analysis   ###########
  ######################################
  ######################################
  observeEvent(chooseData$alpha.div,{
    
    ######################################
    ######################################
    ######### Alpha diversity ############
    ######################################
    ######################################
    
    ######################################
    ### Cross-sectional Data Analysis ####
    ######################################
    output$primvars <- renderUI({
      tagList(
        selectInput("primvar", label = h4(strong("Primary Variable?", style = "color:black")),
                    c("Choose one" = "", chooseData$prim_vars), selected = chooseData$prim_vars[1], width = '70%'))
    })
    
    output$prim_vars_types <- renderUI({
      tagList(
        uiOutput("morePrimvar_opt"))
    })
    
    observeEvent(input$primvar,{
      ## user selects whether to use rarefied or non rarefied biom data
      
      if (input$primvar %in% chooseData$prim_vars) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar)
        # if variable is binary, perform data analysis
        if (is.results$result == "Binary") {
          alpha.categos$cat1 = alpha.bin.cat.func(chooseData$sam.dat, input$primvar)[1]
          alpha.categos$cat2 = alpha.bin.cat.func(chooseData$sam.dat, input$primvar)[2]
          
          output$morePrimvar_opt <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")),
              p("You can rename the categories of primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("alphaCat1", label = (paste0("Reference: ",alpha.categos$cat1)), value = alpha.categos$cat1, width = '80%'),
              textInput("alphaCat2", label = (paste0("Comparison: ",alpha.categos$cat2)), value = alpha.categos$cat2, width = '80%'))
          }) 
          
          ntselected.prim_vars = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar)
          output$covariates <- renderUI({
            tagList(
              prettyRadioButtons("covariates",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variables", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_vars, width = '70%'))),
              
              uiOutput("chooseTest"),
              
              prettyRadioButtons("chooseAdjustment", label = h4(strong("Multiple Testing Adjustment?", style = "color:black")), c("Yes", "No (Default)"), selected = "No (Default)",
                                 icon = icon("check"), animation = "jelly", status = "primary", width = '80%'),
              
              actionButton("runbtn_bin", (strong("Run!")), class = "btn-info")
            )
          })
          
          observeEvent(input$covariates,{
            if (input$covariates == "Covariate(s)") {
              
              shinyjs::show("covariates_variables")
              
              observeEvent(input$covariatesOptions,{
                if (!is.null(input$covariatesOptions)) {
                  output$chooseTest <- renderUI({
                    tagList(
                      selectInput("chooseMethod", label = h4(strong("Method?", style = "color:black")),
                                  c("Choose one" = "", "Linear regression", "Logistic regression"), selected = "Linear regression", width = '80%'),
                      div(id = "tauopt",
                          uiOutput("tauSlider")))
                  })
                } else {
                  output$chooseTest <- renderUI({
                    tagList(
                      selectInput("chooseMethod", label = h4(strong("Method?", style = "color:black")),
                                  c("Choose one" = "", "Welch t-test", "Wilcoxon rank-sum test"), selected = "Welch t-test", width = '80%'))
                  })
                }
              })
            } else if (input$covariates == "None") {
              shinyjs::hide("covariates_variables")
              output$chooseTest <- renderUI({
                tagList(
                  selectInput("chooseMethod", label = h4(strong("Method?", style = "color:black")), 
                              c("Choose one" = "", "Welch t-test", "Wilcoxon rank-sum test"), selected = "Welch t-test", width = '80%'))
              })
            }
          })
        } else if (is.results$result == "Continuous") {
          
          #if variable is continous, perform data analysis
          output$morePrimvar_opt <- renderUI({
            tagList(
              h4(strong("Rename Primary Variable?", style = "color:black")),
              p("You can rename the primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"),
              textInput("rename.con.var", label = NULL, value = input$primvar, width = '80%'))
          })
          
          ntselected.prim_vars = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar)
          output$covariates <- renderUI({
            tagList(
              prettyRadioButtons("covariates_cont",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variables", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions_cont"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_vars, width = '70%'))),
              
              selectInput("chooseMethod_cont", label = h4(strong("Method?", style = "color:black")),
                          c("Choose one" = "", "Linear regression"), selected = "Linear regression", width = '80%'),
              
              prettyRadioButtons("chooseAdjustment_cont", label = h4(strong("Multiple Testing Adjustment?", style = "color:black")), c("Yes", "No (Default)"), selected = "No (Default)",
                                 icon = icon("check"), animation = "jelly", status = "primary", width = '80%'),
              
              actionButton("runbtn_cont", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$covariates_cont,{
            if (input$covariates_cont == "Covariate(s)") {shinyjs::show("covariates_variables")
            } else if (input$covariates_cont == "None") {shinyjs::hide("covariates_variables")
            }
          })
        }
      }
    })
    
    ######################################
    ###### Longitudinal Data Analysis ####
    ######################################
    output$primvars_long <- renderUI({
      tagList(
        selectInput("primvarlong", label = h4(strong("Primary Variable?", style = "color:black")),
                    c("Choose one" = "", chooseData$prim_vars), selected = chooseData$prim_vars[1], width = '70%'))
    })
    
    output$prim_vars_types_long <- renderUI({
      tagList(
        uiOutput("morePrimvar_optlong"))
    })
    
    observeEvent(input$primvarlong,{
      ## user selects whether to use rarefied or non rarefied biom data
      if (input$primvarlong %in% chooseData$prim_vars) {
        is.results.long$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvarlong)
        # if variable is binary, perform data analysis
        if (is.results.long$result == "Binary") {
          
          alpha.categos.long$cat1 = alpha.bin.cat.func(chooseData$sam.dat, input$primvarlong)[1]
          alpha.categos.long$cat2 = alpha.bin.cat.func(chooseData$sam.dat, input$primvarlong)[2]
          
          ntselected.prim_varslong = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvarlong)
          
          output$morePrimvar_optlong <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")), 
              p("You can rename the categories of primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"),
              textInput("alphaCat1long", label = (paste0("Reference: ",alpha.categos.long$cat1)), value = alpha.categos.long$cat1, width = '80%'),
              textInput("alphaCat2long", label = (paste0("Comparison: ",alpha.categos.long$cat2)), value = alpha.categos.long$cat2, width = '80%'),
              
              br(),
              p(" ", style = "margin-bottom: -20px;"),
              
              h4(strong("Cluster Variable?", style = "color:black")), 
              p("Please select the cluster (group) variable. An example cluster variable contains subject IDs for repeated measures
              designs or family IDs for family-based studies.", style = "font-size:11pt"),
              prettyRadioButtons("clustervar",label = NULL, status = "primary", icon = icon("check"), animation = "jelly",
                                 ntselected.prim_varslong, selected = ntselected.prim_varslong[1], width = '70%'))
          })
          
          output$covariates_long <- renderUI({
            ntselected.clustervar = ntselected.prim_varslong[which(ntselected.prim_varslong != input$clustervar)]
            
            tagList(
              prettyRadioButtons("covariateslong",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variableslong", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptionslong"," Please select covariate(s)", status = "primary",
                                               ntselected.clustervar, width = '70%'))),
              
              uiOutput("chooseTestLong"),
              
              prettyRadioButtons("chooseAdjustmentlong", label = h4(strong("Multiple Testing Adjustment?", style = "color:black")), c("Yes", "No (Default)"), selected = "No (Default)",
                                 icon = icon("check"), animation = "jelly", status = "primary", width = '80%'),
              
              actionButton("runbtnlong_bin", (strong("Run!")), class = "btn-info")
            )
          })
          
          observeEvent(input$covariateslong,{
            if (input$covariateslong == "Covariate(s)") {
              shinyjs::show("covariates_variableslong")
              observeEvent(input$covariatesOptionslong,{
                if (!is.null(input$covariatesOptionslong)) {
                  output$chooseTestLong <- renderUI({
                    tagList(
                      selectInput("chooseMethodlong", label = h4(strong("Method?", style = "color:black")),
                                  c("Choose one" = "", "LMM", "GEE (Binomial)", "GLMM (Binomial)"), selected = "LMM", width = '80%'))
                  })
                } else {
                  output$chooseTestLong <- renderUI({
                    tagList(
                      selectInput("chooseMethodlong", label = h4(strong("Method?", style = "color:black")),
                                  c("Choose one" = "", "LMM", "GEE (Binomial)", "GLMM (Binomial)"), selected = "LMM", width = '80%'))
                  })
                }
              })
            } else if (input$covariateslong == "None") {
              shinyjs::hide("covariates_variableslong")
              output$chooseTestLong <- renderUI({
                tagList(
                  selectInput("chooseMethodlong", label = h4(strong("Method?", style = "color:black")),
                              c("Choose one" = "", "LMM", "GEE (Binomial)", "GLMM (Binomial)"), selected = "LMM", width = '80%'))
              })
            }
          })
        } else if (is.results.long$result == "Continuous") {
          
          ntselected.prim_vars.con_long = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvarlong)
          
          output$morePrimvar_optlong <- renderUI({
            tagList(
              h4(strong("Rename Primary Variable?", style = "color:black")),
              p("You can rename the primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("rename.con.varlong", label = NULL, value = input$primvarlong, width = '80%'),
              
              br(),
              p(" ", style = "margin-bottom: -20px;"),
              
              h4(strong("Cluster Variable?", style = "color:black")),
              p("Please select the cluster (group) variable. An example cluster variable contains subject IDs for repeated measures
              designs or family IDs for family-based studies.", style = "font-size:11pt"),
              prettyRadioButtons("clustervar.cont",label = NULL, status = "primary", icon = icon("check"), animation = "jelly",
                                 ntselected.prim_vars.con_long, selected = ntselected.prim_vars.con_long[1], width = '70%'))
          })
          
          output$covariates_long <- renderUI({
            
            ntselected.clustervar.cont = ntselected.prim_vars.con_long[which(ntselected.prim_vars.con_long != input$clustervar.cont)]
            
            tagList(
              prettyRadioButtons("covariates_contlong",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variableslong", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions_contlong"," Please select covariate(s)", status = "primary",
                                               ntselected.clustervar.cont, width = '70%'))),
              
              selectInput("chooseMethodlong_cont", label = h4(strong("Method?", style = "color:black")),
                          c("Choose one" = "", "LMM"), selected = "LMM", width = '80%'),
              
              prettyRadioButtons("chooseAdjustment_contlong", label = h4(strong("Multiple Testing Adjustment?", style = "color:black")), c("Yes", "No (Default)"), selected = "No (Default)",
                                 icon = icon("check"), animation = "jelly", status = "primary", width = '80%'),
              
              actionButton("runbtnlong_cont", (strong("Run!")), class = "btn-info"))
          })
          
          ntselected.prim_vars = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvarlong)
          
          observeEvent(input$covariates_contlong,{
            if (input$covariates_contlong == "Covariate(s)") {shinyjs::show("covariates_variableslong")
            } else if (input$covariates_contlong == "None") {shinyjs::hide("covariates_variableslong")
            }
          })
        }
      }
    })
    ################################################################################################
    
    ######################################
    ######################################
    ########## Beta diversity ############
    ######################################
    ######################################
    
    ######################################
    #### Cross-Sectional Data Analysis ###
    ######################################
    
    ## Beta Cross sectional UI options
    ## Beta primary variable UI option
    output$beta_primvar_cross <- renderUI({
      tagList(
        selectInput("beta.primvar_cross", label = h4(strong("Primary Variable?", style = "color:black")),
                    c("Choose one" = "", chooseData$prim_vars), selected = chooseData$prim_vars[1], width = '70%'))
    })
    
    output$beta_prim_vars_types_cross <- renderUI({
      tagList(
        uiOutput("beta_morePrimvar_optcross"))
    })
    
    observeEvent(input$beta.primvar_cross,{
      
      if (input$beta.primvar_cross %in% chooseData$prim_vars) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_cross)
        
        if (is.results$result == "Binary") {
          beta.categos$cat1 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar_cross)[1]
          beta.categos$cat2 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar_cross)[2]
          
          output$beta_morePrimvar_optcross <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")),
              p("You can rename the categories of primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("betaCat1", label = (paste0("Reference: ",beta.categos$cat1)), value = beta.categos$cat1, width = '80%'),
              textInput("betaCat2", label = (paste0("Comparison: ",beta.categos$cat2)), value = beta.categos$cat2, width = '80%'))
          }) 
          
          ntselected.prim_varsbin = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_cross)
          
          output$beta_covariates_cross <- renderUI({
            tagList(
              prettyRadioButtons("beta_covariates_bin",label = h4(strong("Covariate(s)?", style = "color:black")), 
                                 status = "primary", icon = icon("check"), animation = "jelly",
                                 c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "beta_covariates_variables_bin", style = "margin-left: 2%",
                           prettyCheckboxGroup("beta_covariatesOptions_bin"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_varsbin, width = '70%'))),
              
              selectInput("beta_chooseMethod_bin", label = h4(strong("Method?", style = "color:black")),
                          c("Choose one" = "", "MiRKAT"), selected = "MiRKAT", width = '80%'),
              
              actionButton("beta_runbtn_cross_bin", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$beta_covariates_bin,{
            if (input$beta_covariates_bin == "Covariate(s)") {
              shinyjs::show("beta_covariates_variables_bin")
            } else if (input$beta_covariates_bin == "None") {
              shinyjs::hide("beta_covariates_variables_bin")
            }
          })
        } else if (is.results$result == "Continuous") {
          #if variable is continous, perform data analysis
          output$beta_morePrimvar_optcross <- renderUI({
            tagList(
              h4(strong("Rename Primary Variable?", style = "color:black")),
              p("You can rename the primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"),
              textInput("beta.rename.con.var", label = NULL, value = input$beta.primvar_cross, width = '80%'))
          })
          
          ntselected.prim_vars = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_cross)
          
          output$beta_covariates_cross <- renderUI({
            tagList(
              prettyRadioButtons("beta.covariates_cont",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "beta_covariates_variables_cont", style = "margin-left: 2%",
                           prettyCheckboxGroup("beta.covariatesOptions_cont"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_vars, width = '70%'))),
              
              selectInput("beta.chooseMethod_cont", label = h4(strong("Method?", style = "color:black")),
                          c("Choose one" = "", "MiRKAT"), selected = "MiRKAT", width = '80%'),
              
              actionButton("beta_runbtn_cross_cont", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$beta.covariates_cont,{
            if (input$beta.covariates_cont == "Covariate(s)") {
              shinyjs::show("beta_covariates_variables_cont")
            } else if (input$beta.covariates_cont == "None") {shinyjs::hide("beta_covariates_variables_cont")}
          })
        }
      }
    })
    
    ######################################
    ###### Longitudinal Data Analysis ####
    ######################################
    
    ## Beta Longitudinal UI options
    ## Beta primary variable UI option
    output$beta_primvars_long <- renderUI({
      tagList(
        selectInput("beta.primvar_long", label = h4(strong("Primary Variable?", style = "color:black")),
                    c("Choose one" = "", chooseData$prim_vars), selected = chooseData$prim_vars[1], width = '70%'))
    })
    
    output$beta_prim_vars_types_long <- renderUI({
      tagList(
        uiOutput("beta_morePrimvar_optlong"))
    })
    
    observeEvent(input$beta.primvar_long,{
      
      if (input$beta.primvar_long %in% chooseData$prim_vars) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_long)
        
        # if variable is binary, perform data analysis
        if (is.results$result == "Binary") {
          
          beta.categos.long$cat1 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar_long)[1]
          beta.categos.long$cat2 = beta.bin.cat.func(chooseData$sam.dat, input$beta.primvar_long)[2]
          
          ntselected.prim_varsbin = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_long)
          
          output$beta_morePrimvar_optlong <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")), 
              p("You can rename the categories of primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"),
              textInput("betaCat1.long", label = (paste0("Reference: ",beta.categos.long$cat1)), value = beta.categos.long$cat1, width = '80%'),
              textInput("betaCat2.long", label = (paste0("Comparison: ",beta.categos.long$cat2)), value = beta.categos.long$cat2, width = '80%'),
              
              br(),
              p(" ", style = "margin-bottom: -20px;"),
              
              h4(strong("Cluster Variable?", style = "color:black")), 
              p("Please select the cluster (group) variable. Example cluster variables are subject IDs for repeated measures designs 
                and family IDs for family-based studies.", style = "font-size:11pt"),
              prettyRadioButtons("clustervar_bin.long",label = NULL, status = "primary", icon = icon("check"), animation = "jelly",
                                 ntselected.prim_varsbin, selected = ntselected.prim_varsbin[1], width = '70%'))
          })
          
          output$beta_covariates_long <- renderUI({
            ntselected.clustervar_long = ntselected.prim_varsbin[which(ntselected.prim_varsbin != input$clustervar_bin.long)]
            
            tagList(
              prettyRadioButtons("beta_covariates_bin.long",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "beta_covariates_variables_binLong", style = "margin-left: 2%",
                           prettyCheckboxGroup("beta_covariatesOptions_bin.long"," Please select covariate(s)", status = "primary",
                                               ntselected.clustervar_long, width = '70%'))),
              
              selectInput("beta_chooseMethod_bin.long", label = h4(strong("Method?", style = "color:black")),
                          c("Choose one" = "", "GLMM-MiRKAT"), selected = "GLMM-MiRKAT", width = '80%'),
              
              actionButton("beta_runbtn_bin.long", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$beta_covariates_bin.long,{
            if (input$beta_covariates_bin.long == "Covariate(s)") {
              shinyjs::show("beta_covariates_variables_binLong")
            } else if (input$beta_covariates_bin.long == "None") {shinyjs::hide("beta_covariates_variables_binLong")}
          })
        } else if (is.results$result == "Continuous") {
          #if variable is continous, perform data analysis
          ntselected.prim_vars_con =  cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$beta.primvar_long)
          
          output$beta_morePrimvar_optlong <- renderUI({
            tagList(
              h4(strong("Rename Primary Variable?", style = "color:black")),
              p("You can rename the primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("beta.rename.con.long", label = NULL, value = input$beta.primvar_long, width = '80%'),
              
              br(),
              p(" ", style = "margin-bottom: -20px;"),
              
              h4(strong("Cluster Variable?", style = "color:black")),
              p("Please select the cluster (group) variable. Example cluster variables are subject IDs for 
                repeated measures designs and family IDs for family-based studies.", style = "font-size:11pt"),
              prettyRadioButtons("clustervar_con.long",label = NULL, status = "primary", icon = icon("check"), animation = "jelly",
                                 ntselected.prim_vars_con, selected = ntselected.prim_vars_con[1], width = '70%'))
          })
          
          output$beta_covariates_long <- renderUI({
            
            ntselected.clustervarLong = ntselected.prim_vars_con[which(ntselected.prim_vars_con != input$clustervar_con.long)]
            
            tagList(
              prettyRadioButtons("beta.covariates_contLong",label = h4(strong("Covariate(s)?", style = "color:black")), status = "primary", icon = icon("check"),
                                 animation = "jelly", c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "beta_covariates_variables_contLong", style = "margin-left: 2%",
                           prettyCheckboxGroup("beta.covariatesOptions_contLong"," Please select covariate(s)", status = "primary",
                                               ntselected.clustervarLong, width = '70%'))),
              
              selectInput("beta.chooseMethod_contLong", label = h4(strong("Method?", style = "color:black")),
                          c("Choose one" = "", "GLMM-MiRKAT"), selected = "GLMM-MiRKAT", width = '80%'),
              
              actionButton("beta.runbtn_contLong", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$beta.covariates_contLong,{
            if (input$beta.covariates_contLong == "Covariate(s)") {
              shinyjs::show("beta_covariates_variables_contLong")
            } else if (input$beta.covariates_contLong == "None") {shinyjs::hide("beta_covariates_variables_contLong")}
          })
        }
      }
      
    })
  })
  observeEvent(chooseData$taxa.out,{
    
    ################################################################################################
    
    ######################################
    ######################################
    ######## Taxonomic Analysis ##########
    ######################################
    ######################################
    
    output$primvars_taxa <- renderUI({
      tagList(
        prettyRadioButtons("dataType_taxa", label = h4(strong("Data Format?", style = "color:black")), animation = "jelly",
                           c("CLR (Default)", "Count", "Proportion"), selected = "CLR (Default)",width = '70%'),
        selectInput("primvar_taxa", label = h4(strong("Primary Variable?", style = "color:black")),
                    choices = chooseData$prim_vars, selected = chooseData$prim_vars[1], width = '70%'))
    })
    
    output$primvars_taxa.long <- renderUI({
      tagList(
        prettyRadioButtons("dataType_taxa.long", label = h4(strong("Data Format?", style = "color:black")), animation = "jelly",
                           c("CLR (Default)", "Count", "Proportion"), selected = "CLR (Default)",width = '70%'),
        selectInput("primvar_taxa.long", label = h4(strong("Primary Variable?", style = "color:black")),
                    choices = chooseData$prim_vars, selected = chooseData$prim_vars[1], width = '70%'))
    })
    
    observeEvent(input$dataType_taxa,{
      if (input$dataType_taxa == "Count") {
        taxa.types$dataType = "rare.count"
        taxa.types$regression = "Negative binomial regression"
      } else if (input$dataType_taxa == "Proportion") {
        taxa.types$dataType = "imp.prop"
        taxa.types$regression = "Beta regression"
      } else if (input$dataType_taxa == "CLR (Default)") {
        taxa.types$dataType = "clr"
        taxa.types$regression = "Linear regression"
      }
    })
    
    observeEvent(input$dataType_taxa.long,{
      if (input$dataType_taxa.long == "Count") {
        taxa.types$dataType = "rare.count"
        taxa.types$regression.long = "GLMM (Negative Binomial)"
      } else if (input$dataType_taxa.long == "Proportion") {
        taxa.types$dataType = "imp.prop"
        taxa.types$regression.long = "GLMM (Beta)"
      } else if (input$dataType_taxa.long == "CLR (Default)") {
        taxa.types$dataType = "clr"
        taxa.types$regression.long = "LMM"
      }
    })
    
    observeEvent(input$primvar_taxa,{
      if (input$primvar_taxa %in% chooseData$prim_vars) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa)
        
        if (is.results$result == "Binary") {
          taxa.categos$cat1 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa)[1]
          taxa.categos$cat2 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa)[2]
          
          output$morePrimvar_opt_taxa <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")), 
              p("You can rename the categories of primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"),
              textInput("taxaCat1", label = (paste0("Reference: ",taxa.categos$cat1)), value = taxa.categos$cat1, width = '80%'),
              textInput("taxaCat2", label = (paste0("Comparison: ",taxa.categos$cat2)), value = taxa.categos$cat2, width = '80%'))
          }) 
          
          ntselected.prim_vars_taxa = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa)
          
          output$covariates_taxa <- renderUI({
            tagList(
              prettyRadioButtons("covariates_taxa", label = h4(strong("Covariate(s)?", style = "color:black")), animation = "jelly",
                                 c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variables_taxa", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions_taxa"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_vars_taxa, width = '70%'))),
              
              uiOutput("chooseTest_taxa"),
              
              prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%'),
              
              actionButton("runbtn_bin_taxa", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$covariates_taxa,{
            if (input$covariates_taxa == "Covariate(s)") {
              
              shinyjs::show("covariates_variables_taxa")
              
              observeEvent(input$covariatesOptions_taxa,{
                if (!is.null(input$covariatesOptions_taxa)) {
                  output$chooseTest_taxa <- renderUI({
                    tagList(
                      selectInput("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")),
                                  c("Choose one" = "", taxa.types$regression, "Logistic regression"), selected = taxa.types$regression,
                                  width = '80%'))
                  })
                } else {
                  shinyjs::hide("covariates_variables_taxa")
                  output$chooseTest_taxa <- renderUI({
                    tagList(
                      selectInput("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), 
                                  c("Choose one" = "", "Welch t-test", "Wilcoxon rank-sum test",taxa.types$regression, "Logistic regression"),
                                  selected = "Welch t-test", width = '80%'))
                  })
                }
              })
              
            } else if (input$covariates_taxa == "None") {
              shinyjs::hide("covariates_variables_taxa")
              output$chooseTest_taxa <- renderUI({
                tagList(
                  selectInput("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")), 
                              c("Choose one" = "", "Welch t-test", "Wilcoxon rank-sum test",taxa.types$regression, "Logistic regression"),
                              selected = "Welch t-test", width = '80%'))
              })
            }
          })
          
        } else if (is.results$result == "Continuous") {
          
          output$morePrimvar_opt_taxa <- renderUI({
            tagList(
              h4(strong("Rename Primary Variable?", style = "color:black")),
              p("You can rename the primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("rename.con.var_taxa", label = NULL, value = input$primvar_taxa))
          })
          
          ntselected.prim_vars_taxa = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa)
          
          output$covariates_taxa <- renderUI({
            tagList(
              prettyRadioButtons("covariates_taxa", label = h4(strong("Covariate(s)?", style = "color:black")), animation = "jelly",
                                 c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variables_taxa", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions_taxa"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_vars_taxa, width = '70%'))),
              
              uiOutput("chooseTest_taxa.cont"),
              
              prettyRadioButtons("include_species.dend", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%'),
              
              actionButton("runbtn_cont_taxa", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$covariates_taxa,{
            if (input$covariates_taxa == "Covariate(s)") {
              shinyjs::show("covariates_variables_taxa")
            } else if (input$covariates_taxa == "None") {
              shinyjs::hide("covariates_variables_taxa")
            }
            output$chooseTest_taxa.cont <- renderUI({
              tagList(
                selectInput("chooseMethod_taxa", label = h4(strong("Method?", style = "color:black")),
                            c("Choose one" = "", taxa.types$regression), selected = taxa.types$regression, width = '80%'))
            })
          }) 
        }
      }
    })
    
    observeEvent(input$primvar_taxa.long,{
      if (input$primvar_taxa.long %in% chooseData$prim_vars) {
        is.results$result = is.bin.con.pri(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa.long)
        
        if (is.results$result == "Binary") {
          
          taxa.categos$cat1 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa.long)[1]
          taxa.categos$cat2 = taxa.bin.cat.func(chooseData$sam.dat, input$primvar_taxa.long)[2]
          
          ntselected.prim_vars_taxa = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa.long)
          
          output$morePrimvar_opt_taxa.long <- renderUI({
            tagList(
              h4(strong("Rename Categories?", style = "color:black")), 
              p("You can rename the categories of primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size: 11pt"),
              textInput("taxaCat1", label = (paste0("Reference: ", taxa.categos$cat1)), value = taxa.categos$cat1, width = '80%'),
              textInput("taxaCat2", label = (paste0("Comparison: ", taxa.categos$cat2)), value = taxa.categos$cat2, width = '80%'),
              
              br(),
              p(" ", style = "margin-bottom: -20px;"),
              
              h4(strong("Cluster Variable?", style = "color:black")), 
              p("Please select the cluster (group) variable. An example cluster variable contains subject IDs for repeated measures 
                  designs or family IDs for family-based studies.", style = "font-size:11pt"),
              prettyRadioButtons("clustervar_taxa",label = NULL, status = "primary", icon = icon("check"), animation = "jelly",
                                 ntselected.prim_vars_taxa, selected = ntselected.prim_vars_taxa[1], width = '70%'))
          })
          
          output$covariates_taxa.long <- renderUI({
            ntselected.prim_vars_taxa2 <- ntselected.prim_vars_taxa[which(ntselected.prim_vars_taxa != input$clustervar_taxa)]
            
            tagList(
              prettyRadioButtons("covariates_taxa.long", label = h4(strong("Covariate(s)?", style = "color:black")), animation = "jelly",
                                 c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variables_taxa.long", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions_taxa.long"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_vars_taxa2, width = '70%'))),
              
              uiOutput("chooseTest_taxa.long"),
              
              prettyRadioButtons("include_species.dend.long", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%'),
              
              actionButton("runbtn_bin_taxa.long", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$covariates_taxa.long,{
            if (input$covariates_taxa.long == "Covariate(s)") {
              shinyjs::show("covariates_variables_taxa.long")
            } else {
              shinyjs::hide("covariates_variables_taxa.long")
            }
            output$chooseTest_taxa.long <- renderUI({
              tagList(
                selectInput("chooseMethod_taxa.long", label = h4(strong("Method?", style = "color:black")),
                            c("Choose one" = "",
                              taxa.types$regression.long, "GLMM (Binomial)", "GEE (Binomial)"),
                            selected = taxa.types$regression.long,
                            width = '80%')
              )
            })
          })
          
        } else if (is.results$result == "Continuous") {
          
          ntselected.prim_vars_taxa = cov.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con, input$primvar_taxa.long)
          
          output$morePrimvar_opt_taxa.long <- renderUI({
            tagList(
              h4(strong("Rename Primary Variable?", style = "color:black")),
              p("You can rename the primary variable. MiCloud keeps up to 8 characters on graphs.", style = "font-size:11pt"), 
              textInput("rename.con.var_taxa", label = NULL, value = input$primvar_taxa.long),
              
              br(),
              p(" ", style = "margin-bottom: -20px;"),
              
              h4(strong("Cluster Variable?", style = "color:black")), 
              p("Please select the cluster (group) variable. An example cluster variable contains subject IDs for repeated measures 
                  designs or family IDs for family-based studies.", style = "font-family: 'Verdana'; font-size:10pt"),
              prettyRadioButtons("clustervar_taxa",label = NULL, status = "primary", icon = icon("check"), animation = "jelly",
                                 ntselected.prim_vars_taxa, selected = ntselected.prim_vars_taxa[1], width = '70%'))
          })
          
          output$covariates_taxa.long <- renderUI({
            ntselected.prim_vars_taxa2 <- ntselected.prim_vars_taxa[which(ntselected.prim_vars_taxa != input$clustervar_taxa)]
            
            tagList(
              prettyRadioButtons("covariates_taxa.long", label = h4(strong("Covariate(s)?", style = "color:black")), animation = "jelly",
                                 c("None", "Covariate(s)"), selected = "None",width = '70%'),
              
              shinyjs::hidden(
                shiny::div(id = "covariates_variables_taxa.long", style = "margin-left: 2%",
                           prettyCheckboxGroup("covariatesOptions_taxa.long"," Please select covariate(s)", status = "primary",
                                               ntselected.prim_vars_taxa2, width = '70%'))),
              
              uiOutput("chooseTest_taxa.cont.long"),
              
              prettyRadioButtons("include_species.dend.long", label = h4(strong("Taxonomic Ranks?", style = "color:black")), animation = "jelly",
                                 c("Phylum - Genus (Default)", "Phylum - Species"), selected = "Phylum - Genus (Default)",
                                 icon = icon("check"), width = '80%'),
              
              actionButton("runbtn_cont_taxa.long", (strong("Run!")), class = "btn-info"))
          })
          
          observeEvent(input$covariates_taxa.long,{
            if (input$covariates_taxa.long == "Covariate(s)") {
              shinyjs::show("covariates_variables_taxa.long")
            } else if (input$covariates_taxa.long == "None") {
              shinyjs::hide("covariates_variables_taxa.long")
            }
            output$chooseTest_taxa.cont.long <- renderUI({
              tagList(
                selectInput("chooseMethod_taxa.long", label = h4(strong("Method?", style = "color:black")),
                            c("Choose one" = "", taxa.types$regression.long),
                            selected = taxa.types$regression.long, width = '80%'))
            })
          })
        }
      }
    })
    
  })
  #########################################################################################################
  
  ######################################
  ######################################
  ##########  RUN BUTTONS   ############
  ######################################
  ######################################
  
  ##########################################
  ##                  QC                  ##
  ##########################################
  observeEvent(input$run, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Data Trimming in progress")
        
        if (nchar(input$part.rem.str) == 0) {
          rem.tax.complete <- rem.tax.d
          rem.tax.partial <- rem.tax.str.d
        } else {
          rem.tax.complete <- unique(c(unlist(strsplit(input$rem.str, split = ",")), rem.tax.d))
          rem.tax.partial <- unique(c(unlist(strsplit(input$part.rem.str, split = ",")), rem.tax.str.d))
        }
        
        tax.tab <- tax_table(infile$biom)
        
        if (input$kingdom != "all") {
          ind <- is.element(tax.tab[,1], input$kingdom)
          validate(
            if (sum(ind) == 0) {
              showNotification(h4(paste("Error: Please select valid Kingdom. Available kingdoms are:", 
                                        paste(c(na.omit(unique(tax.tab[,1])) ,"and all"), collapse = ", "))),
                               type = "error")
            } else {
              NULL
            }
          )
        }
        
        shinyjs::disable("run")
        shinyjs::disable("slider1")
        shinyjs::disable("slider2")
        shinyjs::disable("kingdom")
        shinyjs::disable("skip")
        shinyjs::disable("binwidth")
        shinyjs::disable("binwidth2")
        
        rcol$selected = "rgba(255, 0, 0, 0.6)"
        
        infile$qc_biom = biom.clean(infile$biom, 
                                    input$kingdom, 
                                    lib.size.cut.off = input$slider1, 
                                    mean.prop.cut.off = input$slider2/100,
                                    rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial)
        
        infile$qc_biomNA = biom.cleanNA(infile$biom, 
                                    input$kingdom, 
                                    lib.size.cut.off = input$slider1, 
                                    mean.prop.cut.off = input$slider2/100,
                                    rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial)
        
        incProgress(3/10, message = "Rarefying in progress")
        lib_size.sum = lib.size.func(infile$qc_biom)$lib.size.sum
        infile$rare_biom = rarefy.func(infile$qc_biom, 
                                       cut.off = lib_size.sum["Minimum"],
                                       multi.rarefy = 1)
        infile$rare_biomNA = rarefy.func(infile$qc_biomNA, 
                                         cut.off = lib_size.sum["Minimum"],
                                         multi.rarefy = 1)
        
        incProgress(2/10, message = "Saving File in progress")
        
        chooseData$sam.dat = sample_data(infile$qc_biom)
        chooseData$mon.sin.rev.bin.con = is.mon.sin.rev.bin.con(chooseData$sam.dat)
        chooseData$prim_vars = pri.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
        chooseData$tax.tab = tax_table(infile$rare_biom)
        chooseData$tax.tabNA = tax_table(infile$rare_biomNA)
        
        #library.size <- library.size[names(library.size) %in% rownames(rare.sam.dat)]
        
        output$moreControls <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:13pt"),
                h5("Data after Quality Control"),
                downloadButton("downloadData2", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                h5("Data after Quality Control and Rarefaction"),
                downloadButton("downloadData3", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                p("For your reference, you can download the data files above for the phyloseq object (biom.after.qc) after QC and
                      (rare.biom.after.qc) after QC and rarefaction.",
                  style = "font-size:11pt")
            )
          )
        })
        
        output$text <- renderText({"You are all set! You can proceed to data analysis!"})
        
        qc_biom = infile$qc_biom
        output$downloadData2 <- downloadHandler(
          filename = function() {
            paste("biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(qc_biom, file = file1)
          })
        
        rare_biom = infile$rare_biom
        output$downloadData3 <- downloadHandler(
          filename = function() {
            paste("rare.biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(rare_biom, file = file1)
          })
        
        incProgress(1/10, message = "Done")
        shinyjs::enable("run")
        shinyjs::enable("slider1")
        shinyjs::enable("slider2")
        shinyjs::enable("kingdom")
        shinyjs::enable("skip")
        shinyjs::enable("binwidth")
        shinyjs::enable("binwidth2")
      })
  })
  
  ##########################################
  ## Data Calculation(Alpha, Beta, Taxa)  ##
  ##########################################
  observeEvent(input$divCalcRun, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("divCalcRun")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        
        incProgress(3/10, message = "Calculating Diversity")
        
        chooseData$alpha.div.rare = alpha.v1.func(infile$rare_biom)
        chooseData$alpha.div.qc = alpha.v1.func(infile$qc_biom)
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        incProgress(3/10, message = "Calculating Distance")
        
        ds.Ks$res = Ds.Ks.func(infile$rare_biom, infile$qc_biom)
        
        
        output$divCalcDownload <- renderUI({
          tagList(
            box(title = strong("Download Diversity Data", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                p("You can download alpha- and beta-diversity data.",
                  style = "font-size:11pt"), 
                h5("Alpha Diversity"),
                downloadButton("alphaDiv", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                h5("Beta Diversity"),
                downloadButton("betaDiv", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        alpha.div = chooseData$alpha.div
        
        output$alphaDiv <- downloadHandler(
          filename = function() {
            paste("Alpha.Diversity.txt")
          },
          content = function(alpha.file) {
            write.table(chooseData$alpha.div, file = alpha.file, row.names = TRUE, col.names = TRUE, sep = "\t")
          })
        
        output$betaDiv <- downloadHandler(
          filename = function() {
            paste("Beta.Diversity.zip")
          },
          content <- function(fname) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Jaccard.txt", "Bray.Curtis.txt", "U.UniFrac.txt" ,"G.UniFrac.txt", "W.UniFrac.txt")
            
            write.table(as.data.frame(ds.Ks$res$Ds$Jaccard), file = "Jaccard.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$Bray.Curtis), file = "Bray.Curtis.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$U.UniFrac), file = "U.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$G.UniFrac), file = "G.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$W.UniFrac), file = "W.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=fname, files=dataFiles)
          }
        )
        incProgress(1/10, message = "Done")
        shinyjs::enable("divCalcRun")
      })
  })
  observeEvent(input$datTransRun, {
    
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("divCalcRun")
        incProgress(3/10, message = "Transformation in progress")
        
        rare.otu.tab <- otu_table(infile$rare_biom)
        rare.tax.tab <- tax_table(infile$rare_biom)
        no.rare.otu.tab <- otu_table(infile$qc_biom)
        no.rare.tax.tab <- tax_table(infile$qc_biom)
        
        rare.otu.tabNA <- otu_table(infile$rare_biomNA)
        rare.tax.tabNA <- tax_table(infile$rare_biomNA)
        no.rare.otu.tabNA <- otu_table(infile$qc_biomNA)
        no.rare.tax.tabNA <- tax_table(infile$qc_biomNA)
        
        chooseData$taxa.out = tax.trans(no.rare.otu.tab, no.rare.tax.tab, rare.otu.tab, rare.tax.tab)
        chooseData$taxa.outNA = tax.trans.na(no.rare.otu.tabNA, no.rare.tax.tabNA, rare.otu.tabNA, rare.tax.tabNA)
        chooseData$taxa.names.out = taxa.names.rank(chooseData$taxa.out[[1]])
        chooseData$tax.tab = rare.tax.tab
        chooseData$tax.tabNA = rare.tax.tabNA
        
        chooseData$NAadded <- add_NA(chooseData$taxa.outNA, chooseData$tax.tabNA)
        
        incProgress(3/10, message = "Transformation in progress")
        
        output$datTransDownload <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:15pt"),
                p("You can download taxonomic abundance data.",
                  style = "font-size:11pt"),
                h5("Count", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), "Count (Rarefied)"),
                downloadButton("taxadataCount", "Download", width = '50%', style = " background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataRareCount", "Download", width = '50%', style = " background-color: red3"), br(), 
                h5("Proportion", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&nbsp;'), "CLR"),
                downloadButton("taxadataProp", "Download", width = '50%', style = " background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataCLR", "Download", width = '50%', style = " background-color: red3"), br(), br(),
            )
          )
        })
        
        count_biom = chooseData$taxa.out$count
        rare_biom = chooseData$taxa.out$rare.count
        prop_biom = chooseData$taxa.out$prop
        clr_biom = chooseData$taxa.out$clr
        
        output$taxadataCount <- downloadHandler(
          
          filename = function() {
            paste("Count.Data.zip")
          },
          content = function(count.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(count_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=count.file, files=dataFiles)
          }
        )
        
        incProgress(3/10, message = "Saving")
        
        output$taxadataRareCount <- downloadHandler(
          
          filename = function() {
            paste("Rarefied.Count.Data.zip")
          },
          content = function(rare.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(rare_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=rare.file, files=dataFiles)
          }
        )
        output$taxadataProp <- downloadHandler(
          filename = function() {
            paste("Proportion.Data.zip")
          },
          content = function(prop.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(prop_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=prop.file, files=dataFiles)
          }
        )
        output$taxadataCLR <- downloadHandler(
          filename = function() {
            paste("CLR.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(clr_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        incProgress(1/10, message = "Done")
        
        shinyjs::enable("divCalcRun")
      })
  })
  
  ##########################################
  ## Alpha Cross-Sectional Data Analysis ##
  ##########################################
  observeEvent(input$runbtn_bin,{
    validate(
      if (input$covariates == "Covariate(s)" & is.null(input$covariatesOptions)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_bin")
    shinyjs::disable("chooseAdjustment")
    shinyjs::disable("primvar")
    shinyjs::disable("chooseMethod")
    shinyjs::disable("covariates")
    shinyjs::disable("covariatesOptions")
    shinyjs::disable("alphaCat1")
    shinyjs::disable("alphaCat2")
    shinyjs::disable("chooseData")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Renaming")
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        alpha.categors = c(alpha.categos$cat1, alpha.categos$cat2)
        alpha.bin_categos = c(input$alphaCat1, input$alphaCat2)
        
        rename.cats_ref = alpha.bin_categos[which(alpha.categors == alpha.categos$cat1)]
        rename.cats_com = alpha.bin_categos[which(alpha.categors != alpha.categos$cat1)]
        
        alpha.bin.cat.ref.ori.out <- alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar)
        sam_dat <- alpha.bin.cat.recode.func(chooseData$sam.dat, input$primvar, alpha.bin.cat.ref.ori.out,
                                             rename.cats_ref, rename.cats_com)
        
        incProgress(3/10, message = "Calculating")
        alpha.div.bin.out <- alpha.bin.cat.ref.func(input$primvar, rename.cats_ref,
                                                    rename.cats_com, sam_dat,
                                                    chooseData$alpha.div)
        
        alpha.results$bin.var <- alpha.div.bin.out$bin.var
        alpha.results$alpha_div <- alpha.div.bin.out$alpha.div
        alpha.results$alpha.bin.sum.out <- alpha.bin.sum.func(bin.var = alpha.div.bin.out$bin.var,
                                                              alpha.div = alpha.div.bin.out$alpha.div)
        
        if (input$chooseMethod == "Welch t-test" | input$chooseMethod == "Wilcoxon rank-sum test") {
          if (input$chooseMethod == "Welch t-test") {
            incProgress(3/10, message = "T test")
            t.test.out <- try(alpha.bin.t.test(alpha.results$bin.var, alpha.results$alpha_div), silent = TRUE)
            alpha.data.results$table.p.out = t.test.out
            alpha.data.results$data.q.out = try(q.func(t.test.out, method = "BH"), silent = TRUE)
          }
          else if (input$chooseMethod == "Wilcoxon rank-sum test") {
            incProgress(3/10, message = "Wilcoxon test")
            wilcox.test.out <- try(alpha.bin.wilcox.test(alpha.results$bin.var, alpha.results$alpha_div), silent = TRUE)
            alpha.data.results$table.p.out = wilcox.test.out
            alpha.data.results$data.q.out = try(q.func(wilcox.test.out, method = "BH"), silent = TRUE)
          }
          
          if (input$chooseAdjustment == "Yes") {
            multi.test$boolval = TRUE
            alpha.data.results$table.out = alpha.data.results$data.q.out
          }
          else if (input$chooseAdjustment == "No (Default)") {
            multi.test$boolval = FALSE
            alpha.data.results$table.out = alpha.data.results$table.p.out
          }
          
          incProgress(3/10, message = "Graph")
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Graphs for Alpha Diversity", style = "color:black"), 
                  align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                  plotOutput("box_plots", height = 850, width = 500)
              )
            )
          })
          output$box_plots = renderPlot({
            try(alpha.bin.hist(alpha.results$bin.var, alpha.results$alpha_div, alpha.data.results$data.q.out, multi.test$boolval), silent = TRUE)
          })
        }
        else if (input$chooseMethod == "Linear regression" | input$chooseMethod == "Logistic regression") {
          alpha.bin.cov.out <- try(alpha.bin.cov.cat.ref.func(input$primvar, rename.cats_ref,
                                                          rename.cats_com, input$covariatesOptions,
                                                          sam_dat, chooseData$alpha.div), silent = TRUE)
          alpha.reg.results$bin.var <- alpha.bin.cov.out$bin.var
          alpha.reg.results$cov.var <- alpha.bin.cov.out$cov.var
          alpha.reg.results$alpha.div <- alpha.bin.cov.out$alpha.div
          
          if (input$chooseMethod == "Linear regression") {
            incProgress(3/10, message = "Linear regression with Covariate(s)")
            
            alpha.lm.bin.cov.out <- try(alpha.lm.bin.cov.func(bin.var = alpha.reg.results$bin.var,
                                                          cov.var = alpha.reg.results$cov.var,
                                                          alpha.div = alpha.reg.results$alpha.div,
                                                          scale = TRUE), silent = TRUE)
            alpha.data.results$table.p.out = alpha.lm.bin.cov.out
            alpha.data.results$data.q.out = try(q.func(alpha.lm.bin.cov.out, method = "BH"), silent = TRUE)
            alpha.data.results$table.out = alpha.data.results$data.q.out
          }
          
          else if (input$chooseMethod == "Logistic regression") {
            incProgress(3/10, message = "Logistic regression with Covariate(s)")
            
            alpha.logit.bin.cov.out <- try(alpha.logit.bin.cov.func(bin.var = alpha.reg.results$bin.var,
                                                                cov.var = alpha.reg.results$cov.var,
                                                                alpha.div = alpha.reg.results$alpha.div, 
                                                                scale = TRUE), silent = TRUE)
            alpha.logit.bin.cov.q.out = try(q.func(alpha.logit.bin.cov.out, method = "BH"), silent = TRUE)
            
            alpha.logit.reg.coef.bin.cov.out <- try(alpha.logit.reg.coef.bin.cov.func(bin.var = alpha.reg.results$bin.var,
                                                                                  cov.var = alpha.reg.results$cov.var,
                                                                                  alpha.div = alpha.reg.results$alpha.div, 
                                                                                  scale = TRUE), silent = TRUE)
            alpha.logit.reg.coef.bin.cov.q.out <- try(q.func(alpha.logit.reg.coef.bin.cov.out, method = "BH"), silent = TRUE)
            
            alpha.data.results$table.p.out = alpha.logit.reg.coef.bin.cov.out
            alpha.data.results$data.q.out = try(q.func(alpha.logit.reg.coef.bin.cov.out, method = "BH"), silent = TRUE)
            alpha.data.results$table.out = alpha.data.results$data.q.out
          }
          
          if (input$chooseAdjustment == "Yes") {
            multi.test$boolval = TRUE
            alpha.data.results$table.out = alpha.data.results$data.q.out
          }
          else if (input$chooseAdjustment == "No (Default)") {
            multi.test$boolval = FALSE
            alpha.data.results$table.out = alpha.data.results$table.p.out
          }
          
          incProgress(3/10, message = "Graph")
          
          if (input$chooseMethod == "Logistic regression") {
            
            output$alpha_display_results = renderUI({
              tagList(
                tabBox(title = strong("Forest plot", style = "color:black"), width = NULL,
                       tabPanel("Coefficient", align = "center",
                                plotOutput("forest_plots", height = 350, width = 500),
                       )
                       ,
                       tabPanel("Odds Ratio", align = "center",
                                plotOutput("forest_plots.or", height = 350, width = 500),
                       )
                )
              )
            })
            
            output$forest_plots = renderPlot({
              try(alpha.forest.plot(alpha.data.results$data.q.out, multi.test$boolval), silent = TRUE)
            })
            
            output$forest_plots.or = renderPlot({
              try(alpha.logit.forest.plot(alpha.logit.bin.cov.q.out, multi.test$boolval), silent = TRUE)
            })
          } else {
            
            output$alpha_display_results = renderUI({
              tagList(
                box(title = strong("Graphs for Alpha Diversity", style = "color:black"), 
                    align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                    plotOutput("forest_plots", height = 350, width = 500)
                )
              )
            })
            
            output$forest_plots = renderPlot({
              try(alpha.forest.plot(alpha.data.results$data.q.out, multi.test$boolval), silent = TRUE)
            })
          }
        }
        
        output$alpha_downloadTable = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Summary Statistics"),
                downloadButton("downloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        output$downloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Alpha.Sum.Table.txt")
          },
          content = function(file) {
            write.table(alpha.results$alpha.bin.sum.out, file, sep="\t")
          }
        )
        output$downloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Alpha.DA.Output.txt")
          },
          content = function(file) {
            write.table(alpha.data.results$table.out, file)
          }
        )
        
        ref_string = REFERENCE_CHECK(method_name = isolate(input$chooseMethod), FDR = isolate(input$chooseAdjustment))
        if (is.null(ref_string)) {
          shinyjs::hide("alpha_references")
        } else {
          shinyjs::show("alpha_references")
          output$alpha_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        shinyjs::enable("runbtn_bin")
        shinyjs::enable("primvar")
        shinyjs::enable("chooseAdjustment")
        shinyjs::enable("chooseMethod")
        shinyjs::enable("covariates")
        shinyjs::enable("covariatesOptions")
        shinyjs::enable("alphaCat1")
        shinyjs::enable("alphaCat2")
        shinyjs::enable("chooseData")
      })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  observeEvent(input$runbtn_cont,{
    validate(
      if (input$covariates_cont == "Covariate(s)" & is.null(input$covariatesOptions_cont)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    
    
    shinyjs::disable("runbtn_cont")
    shinyjs::disable("chooseAdjustment_cont")
    shinyjs::disable("primvar")
    shinyjs::disable("chooseMethod_cont")
    shinyjs::disable("covariates_cont")
    shinyjs::disable("rename.con.var")
    shinyjs::disable("covariatesOptions_cont")
    shinyjs::disable("chooseData")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(1/10, message = "Rename")
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        if (input$covariates_cont == "None") {
          alpha.results.cont$alpha.con.out = alpha.con.recode.func(chooseData$sam.dat, input$primvar, 
                                                                   input$rename.con.var, chooseData$alpha.div)
          alpha.results.cont$alpha.table.out = apply(alpha.results.cont$alpha.con.out$alpha.div, 2, alpha.ind.sum.func)
          
          if (input$chooseMethod_cont =="Linear regression") {
            incProgress(3/10, message = "Linear regression without Covariate(s)")
            alpha.data.results$table.p_out = try(alpha.lm.con.func(alpha.results.cont$alpha.con.out$con.var, 
                                                               alpha.results.cont$alpha.con.out$alpha.div), silent = TRUE)
            alpha.data.results$data.q.out = try(q.func(alpha.data.results$table.p_out, method = "BH"), silent = TRUE)
          }
        }
        else if (input$covariates_cont == "Covariate(s)") {
          if (is.null(input$covariatesOptions_cont)) {
            alpha.results.cont$alpha.con.out = alpha.con.recode.func(chooseData$sam.dat, input$primvar, 
                                                                     input$rename.con.var, chooseData$alpha.div)
            alpha.results.cont$alpha.table.out = apply(alpha.results.cont$alpha.con.out$alpha.div, 2, alpha.ind.sum.func)
            
            if (input$chooseMethod_cont =="Linear regression") {
              incProgress(3/10, message = "Linear regression without Covariate(s)")
              alpha.data.results$table.p_out = try(alpha.lm.con.func(alpha.results.cont$alpha.con.out$con.var, 
                                                                 alpha.results.cont$alpha.con.out$alpha.div), silent = TRUE)
              alpha.data.results$data.q.out = try(q.func(alpha.data.results$table.p_out, method = "BH"), silent = TRUE)
            }
          } else {
            alpha.con.cov.out <- alpha.con.cov.recode.func(chooseData$sam.dat, input$primvar, 
                                                           input$covariatesOptions_cont, input$rename.con.var, chooseData$alpha.div)
            alpha.results.cont$alpha.table.out =  apply(alpha.con.cov.out$alpha.div, 2, alpha.ind.sum.func)
            if (input$chooseMethod_cont =="Linear regression") {
              incProgress(3/10, message = "Linear regression with Covariate(s)")
              alpha.data.results$table.p_out = try(alpha.lm.con.cov.func(alpha.con.cov.out$con.var, alpha.con.cov.out$cov.var,
                                                                     alpha.con.cov.out$alpha.div, scale = TRUE), silent = TRUE)
              alpha.data.results$data.q.out = try(q.func(alpha.data.results$table.p_out, method = "BH"), silent = TRUE)
            }
          }
        }
        
        if (input$chooseAdjustment_cont == "Yes") {
          multi.test$boolval = TRUE
          alpha.data.results$table.out = alpha.data.results$data.q.out
        }
        else if (input$chooseAdjustment_cont == "No (Default)") {
          multi.test$boolval = FALSE
          alpha.data.results$table.out = alpha.data.results$table.p_out
        }
        
        incProgress(3/10, message = "Graph")
        
        output$alpha_display_results = renderUI({
          tagList(
            box(title = strong("Graphs for Alpha Diversity", style = "color:black"), 
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("forest_plots.cont", height = 850, width = 500)
            )
          )
        })
        
        if (input$covariates_cont == "None") {
          output$forest_plots.cont = renderPlot({
            if (input$chooseMethod_cont == "Linear regression") {
              try(alpha.con.plot(alpha.results.cont$alpha.con.out, alpha.data.results$data.q.out, multi.test$boolval), silent = TRUE)
            }
          })
        }
        else if (input$covariates_cont == "Covariate(s)") {
          if (is.null(input$covariatesOptions_cont)) {
            output$forest_plots.cont = renderPlot({
              if (input$chooseMethod_cont == "Linear regression") {
                try(alpha.con.plot(alpha.results.cont$alpha.con.out, alpha.data.results$data.q.out, multi.test$boolval), silent = TRUE)
              }
            })
          } else {
            output$forest_plots.cont = renderPlot({
              try(alpha.forest.plot(alpha.data.results$data.q.out, multi.test$boolval), silent = TRUE)
            })
          }
        }
        
        incProgress(1/10, message = "Save")
        output$downloadTabl_cont <- downloadHandler(
          filename = function() {
            paste("Alpha.Sum.Table.txt")
          },
          content = function(file) {
            write.table(alpha.results.cont$alpha.table.out, file, sep="\t")
          }
        )
        
        output$downloadTabl_cont2 <- downloadHandler(
          filename = function() {
            paste("Alpha.DA.Output.txt")
          },
          content = function(file) {
            write.table(alpha.data.results$table.out, file, sep="\t")
          }
        )
        ref_string = REFERENCE_CHECK(method_name = (input$chooseMethod_cont), FDR = isolate(input$chooseAdjustment_cont))
        if (is.null(ref_string)) {
          shinyjs::hide("alpha_references")
        } else {
          shinyjs::show("alpha_references")
          output$alpha_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      })
    
    
    shinyjs::enable("runbtn_cont")
    shinyjs::enable("chooseAdjustment_cont")
    shinyjs::enable("primvar")
    shinyjs::enable("chooseMethod_cont")
    shinyjs::enable("covariates_cont")
    shinyjs::enable("rename.con.var")
    shinyjs::enable("covariatesOptions_cont")
    shinyjs::enable("chooseData")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ######################################
  ## Alpha Longitudinal Data Analysis ##
  ######################################
  observeEvent(input$runbtnlong_bin,{
    validate(
      if (input$covariateslong == "Covariate(s)" & is.null(input$covariatesOptionslong)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    validate(
      if (is.null(input$clustervar) & is.null(input$clustervar.cont)) {
        showNotification("Error: No cluster variable available.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtnlong_bin")
    shinyjs::disable("chooseAdjustmentlong")
    shinyjs::disable("primvarlong")
    shinyjs::disable("chooseMethodlong")
    shinyjs::disable("covariateslong")
    shinyjs::disable("covariatesOptionslong")
    shinyjs::disable("clustervar")
    shinyjs::disable("alphaCat1long")
    shinyjs::disable("alphaCat2long")
    shinyjs::disable("chooseDatalong")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Renaming")
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        alpha.categors.long = c(alpha.categos.long$cat1, alpha.categos.long$cat2)
        alpha.bin_categos.long = c(input$alphaCat1long, input$alphaCat2long)
        
        rename.catslong_ref = alpha.bin_categos.long[which(alpha.categors.long == alpha.categos.long$cat1)]
        rename.catslong_com = alpha.bin_categos.long[which(alpha.categors.long != alpha.categos.long$cat1)]
        
        sam_dat.long <- alpha.bin.cat.recode.func(chooseData$sam.dat, input$primvarlong, 
                                                  alpha.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvarlong),
                                                  rename.catslong_ref, rename.catslong_com)
        
        incProgress(2/10, message = "Calculating")
        
        alpha.bin.out.long <- alpha.bin.cat.ref.func(input$primvarlong, rename.catslong_ref,
                                                     rename.catslong_com, sam_dat.long,
                                                     chooseData$alpha.div)
        
        alpha.bin.id.out.long <- alpha.bin.id.cat.ref.func(input$primvarlong, rename.catslong_ref, rename.catslong_com,
                                                           sel.id.var = input$clustervar, sam_dat.long, chooseData$alpha.div)
        
        alpha.resultslong$bin.var <- alpha.bin.id.out.long$bin.var
        alpha.resultslong$id.var <- alpha.bin.id.out.long$id.var
        alpha.resultslong$alpha.div <- alpha.bin.id.out.long$alpha.div
        
        data.alphaBin_res$alpha.bin.sum.out <- alpha.bin.sum.func(alpha.bin.out.long$bin.var, alpha.bin.out.long$alpha.div)
        
        
        
        if (input$chooseMethodlong == "LMM") {
          if (input$covariateslong == "None") {
            incProgress(3/10, message = "LMM without Covariate(s)")
            
            data.alphaBin_res$table.p_outbin = try(alpha.lmer.bin.id.func(bin.var = alpha.resultslong$bin.var, 
                                                                      id.var = alpha.resultslong$id.var, 
                                                                      alpha.div = alpha.resultslong$alpha.div, scale = TRUE), silent = TRUE)
            data.alphaBin_res$data.output = try(q.func(data.alphaBin_res$table.p_outbin, method = "BH"), silent = TRUE)
            
          } else if (input$covariateslong == "Covariate(s)") {
            incProgress(3/10, message = "LMM with Covariate(s)")
            
            alpha.bin.id.cov.out <- try(alpha.bin.id.cov.cat.ref.func(input$primvarlong, rename.catslong_ref,
                                                                  rename.catslong_com, input$clustervar, 
                                                                  input$covariatesOptionslong,
                                                                  sam_dat.long, chooseData$alpha.div), silent = TRUE)
            
            data.alphaBin_res$table.p_outbin = try(alpha.lmer.bin.id.cov.func(bin.var =  alpha.bin.id.cov.out$bin.var, 
                                                                          id.var =  alpha.bin.id.cov.out$id.var,
                                                                          cov.var =  alpha.bin.id.cov.out$cov.var,
                                                                          alpha.div =  alpha.bin.id.cov.out$alpha.div, scale = TRUE), silent = TRUE)
            
            data.alphaBin_res$data.output = try(q.func(data.alphaBin_res$table.p_outbin, method = "BH"), silent = TRUE)
          }
        } else if (input$chooseMethodlong == "GEE (Binomial)") {
          incProgress(3/10, message = "GEE (Binomial) with Covariate(s)")
          if (input$covariateslong == "Covariate(s)") {
            alpha.bin.id.cov.out <- try(alpha.bin.id.cov.cat.ref.func(input$primvarlong, rename.catslong_ref,
                                                                  rename.catslong_com, input$clustervar, 
                                                                  input$covariatesOptionslong,
                                                                  sam_dat.long, chooseData$alpha.div), silent = TRUE)
            data.alphaBin_res$table.p_outbin = try(alpha.logit.reg.coef.bin.cov.gee.func(bin.var =  alpha.bin.id.cov.out$bin.var, 
                                                                                     id.var =  alpha.bin.id.cov.out$id.var,
                                                                                     cov.var =  alpha.bin.id.cov.out$cov.var,
                                                                                     alpha.div =  alpha.bin.id.cov.out$alpha.div, scale = TRUE), silent = TRUE)
            data.alphaBin_res$data.output = try(q.func(data.alphaBin_res$table.p_outbin, method = "BH"), silent = TRUE)
            
          } else if (input$covariateslong == "None") {
            data.alphaBin_res$table.p_outbin = try(alpha.logit.reg.coef.bin.gee.func(bin.var = alpha.resultslong$bin.var, 
                                                                                 id.var = alpha.resultslong$id.var, 
                                                                                 alpha.div = alpha.resultslong$alpha.div, scale = TRUE), silent = TRUE)
            data.alphaBin_res$data.output = try(q.func(data.alphaBin_res$table.p_outbin, method = "BH"), silent = TRUE)
          }
        } else if (input$chooseMethodlong == "GLMM (Binomial)") {
          
          incProgress(3/10, message = "GLMM (Binomial) with Covariate(s)")
          if (input$covariateslong == "Covariate(s)") {
            alpha.bin.id.cov.out <- try(alpha.bin.id.cov.cat.ref.func(input$primvarlong, rename.catslong_ref,
                                                                  rename.catslong_com, input$clustervar, 
                                                                  input$covariatesOptionslong,
                                                                  sam_dat.long, chooseData$alpha.div), silent = TRUE)
            
            data.alphaBin_res$table.p_outbin = try(alpha.logit.reg.coef.bin.cov.glmm.b.func(bin.var =  alpha.bin.id.cov.out$bin.var, 
                                                                                        id.var =  alpha.bin.id.cov.out$id.var,
                                                                                        cov.var =  alpha.bin.id.cov.out$cov.var,
                                                                                        alpha.div =  alpha.bin.id.cov.out$alpha.div, scale = TRUE), silent = TRUE)
            data.alphaBin_res$data.output = try(q.func(data.alphaBin_res$table.p_outbin, method = "BH"), silent = TRUE)
            
          } else if (input$covariateslong == "None") {
            data.alphaBin_res$table.p_outbin = try(alpha.logit.reg.coef.bin.glmm.b.func(bin.var = alpha.resultslong$bin.var, 
                                                                                    id.var = alpha.resultslong$id.var, 
                                                                                    alpha.div = alpha.resultslong$alpha.div, scale = TRUE), silent = TRUE)
            data.alphaBin_res$data.output = try(q.func(data.alphaBin_res$table.p_outbin, method = "BH"), silent = TRUE)
          }
        }
        
        incProgress(3/10, message = "Graph")
        output$alpha_display_resultslong = renderUI({
          tagList(
            box(title = strong("Graphs for Alpha Diversity", style = "color:black"), 
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("graph_plotslong", height = 350, width = 500)
            )
          )
        })
        
        incProgress(1/10, message = "Save")
        output$alpha_downloadTablelong = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), 
                width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Summary Statistics"),
                downloadButton("downloadTabl1long", "Download", width = '50%', 
                               style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl2long", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        if (input$chooseAdjustmentlong == "Yes") {
          multi.test.long$boolval = TRUE
          data.alphaBin_res$table.output = data.alphaBin_res$data.output
        } else if (input$chooseAdjustmentlong == "No (Default)") {
          multi.test.long$boolval = FALSE
          data.alphaBin_res$table.output = data.alphaBin_res$table.p_outbin
        }
        
        output$downloadTabl1long <- downloadHandler(
          filename = function() {
            paste("Alpha.Sum.Table.txt")
          },
          content = function(file) {
            write.table(data.alphaBin_res$alpha.bin.sum.out, file, sep="\t")
          }
        )
        
        output$downloadTabl2long <- downloadHandler(
          filename = function() {
            paste("Alpha.DA.Output.txt")
          },
          content = function(file) {
            write.table(data.alphaBin_res$table.output, file)
          }
        )
        
        
        output$graph_plotslong = renderPlot({
          try(alpha.forest.lmer.plot(data.alphaBin_res$data.output, multi.test.long$boolval), silent = TRUE)
        })
        
        ref_string = REFERENCE_CHECK(method_name = isolate(input$chooseMethodlong), FDR = isolate(input$chooseAdjustmentlong))
        if (is.null(ref_string)) {
          shinyjs::hide("alpha_references_long")
        } else {
          shinyjs::show("alpha_references_long")
          output$alpha_references_long = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        shinyjs::enable("runbtnlong_bin")
        shinyjs::enable("chooseAdjustmentlong")
        shinyjs::enable("primvarlong")
        shinyjs::enable("chooseMethodlong")
        shinyjs::enable("covariateslong")
        shinyjs::enable("covariatesOptionslong")
        shinyjs::enable("clustervar")
        shinyjs::enable("alphaCat1long")
        shinyjs::enable("alphaCat2long")
        shinyjs::enable("chooseDatalong")
      })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  observeEvent(input$runbtnlong_cont,{
    validate(
      if (input$covariates_contlong == "Covariate(s)" & is.null(input$covariatesOptions_contlong)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtnlong_cont")
    shinyjs::disable("chooseAdjustment_contlong")
    shinyjs::disable("primvarlong")
    shinyjs::disable("chooseMethodlong_cont")
    shinyjs::disable("covariatesOptions_contlong")
    shinyjs::disable("covariates_contlong")
    shinyjs::disable("rename.con.varlong")
    shinyjs::disable("clustervar.cont")
    shinyjs::disable("chooseDatalong")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(1/10, message = "Renaming")
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        alpha.con.id.out <- alpha.con.id.recode.func(chooseData$sam.dat, input$primvarlong, input$clustervar.cont,
                                                     input$rename.con.varlong, chooseData$alpha.div)
        
        
        data.results.cont_long$alpha.table.out = apply(alpha.con.id.out$alpha.div, 2, alpha.ind.sum.func)
        
        alpha.noncovs_res$con.var = alpha.con.id.out$con.var
        alpha.noncovs_res$id.var = alpha.con.id.out$id.var
        alpha.noncovs_res$alpha.div = alpha.con.id.out$alpha.div
        
        if (input$chooseMethodlong_cont == "LMM") {
          if (input$covariates_contlong == "None") {
            incProgress(3/10, message = "LMM without Covariate(s)")
            
            data.results.cont_long$table_p.out = try(alpha.lmer.con.id.func(con.var = alpha.noncovs_res$con.var, 
                                                                        id.var = alpha.noncovs_res$id.var,
                                                                        alpha.div = alpha.noncovs_res$alpha.div, scale = TRUE), silent = TRUE)
            
            data.results.cont_long$data.q.out = try(q.func(data.results.cont_long$table_p.out, method = "BH"), silent = TRUE)
            
          } else if (input$covariates_contlong == "Covariate(s)") {
            incProgress(3/10, message = "LMM with Covariate(s)")
            alpha.cov_res = alpha.con.id.cov.recode.func(chooseData$sam.dat, input$primvarlong, input$clustervar.cont, 
                                                         input$covariatesOptions_contlong, input$rename.con.varlong,
                                                         chooseData$alpha.div)
            
            data.results.cont_long$table_p.out = try(alpha.lmer.con.id.cov.func(con.var = alpha.cov_res$con.var, 
                                                                            cov.var = alpha.cov_res$cov.var, 
                                                                            id.var = alpha.cov_res$id.var,
                                                                            alpha.div = alpha.cov_res$alpha.div, scale = TRUE), silent = TRUE)
            data.results.cont_long$data.q.out = try(q.func(data.results.cont_long$table_p.out, method = "BH"), silent = TRUE)
            
          }
        }
        
        
        incProgress(3/10, message = "Graph")
        if (input$chooseAdjustment_contlong == "Yes") {
          multi.test.long$boolval = TRUE
          data.alphaBin_res$table.output = data.results.cont_long$data.q.out
        } else if (input$chooseAdjustment_contlong == "No (Default)") {
          multi.test.long$boolval = FALSE
          data.alphaBin_res$table.output = data.results.cont_long$table.p_outbin
        }
        
        output$alpha_display_resultslong = renderUI({
          tagList(
            box(title = strong("Graphs for Alpha Diversity", style = "color:black"), 
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("graph_plotslong", height = 350, width = 500)
            )
          )
        })
        
        output$alpha_downloadTablelong = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"),
                width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Summary Statistics"),
                downloadButton("downloadTabl1long", "Download", width = '50%', 
                               style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl2long", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        incProgress(1/10, message = "SAVE")
        output$downloadTabl1long <- downloadHandler(
          filename = function() {
            paste("Alpha.Sum.Table.txt")
          },
          content = function(file) {
            write.table(data.results.cont_long$alpha.table.out, file, sep="\t")
          }
        )
        
        output$downloadTabl2long <- downloadHandler(
          filename = function() {
            paste("Alpha.DA.Output.txt")
          },
          content = function(file) {
            write.table(data.alphaBin_res$table.output, file)
          }
        )
        
        output$graph_plotslong = renderPlot({
          try(alpha.forest.lmer.plot(data.results.cont_long$data.q.out, multi.test.long$boolval), silent = TRUE)
        })
        
        ref_string = REFERENCE_CHECK(method_name = isolate(input$chooseMethodlong_cont), FDR = isolate(input$chooseAdjustment_contlong))
        if (is.null(ref_string)) {
          shinyjs::hide("alpha_references_long")
        } else {
          shinyjs::show("alpha_references_long")
          output$alpha_references_long = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        shinyjs::enable("runbtnlong_cont")
        shinyjs::enable("chooseAdjustment_contlong")
        shinyjs::enable("primvarlong")
        shinyjs::enable("chooseMethodlong_cont")
        shinyjs::enable("covariatesOptions_contlong")
        shinyjs::enable("covariates_contlong")
        shinyjs::enable("rename.con.varlong")
        shinyjs::enable("clustervar.cont")
        shinyjs::enable("chooseDatalong")
      })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ##########################################
  ## Beta Cross-Sectional Data Analysis ###
  ##########################################
  observeEvent(input$beta_runbtn_cross_bin,{
    validate(
      if (input$beta_covariates_bin == "Covariate(s)" & is.null(input$beta_covariatesOptions_bin)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    
    shinyjs::disable("beta_runbtn_cross_bin")
    shinyjs::disable("beta.primvar_cross")
    shinyjs::disable("betaCat1")
    shinyjs::disable("betaCat2")
    shinyjs::disable("beta_covariates_bin")
    shinyjs::disable("beta_covariatesOptions_bin")
    shinyjs::disable("beta_chooseMethod_bin")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        incProgress(3/10, message = "Rename")
        
        beta.categors = c(beta.categos$cat1, beta.categos$cat2)
        beta.bin_categos = c(input$betaCat1, input$betaCat2)
        
        rename.catsbin_ref = beta.bin_categos[which(beta.categors == beta.categos$cat1)]
        rename.catsbin_com = beta.bin_categos[which(beta.categors != beta.categos$cat1)]
        
        beta.bin.cat.ref.ori.out <- beta.bin.cat.ref.ori.func(chooseData$sam.dat, input$beta.primvar_cross)
        beta.sam_dat.bin <- beta.bin.cat.recode.func(chooseData$sam.dat, input$beta.primvar_cross,
                                                     beta.bin.cat.ref.ori.out,
                                                     rename.catsbin_ref, rename.catsbin_com)
        if (input$beta_covariates_bin == "None") {
          if (input$beta_chooseMethod_bin == "MiRKAT") {
            incProgress(3/10, message = "MiRKAT without Covariate(s)")
            beta.bin.out <- beta.bin.cat.ref.func(input$beta.primvar_cross,
                                                  rename.catsbin_ref, rename.catsbin_com,
                                                  beta.sam_dat.bin, Ds.Ks = ds.Ks$res)
            beta.data.results$data.q.out <- beta.bin.out
          }
        } else if (input$beta_covariates_bin == "Covariate(s)") {
          if (is.null(input$beta_covariatesOptions_bin)) {
            if (input$beta_chooseMethod_bin == "MiRKAT") {
              incProgress(3/10, message = "MiRKAT without Covariate(s)")
              beta.bin.out <- beta.bin.cat.ref.func(input$beta.primvar_cross,
                                                    rename.catsbin_ref, rename.catsbin_com,
                                                    beta.sam_dat.bin, Ds.Ks = ds.Ks$res)
              beta.data.results$data.q.out <- beta.bin.out
            }
          } else {
            beta.bin.cov.out <- beta.bin.cov.cat.ref.func(input$beta.primvar_cross,
                                                          rename.catsbin_ref, rename.catsbin_com,
                                                          input$beta_covariatesOptions_bin, beta.sam_dat.bin,
                                                          Ds.Ks = ds.Ks$res)
            
            if (input$beta_chooseMethod_bin == "MiRKAT") {
              incProgress(3/10, message = "MiRKAT with Covariate(s)")
              beta.data.results$data.q.out <- beta.bin.cov.out
            }
          }
        }
        
        
        incProgress(3/10, message = "Plotting Graph")
        output$beta_display_results_cross = renderUI({
          tagList(
            box(title = strong("Graphs for Beta Diversity", style = "color:black"), 
                align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("beta_graph_plots.bin", height = 850, width = 500)
            )
          )
        })
        
        if (input$beta_covariates_bin == "None") {
          beta.down.results$CS = try(mirkat.bin(beta.data.results$data.q.out), silent = TRUE)
          output$beta_graph_plots.bin = renderPlot({
            try(isolate(mirkat.bin.plot(beta.down.results$CS, beta.data.results$data.q.out)), silent = TRUE)
          })
        } else if (input$beta_covariates_bin == "Covariate(s)") {
          beta.down.results$CS = try(mirkat.bin.cov(beta.data.results$data.q.out), silent = TRUE)
          output$beta_graph_plots.bin = renderPlot({
            try(isolate(mirkat.bin.cov.plot(beta.down.results$CS, beta.data.results$data.q.out)), silent = TRUE)
          })
        }
        
        output$beta_downloadTable = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("beta_downloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        incProgress(3/10, message = "SAVE")
        
        output$beta_downloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Beta.DA.Output.txt")
          },
          content = function(file) {
            out_temp = as.data.frame(unlist(beta.down.results$CS))
            rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac", "omnibus_p")
            colnames(out_temp) = "p-value"
            beta.down.results$CS = out_temp
            write.table(beta.down.results$CS, file, sep="\t")
          }
        )
        
        ref_string = REFERENCE_CHECK(method_name = isolate(input$beta_chooseMethod_bin))
        if (is.null(ref_string)) {
          shinyjs::hide("beta_references")
        } else {
          shinyjs::show("beta_references")
          output$beta_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      }
    )
    
    
    shinyjs::enable("beta_runbtn_cross_bin")
    shinyjs::enable("beta.primvar_cross")
    shinyjs::enable("betaCat1")
    shinyjs::enable("betaCat2")
    shinyjs::enable("beta_covariates_bin")
    shinyjs::enable("beta_covariatesOptions_bin")
    shinyjs::enable("beta_chooseMethod_bin")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  observeEvent(input$beta_runbtn_cross_cont,{
    validate(
      if (input$beta.covariates_cont == "Covariate(s)" & is.null(input$beta.covariatesOptions_cont)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("beta_runbtn_cross_cont")
    shinyjs::disable("beta.primvar_cross")
    shinyjs::disable("beta.rename.con.var")
    shinyjs::disable("beta.covariates_cont")
    shinyjs::disable("beta.covariatesOptions_cont")
    shinyjs::disable("beta_chooseMethod_cont")
    
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
                   incProgress(1/10, message = "Rename")
                   
                   if (input$beta.covariates_cont == "None") {
                     if (input$beta.chooseMethod_cont == "MiRKAT") {
                       incProgress(3/10, message = "MiRKAT without Covariate(s)")
                       beta.con.out <- beta.con.recode.func(sam.dat = chooseData$sam.dat, input$beta.primvar_cross,
                                                            input$beta.rename.con.var, Ds.Ks = ds.Ks$res)
                       beta.resultscont$beta.cont.out = beta.con.out
                       
                     }
                   } else if (input$beta.covariates_cont == "Covariate(s)") {
                     if (is.null(input$beta.covariatesOptions_cont)) {
                       if (input$beta.chooseMethod_cont == "MiRKAT") {
                         incProgress(3/10, message = "MiRKAT without Covariate(s)")
                         beta.con.out <- beta.con.recode.func(sam.dat = chooseData$sam.dat, input$beta.primvar_cross,
                                                              input$beta.rename.con.var, Ds.Ks = ds.Ks$res)
                         beta.resultscont$beta.cont.out = beta.con.out
                       }
                     } else {
                       req(input$beta.covariatesOptions_cont)
                       if (input$beta.chooseMethod_cont == "MiRKAT") {
                         incProgress(3/10, message = "MiRKAT with Covariate(s)")
                         beta.cont.cov.out <- beta.con.cov.recode.func(sam.dat = chooseData$sam.dat, input$beta.primvar_cross,
                                                                       input$beta.covariatesOptions_cont,
                                                                       input$beta.rename.con.var, Ds.Ks = ds.Ks$res)
                         beta.resultscont$beta.cont.out = beta.cont.cov.out
                       }
                     }
                   }
                   
                   incProgress(3/10, message = "Plotting Graph")
                   output$beta_display_results_cross = renderUI({
                     tagList(
                       box(title = strong("Graphs for Beta Diversity", style = "color:black"), align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                           plotOutput("beta_graph_plots.cont", height = 850, width = 550)
                       )
                     )
                   })
                   
                   if (input$beta.covariates_cont == "None") {
                     beta.down.results$CS = try(mirkat.con(beta.resultscont$beta.cont.out), silent = TRUE)
                     output$beta_graph_plots.cont = renderPlot({
                       try(isolate(mirkat.con.plot(beta.down.results$CS, beta.resultscont$beta.cont.out)), silent = TRUE)
                     })
                   } else if (input$beta.covariates_cont == "Covariate(s)") {
                     beta.down.results$CS = try(mirkat.con.cov(beta.resultscont$beta.cont.out), silent = TRUE)
                     output$beta_graph_plots.cont = renderPlot({
                       try(isolate(mirkat.con.cov.plot(beta.down.results$CS, beta.resultscont$beta.cont.out)), silent = TRUE)
                     })
                   }
                   
                   output$beta_downloadTable = renderUI({
                     tagList(
                       box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                           p("You can download the data analysis outputs.",
                             style = "font-size:11pt"), 
                           h5("Data Analysis Outputs"),
                           downloadButton("beta_downloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3")
                       )
                     )
                   })
                   incProgress(3/10, message = "SAVE")
                   output$beta_downloadTabl1 <- downloadHandler(
                     filename = function() {
                       paste("Beta.DA.Output.txt")
                     },
                     content = function(file) {
                       out_temp = as.data.frame(unlist(beta.down.results$CS))
                       rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac", "omnibus_p")
                       colnames(out_temp) = "p-value"
                       beta.down.results$CS = out_temp
                       write.table(beta.down.results$CS, file, sep="\t")
                     }
                   )
                   
                   ref_string = REFERENCE_CHECK(method_name = isolate(input$beta.chooseMethod_cont))
                   if (is.null(ref_string)) {
                     shinyjs::hide("beta_references")
                   } else {
                     shinyjs::show("beta_references")
                     output$beta_references = renderUI({
                       tagList(
                         box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                             HTML(paste(ref_string, collapse="<br/>"))
                         )
                       )
                     })
                   }
                 })
    
    shinyjs::enable("beta_runbtn_cross_cont")
    shinyjs::enable("beta.primvar_cross")
    shinyjs::enable("beta.rename.con.var")
    shinyjs::enable("beta.covariates_cont")
    shinyjs::enable("beta.covariatesOptions_cont")
    shinyjs::enable("beta_chooseMethod_cont")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ######################################
  # Beta Longitudinal Data Analysis ####
  ######################################
  observeEvent(input$beta_runbtn_bin.long,{
    validate(
      if (input$beta_covariates_bin.long == "Covariate(s)" & is.null(input$beta_covariatesOptions_bin.long)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    validate(
      if (is.null(input$clustervar_bin.long) & is.null(input$clustervar_con.long)) {
        showNotification("Error: No cluster variable available.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("beta.primvar_long")
    shinyjs::disable("beta_runbtn_bin.long")
    shinyjs::disable("betaCat1.long")
    shinyjs::disable("betaCat2.long")
    shinyjs::disable("beta_chooseRef_bin.long")
    shinyjs::disable("clustervar_bin.long")
    shinyjs::disable("beta_covariates_bin.long")
    shinyjs::disable("beta_covariatesOptions_bin.long")
    shinyjs::disable("beta_chooseMethod_bin.long")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(1/10, message = "Rename")
        beta.categors = c(beta.categos.long$cat1, beta.categos.long$cat2)
        beta.bin_categos = c(input$betaCat1.long, input$betaCat2.long)
        
        rename.catsbin_long.ref = beta.bin_categos[which(beta.categors == beta.categos.long$cat1)]
        rename.catsbin_long.com = beta.bin_categos[which(beta.categors != beta.categos.long$cat1)]
        
        beta.bin.cat.ref.ori.out <- beta.bin.cat.ref.ori.func(chooseData$sam.dat, input$beta.primvar_long)
        beta.sam_dat.bin <- beta.bin.cat.recode.func(chooseData$sam.dat, input$beta.primvar_long,
                                                     beta.bin.cat.ref.ori.out,
                                                     rename.catsbin_long.ref, rename.catsbin_long.com)
        
        if (input$beta_covariates_bin.long == "None") {
          if (input$beta_chooseMethod_bin.long == "GLMM-MiRKAT") {
            incProgress(3/10, message = "GLMM-MiRKAT without Covariate(s)")
            beta.bin.out.res <-  beta.bin.id.cat.ref.func(input$beta.primvar_long,
                                                          rename.catsbin_long.ref, rename.catsbin_long.com,
                                                          input$clustervar_bin.long,
                                                          beta.sam_dat.bin, Ds.Ks = ds.Ks$res)
            
            beta.data.results_long$beta.bin.out <- beta.bin.out.res
          }
        } else if (input$beta_covariates_bin.long == "Covariate(s)") {
          if (is.null(input$beta_covariatesOptions_bin.long)) {
            if (input$beta_chooseMethod_bin.long == "GLMM-MiRKAT") {
              incProgress(3/10, message = "GLMM-MiRKAT without Covariate(s)")
              beta.bin.out.res <-  beta.bin.id.cat.ref.func(input$beta.primvar_long,
                                                            rename.catsbin_long.ref, rename.catsbin_long.com,
                                                            input$clustervar_bin.long,
                                                            beta.sam_dat.bin, Ds.Ks = ds.Ks$res)
              
              beta.data.results_long$beta.bin.out <- beta.bin.out.res
            }
          } else {
            if (input$beta_chooseMethod_bin.long == "GLMM-MiRKAT") {
              incProgress(3/10, message = "GLMM-MiRKAT with Covariate(s)")
              beta.bin.cov.out <- beta.bin.id.cov.cat.ref.func(input$beta.primvar_long,
                                                               rename.catsbin_long.ref, rename.catsbin_long.com,
                                                               input$clustervar_bin.long,
                                                               input$beta_covariatesOptions_bin.long, beta.sam_dat.bin,
                                                               Ds.Ks = ds.Ks$res)
              
              beta.data.results_long$beta.bin.out <- beta.bin.cov.out
            }
          }
        }
        
        incProgress(3/10, message = "Plotting Graph")
        output$beta_display_resultslong = renderUI({
          tagList(
            box(title = strong("Graphs for Beta Diversity", style = "color:black"), align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("beta_graph_plots.bin_long", height = 850, width = 550)
            )
          )
        })
        
        if (input$beta_covariates_bin.long == "None") {
          beta.down.results$LONG = try(glmm.mirkat.bin(beta.data.results_long$beta.bin.out), silent = TRUE)
          output$beta_graph_plots.bin_long = renderPlot({
            try(isolate(glmm.mirkat.bin.plot(beta.down.results$LONG,beta.data.results_long$beta.bin.out)), silent = TRUE)
          })
        } else if (input$beta_covariates_bin.long == "Covariate(s)") {
          beta.down.results$LONG = try(glmm.mirkat.bin.cov(beta.data.results_long$beta.bin.out), silent = TRUE)
          output$beta_graph_plots.bin_long = renderPlot({
            try(isolate(glmm.mirkat.bin.cov.plot(beta.down.results$LONG,beta.data.results_long$beta.bin.out)), silent = TRUE)
          })
        }
        
        output$beta_downloadTablelong = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("beta_downloadTabl1long", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        incProgress(3/10, message = "SAVE")
        output$beta_downloadTabl1long <- downloadHandler(
          filename = function() {
            paste("Beta.DA.Output.txt")
          },
          content = function(file) {
            out_temp = as.data.frame(unlist(beta.down.results$LONG))
            rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac", "aGLMM-MiRKAT")
            colnames(out_temp) = "p-value"
            beta.down.results$LONG = out_temp
            write.table(beta.down.results$LONG, file, sep="\t")
          }
        )
        ref_string = REFERENCE_CHECK(method_name = isolate(input$beta_chooseMethod_bin.long))
        if (is.null(ref_string)) {
          shinyjs::hide("beta_references_long")
        } else {
          shinyjs::show("beta_references_long")
          output$beta_references_long = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      })
    
    shinyjs::enable("beta.primvar_long")
    shinyjs::enable("beta_runbtn_bin.long")
    shinyjs::enable("betaCat1.long")
    shinyjs::enable("betaCat2.long")
    shinyjs::enable("beta_chooseRef_bin.long")
    shinyjs::enable("clustervar_bin.long")
    shinyjs::enable("beta_covariates_bin.long")
    shinyjs::enable("beta_covariatesOptions_bin.long")
    shinyjs::enable("beta_chooseMethod_bin.long")
    
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  observeEvent(input$beta.runbtn_contLong,{
    validate(
      if (input$beta.covariates_contLong == "Covariate(s)" & is.null(input$beta.covariatesOptions_contLong)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("beta.runbtn_contLong")
    shinyjs::disable("beta.primvar_long")
    shinyjs::disable("beta.rename.con.long")
    shinyjs::disable("clustervar_con.long")
    shinyjs::disable("beta.covariates_contLong")
    shinyjs::disable("beta.covariatesOptions_contLong")
    shinyjs::disable("beta.chooseMethod_contLong")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        if (input$beta.covariates_contLong == "None") {
          if (input$beta.chooseMethod_contLong == "GLMM-MiRKAT") {
            incProgress(3/10, message = "GLMM-MiRKAT without Covariate(s)")
            beta.con.id.out <-  beta.con.id.recode.func(sam.dat = chooseData$sam.dat, 
                                                        sel.con.var = input$beta.primvar_long, 
                                                        sel.id.var = input$clustervar_con.long, 
                                                        rename.con.var = input$beta.rename.con.long, Ds.Ks = ds.Ks$res)
            beta.resultscon_long$beta.con.out = beta.con.id.out
          }
        } else if (input$beta.covariates_contLong == "Covariate(s)") {
          if (is.null(input$beta.covariatesOptions_contLong)) {
            if (input$beta.chooseMethod_contLong == "GLMM-MiRKAT") {
              incProgress(3/10, message = "GLMM-MiRKAT without Covariate(s)")
              beta.con.id.out <-  beta.con.id.recode.func(sam.dat = chooseData$sam.dat, 
                                                          sel.con.var = input$beta.primvar_long, 
                                                          sel.id.var = input$clustervar_con.long, 
                                                          rename.con.var = input$beta.rename.con.long, Ds.Ks = ds.Ks$res)
              beta.resultscon_long$beta.con.out = beta.con.id.out
            }
          } else {
            if (input$beta.chooseMethod_contLong == "GLMM-MiRKAT") {
              incProgress(3/10, message = "GLMM-MiRKAT with Covariate(s)")
              beta.con.id.cov.out <-  beta.con.id.cov.recode.func(sam.dat = chooseData$sam.dat, 
                                                                  sel.con.var = input$beta.primvar_long, 
                                                                  sel.id.var = input$clustervar_con.long,
                                                                  sel.cov.var = input$beta.covariatesOptions_contLong, 
                                                                  rename.con.var = input$beta.rename.con.long, 
                                                                  Ds.Ks = ds.Ks$res)
              
              beta.resultscon_long$beta.con.out = beta.con.id.cov.out
            }
          }
        }
        
        incProgress(3/10, message = "Plotting Graph")
        output$beta_display_resultslong = renderUI({
          tagList(
            box(title = strong("Graphs for Beta Diversity", style = "color:black"), align = "center", width = NULL, status = "primary", solidHeader = TRUE,
                plotOutput("beta_graph_plots.conLong", height = 850, width = 550)
            )
          )
        })
        
        if (input$beta.covariates_contLong == "None") {
          beta.down.results$LONG = try(glmm.mirkat.con(beta.con.id.out = beta.resultscon_long$beta.con.out), silent = TRUE)
          output$beta_graph_plots.conLong = renderPlot({
            try(isolate(glmm.mirkat.con.plot(beta.down.results$LONG,beta.resultscon_long$beta.con.out)), silent = TRUE)
          })
        } else if (input$beta.covariates_contLong == "Covariate(s)") {
          beta.down.results$LONG = try(glmm.mirkat.con.cov(beta.con.id.cov.out = beta.resultscon_long$beta.con.out) , silent = TRUE)
          output$beta_graph_plots.conLong = renderPlot({
            try(isolate(glmm.mirkat.con.cov.plot(beta.down.results$LONG,beta.resultscon_long$beta.con.out)), silent = TRUE)
          })
        }
        
        output$beta_downloadTablelong = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("beta_downloadTabl1long", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        incProgress(3/10, message = "SAVE")
        output$beta_downloadTabl1long <- downloadHandler(
          filename = function() {
            paste("Beta.DA.Output.txt")
          },
          content = function(file) {
            out_temp = as.data.frame(unlist(beta.down.results$LONG))
            rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac", "aGLMM-MiRKAT")
            colnames(out_temp) = "p-value"
            beta.down.results$LONG = out_temp
            write.table(beta.down.results$LONG, file, sep="\t")
          }
        )
        ref_string = REFERENCE_CHECK(method_name = isolate(input$beta.chooseMethod_contLong))
        if (is.null(ref_string)) {
          shinyjs::hide("beta_references_long")
        } else {
          shinyjs::show("beta_references_long")
          output$beta_references_long = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
      })
    
    shinyjs::enable("beta.runbtn_contLong")
    shinyjs::enable("beta.primvar_long")
    shinyjs::enable("beta.rename.con.long")
    shinyjs::enable("clustervar_con.long")
    shinyjs::enable("beta.covariates_contLong")
    shinyjs::enable("beta.covariatesOptions_contLong")
    shinyjs::enable("beta.chooseMethod_contLong")
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ######################################
  # Taxa Cross-Sectional Data Analysis #
  ######################################
  observeEvent(input$runbtn_bin_taxa,{
    
    validate(
      if (input$covariates_taxa == "Covariate(s)" & is.null(input$covariatesOptions_taxa)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_bin_taxa")
    shinyjs::disable("dataType_taxa")
    shinyjs::disable("primvar_taxa")
    shinyjs::disable("chooseMethod_taxa")
    shinyjs::disable("covariates_taxa")
    shinyjs::disable("chooseRef_taxa")
    shinyjs::disable("taxaCat1")
    shinyjs::disable("taxaCat2")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        taxa.categors = c(taxa.categos$cat1, taxa.categos$cat2)
        taxa.bin_categos = c(input$taxaCat1, input$taxaCat2)
        
        rename.cats_ref = taxa.bin_categos[which(taxa.categors == taxa.categos$cat1)]
        rename.cats_com = taxa.bin_categos[which(taxa.categors != taxa.categos$cat1)]
        
        taxa.bin.cat.ref.ori.out <- taxa.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar_taxa)
        taxa.results$lib.size <- taxa.results$lib.size[names(taxa.results$lib.size) %in% rownames(chooseData$sam.dat)]
        
        if (input$include_species.dend == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        incProgress(1/10, message = "Examining Data in progress")
        sam_dat <- taxa.bin.cat.recode.func(chooseData$sam.dat, input$primvar_taxa, taxa.bin.cat.ref.ori.out,
                                            rename.cats_ref, rename.cats_com)
        
        if (input$covariates_taxa == "None") {
          
          if (input$chooseMethod_taxa == "Logistic regression") {
            
            taxa.bin.logit.out <- taxa.bin.cat.ref.logit.united.func(input$primvar_taxa, rename.cats_ref,
                                                                     rename.cats_com, sam_dat, taxa = chooseData$taxa.out[[taxa.types$dataType]])
            taxa.results$bin.var <- taxa.bin.logit.out$bin.var
            taxa.results$taxa <- taxa.bin.logit.out$taxa
            
          } else {
            if (input$chooseMethod_taxa == "Negative binomial regression") {
              taxa.bin.out <- taxa.bin.cat.ref.united.func(input$primvar_taxa, rename.cats_ref,
                                                           rename.cats_com, sam_dat, taxa = chooseData$taxa.out[["count"]])
              
            } else {
              taxa.bin.out <- taxa.bin.cat.ref.united.func(input$primvar_taxa, rename.cats_ref,
                                                           rename.cats_com, sam_dat, taxa = chooseData$taxa.out[[taxa.types$dataType]])
            }
            taxa.results$bin.var <- taxa.bin.out$bin.var
            taxa.results$taxa <- taxa.bin.out$taxa
          }
          
        } else if (input$covariates_taxa == "Covariate(s)") {
          if (input$chooseMethod_taxa == "Negative binomial regression") {
            taxa.bin.cov.out <- taxa.bin.cov.cat.ref.united.func(input$primvar_taxa, rename.cats_ref, 
                                                                 rename.cats_com, input$covariatesOptions_taxa, 
                                                                 sam_dat, chooseData$taxa.out[["count"]])
          } else {
            taxa.bin.cov.out <- taxa.bin.cov.cat.ref.united.func(input$primvar_taxa, rename.cats_ref, 
                                                                 rename.cats_com, input$covariatesOptions_taxa, 
                                                                 sam_dat, chooseData$taxa.out[[taxa.types$dataType]])
          }
          taxa.results$bin.var <- taxa.bin.cov.out$bin.var
          taxa.results$taxa <- taxa.bin.cov.out$taxa
          taxa.results$cov.var <- taxa.bin.cov.out$cov.var
          
        }
        taxa.results$taxa.bin.sum.out <- taxa.bin.sum.united.func(taxa.results$bin.var, taxa.results$taxa)
        
        taxa_dataBinvar <- taxa.results$bin.var
        taxa_dataTaxa <- taxa.results$taxa
        
        if (input$chooseMethod_taxa == "Welch t-test" | input$chooseMethod_taxa == "Wilcoxon rank-sum test") {
          if (input$chooseMethod_taxa == "Welch t-test") {
            incProgress(5/10, message = "Welch t-test")
            taxa.t.test.out <- taxa.bin.t.test.united(taxa.results$bin.var, taxa.results$taxa)
            taxa.t.test.q.out <- bin.q.united.func(taxa.t.test.out, method = "BH")
            
            taxa.outputs$DAoutput = taxa.t.test.q.out
            
            nrow = numeric()
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.t.test.q.out[[r]]$Q.value < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              } 
            }
          } else if (input$chooseMethod_taxa == "Wilcoxon rank-sum test") {
            incProgress(5/10, message = "Wilcoxon rank-sum test")
            taxa.wilcox.test.out <- taxa.bin.wilcox.test.united(taxa.results$bin.var, taxa.results$taxa)
            taxa.wilcox.test.q.out <- bin.q.united.func(taxa.wilcox.test.out, method = "BH")
            taxa.wilcox.test.est.added <- taxa.wilcox.test.est.func(taxa.results$bin.var, taxa.results$taxa, rename.cats_ref, rename.cats_com, taxa.wilcox.test.q.out)
            
            taxa.outputs$DAoutput = taxa.wilcox.test.est.added
            
            nrow = numeric()
            for (r in 1:6) {
              row.num <- ceiling(sum(taxa.outputs$DAoutput[[r]]$Q.value < 0.05)/4)
              if (row.num > 0) {
                nrow[r] <- row.num
              } else {
                nrow[r] <- 1
              }
            }
          }
          
          incProgress(3/10, message = "Displaying Results in progress")
          output$taxa_display = renderUI({
            tagList(
              tabBox(title = strong("Box Plot", style = "color:black"), width = NULL,
                     tabPanel("Phylum", align = "center",
                              plotOutput("rank1", height = nrow[1]*250, width = 750),
                     )
                     ,
                     tabPanel("Class", align = "center",
                              plotOutput("rank2", height = nrow[2]*250, width = 750),
                     )
                     ,tabPanel("Order", align = "center",
                               plotOutput("rank3", height = nrow[3]*250, width = 750),
                     )
                     ,tabPanel("Family", align = "center",
                               plotOutput("rank4", height = nrow[4]*250, width = 750),
                     )
                     ,tabPanel("Genus", align = "center",
                               plotOutput("rank5", height = nrow[5]*250, width = 750),
                     )
                     ,tabPanel("Species", align = "center",
                               plotOutput("rank6", height = nrow[6]*250, width = 750),
                     )
              )
            )
          })
          
          output$rank1 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 1, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank2 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 2, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank3 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 3, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank4 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 4, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank5 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 5, TRUE)  ####all.t.test.united should be used
          })
          
          output$rank6 = renderPlot({ 
            taxa.bin.boxplot(taxa_dataBinvar, taxa_dataTaxa, taxa.outputs$DAoutput, chooseData$taxa.names.out, 6, TRUE)  ####all.t.test.united should be used
          })
          
          taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)
          
          flow.text <- taxa.sig$flow.text
          taxon.tab <- taxa.sig$taxon.tab
          ci.tab.all <- taxa.sig$ci.tab.all
          
          if ( include == FALSE){
            ci.tab.all <- ci.tab.all[-(grep("s_", taxon.tab$Taxon, fixed = TRUE)+1)]
            taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
          }
          
          if ( length(ci.tab.all) > 1 ){
            
            for( i in 1:nrow(taxon.tab)){
              if ( ci.tab.all[-1][i] < 0){
                taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
                taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
              }
              else{
                taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
                taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
              }
              
            }
          }
          
          
          N <- dim(taxon.tab)[1]
          itr <- ceiling(N/5)
          tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
          tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
          tab.three <- data.frame( matrix(ncol=2,nrow=0) )
          tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
          tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
          
          colnames(tab.one) <- c("ID", "Taxon")
          colnames(tab.two) <- c("ID", "Taxon")
          colnames(tab.three) <- c("ID", "Taxon")
          colnames(tab.four) <- c("ID", "Taxon")
          colnames(tab.five) <- c("ID", "Taxon")
          
          if ( dim(taxon.tab)[1] > 0 ) {
            
            i = 1
            j = 1
            N.rmd <- N %% 5
            N.fix <- N + 5 - N.rmd
            while ( i <= itr) {
              tab.one  [i,] <- taxon.tab[j,]
              tab.two  [i,] <- taxon.tab[j+1,]
              tab.three[i,] <- taxon.tab[j+2,]
              tab.four [i,] <- taxon.tab[j+3,]
              tab.five [i,] <- taxon.tab[j+4,]
              i <- i + 1  
              j <- j + 5
            }
            row.names(tab.one) <- NULL
            row.names(tab.two) <- NULL
            row.names(tab.three) <- NULL
            row.names(tab.four) <- NULL
            row.names(tab.five) <- NULL
            
            tab.one <- na.omit(tab.one)
            tab.two <- na.omit(tab.two)
            tab.three <- na.omit(tab.three)
            tab.four <- na.omit(tab.four)
            tab.five <- na.omit(tab.five)
          }
          
          output$Ctaxa_display_dend = renderUI({
            
            box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
                
                fluidRow(width = 12, align = "center",
                         div(style = "display: inline-block:vertical-align:top;", grVizOutput("dendrogram", height = 1000, width = 1000)) ),
                br(),
                fluidRow(width = 12, align = "center",
                         tagList(
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("bin1_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("bin2_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("bin3_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("bin4_taxonlist") ),
                           div(style="display: inline-block;vertical-align:top;", htmlOutput("bin5_taxonlist") )
                         )
                )
            )
          })
          
          output$dendrogram = renderGrViz({
            flow.text
          })
          output$bin1_taxonlist <- renderText({
            sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab1
          })
          output$bin2_taxonlist <- renderText({
            
            sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab2
          })
          output$bin3_taxonlist <- renderText({
            
            sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab3
          })
          output$bin4_taxonlist <- renderText({
            
            sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab4
          })
          output$bin5_taxonlistt <- renderText({
            
            sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
              kable_styling(latex_options = c('hold_position'))
            sig.tab5
          })
          
          
        } else if (input$chooseMethod_taxa == "Linear regression" | input$chooseMethod_taxa == "Logistic regression" | input$chooseMethod_taxa == "Negative binomial regression" | input$chooseMethod_taxa == "Beta regression") {
          if (input$covariates_taxa == "None") {
            if (input$chooseMethod_taxa =="Linear regression") {
              incProgress(5/10, message = "Linear regression")
              taxa.lm.bin.out <- taxa.bin.lm.united.func(bin.var = taxa.results$bin.var, taxa = taxa.results$taxa)
              taxa.lm.bin.q.out <- bin.q.united.func(taxa.lm.bin.out, method = "BH")
              
              taxa.outputs$DAoutput = taxa.lm.bin.q.out
              nrow <- taxa.forest.plot.pages(taxa.lm.bin.q.out, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Logistic regression") {
              incProgress(5/10, message = "Logistic regression")
              taxa.logit.reg.coef.out <- all.taxa.logit.reg.coef.bin.func(taxa.results$bin.var, taxa.results$taxa, scale = TRUE)
              taxa.logit.reg.coef.q.out <- bin.q.united.func(taxa.logit.reg.coef.out, method = "BH")
              taxa.logit.out <- all.taxa.logit.bin.func(taxa.results$bin.var, taxa.results$taxa, scale = TRUE)
              taxa.logit.q.out <- bin.q.united.func(taxa.logit.out, method = "BH")
              
              taxa.outputs$DAoutput_or <- taxa.logit.q.out
              taxa.outputs$DAoutput <- taxa.logit.reg.coef.q.out
              
              nrow1 <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
              nrow2 <- taxa.forest.plot.pages(taxa.outputs$DAoutput_or, species.include = include)
              
              forestplot.or.data <- taxa.forest.plot.pages1(taxa.outputs$DAoutput_or, chooseData$taxa.names.out, report.type = "OR", species.include = include)
              
            } else if (input$chooseMethod_taxa == "Negative binomial regression") {
              incProgress(5/10, message = "Negative binomial regression")
              taxa.bin.glm.nb.q.out <- all.taxa.bin.glm.nb(taxa.results$bin.var, taxa.results$taxa, taxa.results$lib.size)
              
              taxa.outputs$DAoutput <- taxa.bin.glm.nb.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Beta regression") {
              incProgress(5/10, message = "Beta regression")
              taxa.bin.beta.q.out <- all.taxa.bin.beta(taxa.results$bin.var, taxa.results$taxa)
              
              taxa.outputs$DAoutput <- taxa.bin.beta.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            }
            
          } else if (input$covariates_taxa == "Covariate(s)") {
            if (input$chooseMethod_taxa =="Linear regression") {
              incProgress(5/10, message = "Linear regression")
              
              taxa.lm.bin.cov.out <- taxa.bin.cov.lm.united.func(bin.var = taxa.results$bin.var, 
                                                                 cov.var = taxa.results$cov.var,
                                                                 taxa = taxa.results$taxa)
              
              taxa.lm.bin.cov.q.out <- bin.q.united.func(taxa.lm.bin.cov.out, method = "BH")
              
              taxa.outputs$DAoutput = taxa.lm.bin.cov.q.out
              nrow <- taxa.forest.plot.pages(taxa.lm.bin.cov.q.out, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Logistic regression") {
              incProgress(5/10, message = "Logistic regression")
              taxa.logit.reg.coef.cov.out <- all.taxa.logit.reg.coef.bin.cov.func(taxa.results$bin.var, taxa.results$cov.var, taxa.results$taxa, scale = TRUE)
              taxa.logit.reg.coef.cov.q.out <- bin.q.united.func(taxa.logit.reg.coef.cov.out, method = "BH")
              
              taxa.logit.cov.out <- all.taxa.logit.bin.cov.func(taxa.results$bin.var, taxa.results$cov.var, taxa.results$taxa, scale = TRUE)
              taxa.logit.cov.q.out <- bin.q.united.func(taxa.logit.cov.out, method = "BH")
              
              taxa.outputs$DAoutput_or <- taxa.logit.cov.q.out
              taxa.outputs$DAoutput <- taxa.logit.reg.coef.cov.q.out
              
              nrow1 <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
              nrow2 <- taxa.forest.plot.pages(taxa.outputs$DAoutput_or, species.include = include)
              
              forestplot.or.data <- taxa.forest.plot.pages1(taxa.outputs$DAoutput_or, chooseData$taxa.names.out, report.type = "OR", species.include = include)
              
            } else if (input$chooseMethod_taxa == "Negative binomial regression") {
              incProgress(5/10, message = "Negative binomial regression")
              taxa.bin.cov.glm.nb.q.out <- all.taxa.bin.cov.glm.nb(taxa.results$bin.var, taxa.results$cov.var, taxa.results$taxa, taxa.results$lib.size)
              
              taxa.outputs$DAoutput <- taxa.bin.cov.glm.nb.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
              
            } else if (input$chooseMethod_taxa == "Beta regression") {
              incProgress(5/10, message = "Beta regression")
              taxa.bin.cov.beta.q.out <- all.taxa.bin.cov.beta(taxa.results$bin.var, taxa.results$cov.var, taxa.results$taxa)
              
              taxa.outputs$DAoutput <- taxa.bin.cov.beta.q.out
              nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            }
          }
          
          forestplot.data <- taxa.forest.plot.pages1(taxa.outputs$DAoutput, chooseData$taxa.names.out, report.type = "Est", species.include = include)
          if (any(!is.na(unlist(chooseData$taxa.names.out$duplicates)))) {
            duplicate.taxa <- sapply(strsplit(unlist(chooseData$taxa.names.out$duplicates), " :"),  "[", 1)
            taxon.inplot <- unlist(lapply(forestplot.data$all.text.tab, `[`, i =, j = 3))
            duplicate.texts <- sum(duplicate.taxa %in% taxon.inplot)
          } else {
            duplicate.texts <- 0
          }
          
          if (input$chooseMethod_taxa == "Logistic regression") {
            if (duplicate.texts >0) {
              output$taxa_display = renderUI({
                tagList(
                  tabBox(title = strong("Forest plot", style = "color:black"), width = NULL,
                         tabPanel("Coefficient", align = "center",
                                  do.call(tabsetPanel, lapply(1:nrow1, function(i) {
                                    tabPanel(title = paste0("Page ", i),
                                             plotOutput(paste0("forest", i), height = 800, width = 750),
                                             plotOutput(paste0("duplicates", i), height = 11.5*duplicate.texts+10, width = 750))
                                  }))
                                  
                         )
                         ,
                         tabPanel("Odds Ratio", align = "center",
                                  do.call(tabsetPanel, lapply(1:nrow2, function(i) {
                                    tabPanel(title = paste0("Page ", i),
                                             plotOutput(paste0("forest_or", i), height = 800, width = 750),
                                             plotOutput(paste0("duplicates_or", i), height = 11.5*duplicate.texts+10, width = 750))
                                  }))
                                  
                         )
                  )
                )
              })
            } else {
              output$taxa_display = renderUI({
                tagList(
                  tabBox(title = strong("Forest plot", style = "color:black"), width = NULL,
                         tabPanel("Coefficient", align = "center",
                                  do.call(tabsetPanel, lapply(1:nrow1, function(i) {
                                    tabPanel(title = paste0("Page ", i),
                                             plotOutput(paste0("forest", i), height = 800, width = 750))
                                  }))
                                  
                         )
                         ,
                         tabPanel("Odds Ratio", align = "center",
                                  do.call(tabsetPanel, lapply(1:nrow2, function(i) {
                                    tabPanel(title = paste0("Page ", i),
                                             plotOutput(paste0("forest_or", i), height = 800, width = 750))
                                  })) 
                         )
                  )
                )
              })
            }
            
            lapply(1:nrow1, function(j) {
              output[[paste0("forest", j)]] <- renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data, page = j)
              })
            })
            
            lapply(1:nrow2, function(k) {
              output[[paste0("forest_or", k)]] <- renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.or.data, page = k)
              })
            })
            
            lapply(1:nrow1, function(j) {
              output[[paste0("duplicates", j)]] <- renderPlot({
                duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
              })
            })
            
            lapply(1:nrow2, function(j) {
              output[[paste0("duplicates_or", j)]] <- renderPlot({
                duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
              })
            })
            
            output$duplicates = renderPlot({
              duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
            })
            
            output$duplicates_or = renderPlot({
              duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
            })
            
            taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)
            
            flow.text <- taxa.sig$flow.text
            taxon.tab <- taxa.sig$taxon.tab
            ci.tab.all <- taxa.sig$ci.tab.all
            
            if ( include == FALSE){
              ci.tab.all <- ci.tab.all[-(grep("s_", taxon.tab$Taxon, fixed = TRUE)+1)]
              taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
            }
            
            if ( length(ci.tab.all) > 1 ){
              
              for( i in 1:nrow(taxon.tab)){
                if ( ci.tab.all[-1][i] < 0){
                  taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
                  taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
                }
                else{
                  taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
                  taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
                }
                
              }
            }
            
            
            N <- dim(taxon.tab)[1]
            itr <- ceiling(N/5)
            tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
            tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
            tab.three <- data.frame( matrix(ncol=2,nrow=0) )
            tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
            tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
            
            colnames(tab.one) <- c("ID", "Taxon")
            colnames(tab.two) <- c("ID", "Taxon")
            colnames(tab.three) <- c("ID", "Taxon")
            colnames(tab.four) <- c("ID", "Taxon")
            colnames(tab.five) <- c("ID", "Taxon")
            
            if ( dim(taxon.tab)[1] > 0 ) {
              
              i = 1
              j = 1
              N.rmd <- N %% 5
              N.fix <- N + 5 - N.rmd
              while ( i <= itr) {
                tab.one  [i,] <- taxon.tab[j,]
                tab.two  [i,] <- taxon.tab[j+1,]
                tab.three[i,] <- taxon.tab[j+2,]
                tab.four [i,] <- taxon.tab[j+3,]
                tab.five [i,] <- taxon.tab[j+4,]
                i <- i + 1  
                j <- j + 5
              }
              row.names(tab.one) <- NULL
              row.names(tab.two) <- NULL
              row.names(tab.three) <- NULL
              row.names(tab.four) <- NULL
              row.names(tab.five) <- NULL
              
              tab.one <- na.omit(tab.one)
              tab.two <- na.omit(tab.two)
              tab.three <- na.omit(tab.three)
              tab.four <- na.omit(tab.four)
              tab.five <- na.omit(tab.five)
            }
            
            output$Ctaxa_display_dend = renderUI({
              
              box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
                  
                  fluidRow(width = 12, align = "center",
                           div(style = "display: inline-block:vertical-align:top;", grVizOutput("dendrogram", height = 1000, width = 1000)) ),
                  br(),
                  fluidRow(width = 12, align = "center",
                           tagList(
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin1_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin2_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin3_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin4_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin5_taxonlist") )
                           )
                  )
              )
            })
            
            output$dendrogram = renderGrViz({
              flow.text
            })
            output$bin1_taxonlist <- renderText({
              sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab1
            })
            output$bin2_taxonlist <- renderText({
              
              sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab2
            })
            output$bin3_taxonlist <- renderText({
              
              sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab3
            })
            output$bin4_taxonlist <- renderText({
              
              sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab4
            })
            output$bin5_taxonlistt <- renderText({
              
              sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab5
            })
            
            output$duplicates = renderPlot({
              duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
            })
            
          } 
          else {
            if (duplicate.texts>0) {
              output$taxa_display = renderUI({
                tagList(
                  do.call(tabsetPanel, lapply(1:nrow, function(i) {
                    tabPanel(title = paste0("Page ", i), align = "center",
                             plotOutput(paste0("forest", i), height = 800, width = 750),
                             plotOutput(paste0("duplicates", i), height = 11.5*duplicate.texts+10, width = 750))
                  })) 
                )
              })
            } else {
              output$taxa_display = renderUI({
                tagList(
                  do.call(tabsetPanel, lapply(1:nrow, function(i) {
                    tabPanel(title = paste0("Page ", i), align = "center",
                             plotOutput(paste0("forest", i), height = 800, width = 750))
                  })) 
                )
              })
            }
            
            
            lapply(1:nrow, function(j) {
              output[[paste0("forest", j)]] <- renderPlot({
                taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data, page = j)
              })
            })
            
            lapply(1:nrow, function(j) {
              output[[paste0("duplicates", j)]] <- renderPlot({
                duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
              })
            })
            
            output$duplicates = renderPlot({
              duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
            })
            
            taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)
            
            flow.text <- taxa.sig$flow.text
            taxon.tab <- taxa.sig$taxon.tab
            ci.tab.all <- taxa.sig$ci.tab.all
            
            if ( include == FALSE){
              ci.tab.all <- ci.tab.all[-(grep("s_", taxon.tab$Taxon, fixed = TRUE)+1)]
              taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
            }
            
            if ( length(ci.tab.all) > 1 ){
              
              for( i in 1:nrow(taxon.tab)){
                if ( ci.tab.all[-1][i] < 0){
                  taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
                  taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
                }
                else{
                  taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
                  taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
                }
                
              }
            }
            
            
            N <- dim(taxon.tab)[1]
            itr <- ceiling(N/5)
            tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
            tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
            tab.three <- data.frame( matrix(ncol=2,nrow=0) )
            tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
            tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
            
            colnames(tab.one) <- c("ID", "Taxon")
            colnames(tab.two) <- c("ID", "Taxon")
            colnames(tab.three) <- c("ID", "Taxon")
            colnames(tab.four) <- c("ID", "Taxon")
            colnames(tab.five) <- c("ID", "Taxon")
            
            if ( dim(taxon.tab)[1] > 0 ) {
              
              i = 1
              j = 1
              N.rmd <- N %% 5
              N.fix <- N + 5 - N.rmd
              while ( i <= itr) {
                tab.one  [i,] <- taxon.tab[j,]
                tab.two  [i,] <- taxon.tab[j+1,]
                tab.three[i,] <- taxon.tab[j+2,]
                tab.four [i,] <- taxon.tab[j+3,]
                tab.five [i,] <- taxon.tab[j+4,]
                i <- i + 1  
                j <- j + 5
              }
              row.names(tab.one) <- NULL
              row.names(tab.two) <- NULL
              row.names(tab.three) <- NULL
              row.names(tab.four) <- NULL
              row.names(tab.five) <- NULL
              
              tab.one <- na.omit(tab.one)
              tab.two <- na.omit(tab.two)
              tab.three <- na.omit(tab.three)
              tab.four <- na.omit(tab.four)
              tab.five <- na.omit(tab.five)
            }
            
            output$Ctaxa_display_dend = renderUI({
              
              box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
                  
                  fluidRow(width = 12, align = "center",
                           div(style = "display: inline-block:vertical-align:top;", grVizOutput("dendrogram", height = 1000, width = 1000)) ),
                  br(),
                  fluidRow(width = 12, align = "center",
                           tagList(
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin1_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin2_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin3_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin4_taxonlist") ),
                             div(style="display: inline-block;vertical-align:top;", htmlOutput("bin5_taxonlist") )
                           )
                  )
              )
            })
            
            output$dendrogram = renderGrViz({
              flow.text
            })
            output$bin1_taxonlist <- renderText({
              sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab1
            })
            output$bin2_taxonlist <- renderText({
              
              sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab2
            })
            output$bin3_taxonlist <- renderText({
              
              sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab3
            })
            output$bin4_taxonlist <- renderText({
              
              sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab4
            })
            output$bin5_taxonlistt <- renderText({
              
              sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
                kable_styling(latex_options = c('hold_position'))
              sig.tab5
            })
            
            output$duplicates = renderPlot({
              duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
            })
          } 
        }
        
        incProgress(1/10, message = "Displaying Results in progress")
        output$downloadTable_taxa = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"),
                h5("Summary Statistics"),
                downloadButton("tdownloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("tdownloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3"),
                h5("Dendrogram"),
                downloadButton("gdownload", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        output$tdownloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Taxa.Sum.Table.zip")
          },
          content = function(sum.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=sum.file, files=dataFiles)
          }
        )
        
        output$tdownloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Output.zip")
          },
          content = function(DA.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.outputs$DAoutput$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=DA.file, files=dataFiles)
          }
        )
        
        output$gdownload <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
          },
          content = function(file) {
            htmlwidgets::saveWidget(as_widget(taxa.sig.dend(taxa.outputs$DAoutput, chooseData$tax.tab, "twopi", include)), file)
          }
        )
        
        ref_string = REFERENCE_CHECK(data_transform = input$dataType_taxa, method_name = isolate(input$chooseMethod_taxa), FDR = "Yes")
        #ref_string = REFERENCE_CHECK(method_name = isolate(input$chooseMethod), FDR = isolate(input$chooseAdjustment))
        if (is.null(ref_string)) {
          shinyjs::hide("taxa_references")
        } else {
          shinyjs::show("taxa_references")
          output$taxa_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        delay(1000, shinyjs::enable("runbtn_bin_taxa"))
        delay(1000, shinyjs::enable("dataType_taxa"))
        delay(1000, shinyjs::enable("primvar_taxa"))
        delay(1000, shinyjs::enable("chooseMethod_taxa"))
        delay(1000, shinyjs::enable("covariates_taxa"))
        delay(1000, shinyjs::enable("chooseRef_taxa"))
        delay(1000, shinyjs::enable("taxaCat1"))
        delay(1000, shinyjs::enable("taxaCat2"))
        
      })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  observeEvent(input$runbtn_cont_taxa,{
    validate(
      if (input$covariates_taxa == "Covariate(s)" & is.null(input$covariatesOptions_taxa)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    shinyjs::disable("runbtn_cont_taxa")
    shinyjs::disable("dataType_taxa")
    shinyjs::disable("primvar_taxa")
    shinyjs::disable("chooseMethod_taxa")
    shinyjs::disable("covariates_taxa")
    shinyjs::disable("rename.con.var_taxa")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        if (input$chooseMethod_taxa == "Negative binomial regression") {
          taxa.types$dataType = "count"
        }
        
        incProgress(1/10, message = "Examining Data in progress")
        if (input$covariates_taxa == "None") {
          taxa.con.out <- taxa.con.recode.func(chooseData$sam.dat, input$primvar_taxa, input$rename.con.var_taxa, 
                                               taxa = chooseData$taxa.out[[taxa.types$dataType]])
          taxa.results$con.var <- taxa.con.out$con.var
          taxa.results$taxa <- taxa.con.out$taxa
          taxa.results$taxa.con.sum.out <- taxa.sum.apply(taxa.con.out,2,taxa.ind.sum.func)
        } else if (input$covariates_taxa == "Covariate(s)") {
          taxa.con.cov.out <- taxa.con.cov.recode.func(chooseData$sam.dat, input$primvar_taxa, input$covariatesOptions_taxa, 
                                                       input$rename.con.var_taxa, chooseData$taxa.out[[taxa.types$dataType]])  
          
          taxa.results$con.var <- taxa.con.cov.out$con.var
          taxa.results$cov.var <- taxa.con.cov.out$cov.var
          taxa.results$taxa <- taxa.con.cov.out$taxa
          taxa.results$taxa.con.sum.out <- taxa.sum.apply(taxa.con.cov.out,2,taxa.ind.sum.func)
        }
        
        taxa_dataConvar <- taxa.results$con.var
        taxa_dataTaxa <- taxa.results$taxa
        taxa.results$lib.size <- taxa.results$lib.size[names(taxa.results$lib.size) %in% rownames(chooseData$sam.dat)]
        
        if (input$include_species.dend == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        #incProgress(5/10, message = "Fitting Regression Model in progress")
        if (input$chooseMethod_taxa == "Linear regression") {
          incProgress(5/10, message = "Linear regression")
          if (input$covariates_taxa == "Covariate(s)") {
            
            taxa.lm.con.cov.out <- taxa.lm.con.cov.func(taxa.results$con.var, taxa.results$cov.var, taxa.results$taxa)
            taxa.lm.con.cov.q.out <- taxa.q.func(taxa.lm.con.cov.out, method = "BH")
            
            taxa.outputs$DAoutput <- taxa.lm.con.cov.q.out
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            
          } else if (input$covariates_taxa == "None") {
            
            taxa.lm.con.out <- taxa.lm.con.func(taxa.results$con.var, taxa.results$taxa)
            taxa.lm.con.q.out <- taxa.q.func(taxa.lm.con.out, method = "BH")
            
            taxa.outputs$DAoutput <- taxa.lm.con.q.out
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            
          }
        } else if (input$chooseMethod_taxa == "Negative binomial regression") {
          incProgress(5/10, message = "Negative binomial regression")
          if (input$covariates_taxa == "Covariate(s)") {
            taxa.con.cov.glm.nb.out <- all.taxa.con.cov.glm.nb(taxa.results$con.var, taxa.results$cov.var, taxa.results$taxa, taxa.results$lib.size, "BH") 
            
            taxa.outputs$DAoutput <- taxa.con.cov.glm.nb.out
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            
          } else if (input$covariates_taxa == "None") {
            taxa.con.glm.nb.out <- all.taxa.con.glm.nb(taxa.results$con.var, taxa.results$taxa, taxa.results$lib.size, "BH") 
            
            taxa.outputs$DAoutput <- taxa.con.glm.nb.out
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
          }
          
        } else if (input$chooseMethod_taxa == "Beta regression") {
          incProgress(5/10, message = "Beta regression")
          if (input$covariates_taxa == "Covariate(s)") {
            taxa.con.cov.beta.q.out <- all.taxa.con.cov.beta(taxa.results$con.var, taxa.results$cov.var, taxa.results$taxa, multi.method = "BH")
            
            taxa.outputs$DAoutput <- taxa.con.cov.beta.q.out
            
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
            
          } else if (input$covariates_taxa == "None") {
            taxa.con.beta.q.out <- all.taxa.con.beta(taxa.results$con.var, taxa.results$taxa, multi.method = "BH")
            
            taxa.outputs$DAoutput <- taxa.con.beta.q.out
            
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutput, species.include = include)
          }
        }
        
        forestplot.data <- taxa.forest.plot.pages1(taxa.outputs$DAoutput, chooseData$taxa.names.out, species.include = include)
        #duplicate.texts <- sum(!is.na(unlist(chooseData$taxa.names.out$duplicates)))
        
        if (any(!is.na(unlist(chooseData$taxa.names.out$duplicates)))) {
          duplicate.taxa <- sapply(strsplit(unlist(chooseData$taxa.names.out$duplicates), " :"),  "[", 1)
          #sum(!is.na(unlist(chooseData$taxa.names.out$duplicates)))
          taxon.inplot <- unlist(lapply(forestplot.data$all.text.tab, `[`, i =, j = 3))
          duplicate.texts <- sum(taxon.inplot %in% duplicate.taxa)
        } else {
          duplicate.texts <- 0
        }
        
        incProgress(2/10, message = "Displaying Results in progress")
        
        if (duplicate.texts>0) {
          output$taxa_display = renderUI({
            do.call(tabsetPanel, lapply(1:nrow, function(i) {
              tabPanel(title = paste0("Page ", i), align = "center",
                       plotOutput(paste0("forest", i), height = 800, width = 750),
                       plotOutput(paste0("duplicates", i), height = 11.5*duplicate.texts+10, width = 750))
            }))
            
          })
        } else {
          output$taxa_display = renderUI({
            do.call(tabsetPanel, lapply(1:nrow, function(i) {
              tabPanel(title = paste0("Page ", i), align = "center",
                       plotOutput(paste0("forest", i), height = 800, width = 750))
            }))
            
          })
        }
        
        
        lapply(1:nrow, function(j) {
          output[[paste0("forest", j)]] <- renderPlot({
            taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data, page = j)
          })
        })
        
        lapply(1:nrow, function(j) {
          output[[paste0("duplicates", j)]] <- renderPlot({
            duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
          })
        })
        
        taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutput, chooseData$NAadded$tax.tab, "twopi", include)
        
        flow.text <- taxa.sig$flow.text
        taxon.tab <- taxa.sig$taxon.tab
        ci.tab.all <- taxa.sig$ci.tab.all
        
        if ( include == FALSE){
          ci.tab.all <- ci.tab.all[-(grep("s_", taxon.tab$Taxon, fixed = TRUE)+1)]
          taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
        }
        
        if ( length(ci.tab.all) > 1 ){
          
          for( i in 1:nrow(taxon.tab)){
            if ( ci.tab.all[-1][i] < 0){
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
            else{
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
            
          }
        }
        
        N <- dim(taxon.tab)[1]
        itr <- ceiling(N/5)
        tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.three <- data.frame( matrix(ncol=2,nrow=0) )
        tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
        tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
        
        colnames(tab.one) <- c("ID", "Taxon")
        colnames(tab.two) <- c("ID", "Taxon")
        colnames(tab.three) <- c("ID", "Taxon")
        colnames(tab.four) <- c("ID", "Taxon")
        colnames(tab.five) <- c("ID", "Taxon")
        
        if ( dim(taxon.tab)[1] > 0 ) {
          
          i = 1
          j = 1
          N.rmd <- N %% 5
          N.fix <- N + 5 - N.rmd
          while ( i <= itr) {
            tab.one  [i,] <- taxon.tab[j,]
            tab.two  [i,] <- taxon.tab[j+1,]
            tab.three[i,] <- taxon.tab[j+2,]
            tab.four [i,] <- taxon.tab[j+3,]
            tab.five [i,] <- taxon.tab[j+4,]
            i <- i + 1  
            j <- j + 5
          }
          row.names(tab.one) <- NULL
          row.names(tab.two) <- NULL
          row.names(tab.three) <- NULL
          row.names(tab.four) <- NULL
          row.names(tab.five) <- NULL
          
          tab.one <- na.omit(tab.one)
          tab.two <- na.omit(tab.two)
          tab.three <- na.omit(tab.three)
          tab.four <- na.omit(tab.four)
          tab.five <- na.omit(tab.five)
        }
        
        output$Ctaxa_display_dend = renderUI({
          
          box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
              
              fluidRow(width = 12, align = "center",
                       div(style = "display: inline-block:vertical-align:top;", grVizOutput("dendrogram", height = 1000, width = 1000)) ),
              br(),
              fluidRow(width = 12, align = "center",
                       tagList(
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("con1_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("con2_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("con3_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("con4_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("con5_taxonlist") )
                       )
              )
          )
        })
        
        output$dendrogram = renderGrViz({
          flow.text
        })
        output$con1_taxonlist <- renderText({
          sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab1
        })
        output$con2_taxonlist <- renderText({
          
          sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab2
        })
        output$con3_taxonlist <- renderText({
          
          sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab3
        })
        output$con4_taxonlist <- renderText({
          
          sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab4
        })
        output$con5_taxonlistt <- renderText({
          
          sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab5
        })
        
        incProgress(1/10, message = "Displaying Results in progress")
        output$downloadTable_taxa = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"),
                h5("Summary Statistics"),
                downloadButton("tdownloadTabl1", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("tdownloadTabl2", "Download", width = '50%', style = "color:black; background-color: red3"),
                h5("Dendrogram"),
                downloadButton("gdownload", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        incProgress(1/10, message = "Displaying Results in progress")
        
        output$tdownloadTabl1 <- downloadHandler(
          filename = function() {
            paste("Taxa.Sum.Table.zip")
          },
          content = function(sum.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=sum.file, files=dataFiles)
          }
        )
        
        output$tdownloadTabl2 <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Output.zip")
          },
          content = function(DA.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.outputs$DAoutput$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutput$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=DA.file, files=dataFiles)
          }
        )
        
        output$gdownload <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
          },
          content = function(file) {
            htmlwidgets::saveWidget(as_widget(taxa.sig.dend(taxa.outputs$DAoutput, chooseData$tax.tab, "twopi", include)), file)
          }
        )
        
        ref_string = REFERENCE_CHECK(data_transform = input$dataType_taxa, method_name = isolate(input$chooseMethod_taxa), FDR = "Yes")
        
        if (is.null(ref_string)) {
          shinyjs::hide("taxa_references")
        } else {
          shinyjs::show("taxa_references")
          output$taxa_references = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        delay(1000, shinyjs::enable("runbtn_cont_taxa"))
        delay(1000, shinyjs::enable("dataType_taxa"))
        delay(1000, shinyjs::enable("primvar_taxa"))
        delay(1000, shinyjs::enable("chooseMethod_taxa"))
        delay(1000, shinyjs::enable("covariates_taxa"))
        delay(1000, shinyjs::enable("rename.con.var_taxa"))
      }
    )
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  ######################################
  ## Taxa Longitudinal Data Analysis ###
  ######################################
  observeEvent(input$runbtn_bin_taxa.long,{
    validate(
      if (input$covariates_taxa.long == "Covariate(s)" & is.null(input$covariatesOptions_taxa.long)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    validate(
      if (is.null(input$clustervar_taxa)) {
        showNotification("Error: No cluster variable available.",
                         type = "error")
      } else {
        NULL
      }
    )
    
    shinyjs::disable("runbtn_bin_taxa.long")
    shinyjs::disable("dataType_taxa.long")
    shinyjs::disable("primvar_taxa.long")
    shinyjs::disable("chooseMethod_taxa.long")
    shinyjs::disable("covariates_taxa.long")
    #shinyjs::disable("chooseRef_taxa")
    shinyjs::disable("taxaCat1")
    shinyjs::disable("taxaCat2")
    shinyjs::disable("clustervar_taxa")
    shinyjs::disable("include_species.dend.long")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        taxa.categors = c(taxa.categos$cat1, taxa.categos$cat2)
        taxa.bin_categos = c(input$taxaCat1, input$taxaCat2)
        
        rename.cats_ref = taxa.bin_categos[which(taxa.categors == taxa.categos$cat1)]
        rename.cats_com = taxa.bin_categos[which(taxa.categors != taxa.categos$cat1)]
        
        taxa.bin.cat.ref.ori.out <- taxa.bin.cat.ref.ori.func(chooseData$sam.dat, input$primvar_taxa.long)          # sel.bin.var = input$primvar_taxa = "ecig_status"
        
        incProgress(1/10, message = "Examining Data in progress")
        sam_dat <- taxa.bin.cat.recode.func(chooseData$sam.dat, input$primvar_taxa.long, taxa.bin.cat.ref.ori.out,
                                            rename.ref = rename.cats_ref, rename.com = rename.cats_com)
        
        if (input$chooseMethod_taxa.long == "GLMM (Negative Binomial)") {
          taxa.types$dataType = "count"
        }
        
        if (input$covariates_taxa.long == "None") {
          taxa.bin.id.out <- taxa.bin.id.cat.ref.united.func(input$primvar_taxa.long, input$clustervar_taxa, rename.cats_ref,
                                                             rename.cats_com, sam_dat, taxa = chooseData$taxa.out[[taxa.types$dataType]])
          taxa.results$bin.var <- taxa.bin.id.out$bin.var
          taxa.results$id.var <- taxa.bin.id.out$id.var
          taxa.results$taxa <- taxa.bin.id.out$taxa
          
        } else if (input$covariates_taxa.long == "Covariate(s)") {
          taxa.bin.id.cov.out <- taxa.bin.id.cov.cat.ref.united.func(input$primvar_taxa.long, input$clustervar_taxa, input$covariatesOptions_taxa.long,
                                                                     rename.cats_ref, rename.cats_com, sam_dat, chooseData$taxa.out[[taxa.types$dataType]])
          taxa.results$bin.var <- taxa.bin.id.cov.out$bin.var
          taxa.results$id.var <- taxa.bin.id.cov.out$id.var
          taxa.results$taxa <- taxa.bin.id.cov.out$taxa
          taxa.results$cov.var <- taxa.bin.id.cov.out$cov.var
          
        }
        
        # check if this is equivalent to the one for longitudinal
        taxa.results$taxa.bin.sum.out <- taxa.bin.id.sum.united.func(taxa.results$bin.var, taxa.results$taxa, sel.ref = rename.cats_ref, sel.com = rename.cats_com)
        
        taxa_dataBinvar <- taxa.results$bin.var
        taxa_dataTaxa <- taxa.results$taxa
        
        if (input$include_species.dend.long == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        #incProgress(5/10, message = "Fitting Regression Model")
        if (input$chooseMethod_taxa.long =="LMM") {
          incProgress(5/10, message = "LMM")
          if (input$covariates_taxa.long == "None") {
            taxa.bin.lmer.q.out <- all.taxa.bin.lmer(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, multi.method = "BH") 
            taxa.outputs$DAoutputlong = taxa.bin.lmer.q.out
            
          } else if (input$covariates_taxa.long == "Covariate(s)") {
            taxa.bin.cov.lmer.q.out <- all.taxa.bin.cov.lmer(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, multi.method = "BH")
            taxa.outputs$DAoutputlong = taxa.bin.cov.lmer.q.out
          }
          nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
          
        } else if (input$chooseMethod_taxa.long == "GLMM (Beta)") {
          incProgress(5/10, message = "GLMM (Beta)")
          if (input$covariates_taxa.long == "None") {
            taxa.bin.beta.q.out <- all.taxa.bin.glmm.beta(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.beta.q.out
            
          } else if (input$covariates_taxa.long == "Covariate(s)") {
            
            taxa.bin.cov.beta.q.out <- all.taxa.bin.glmm.cov.beta(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.cov.beta.q.out
          }
          nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
        } else if (input$chooseMethod_taxa.long == "GLMM (Negative Binomial)") {
          incProgress(5/10, message = "GLMM (Negative Binomial)")
          if (input$covariates_taxa.long == "None") {
            taxa.bin.glmm.nb.q.out <- all.taxa.bin.glmm.nb(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, taxa.results$lib.size, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.glmm.nb.q.out
            
          } else if (input$covariates_taxa.long == "Covariate(s)") {
            
            taxa.bin.cov.glmm.nb.q.out <- all.taxa.bin.cov.glmm.nb(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, taxa.results$lib.size, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.cov.glmm.nb.q.out
          }
          nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
          
        } else if (input$chooseMethod_taxa.long == "GEE (Binomial)") {
          incProgress(5/10, message = "GEE (Binomial)")
          if (input$covariates_taxa.long == "None") {
            
            taxa.bin.logit.coef.gee.out <- all.taxa.bin.logit.gee.reg.coef(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, rare.count = ("rare.count"==taxa.types$dataType), multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.logit.coef.gee.out
            
          } else if (input$covariates_taxa.long == "Covariate(s)") {
            
            taxa.bin.cov.logit.coef.gee.out <- all.taxa.bin.cov.logit.gee.reg.coef(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, rare.count = ("rare.count"==taxa.types$dataType), multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.cov.logit.coef.gee.out
          }
          
          nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
          
        } else if (input$chooseMethod_taxa.long == "GLMM (Binomial)") {
          incProgress(5/10, message = "GLMM (Binomial)")
          if (input$covariates_taxa.long == "None") {
            
            taxa.bin.logit.glmm.b.out <- all.taxa.bin.logit.glmm.b(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, rare.count = ("rare.count"==taxa.types$dataType), multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.logit.glmm.b.out
            
          } else if (input$covariates_taxa.long == "Covariate(s)") {
            
            taxa.bin.cov.logit.glmm.b.out <- all.taxa.bin.cov.logit.glmm.b(bin.var = taxa.results$bin.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, rare.count = ("rare.count"==taxa.types$dataType), multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.bin.cov.logit.glmm.b.out
            
          }
          nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
        }
        
        forestplot.data <- taxa.forest.plot.pages1(taxa.outputs$DAoutputlong, chooseData$taxa.names.out, species.include = include)
        
        if (any(!is.na(unlist(chooseData$taxa.names.out$duplicates)))) {
          duplicate.taxa <- sapply(strsplit(unlist(chooseData$taxa.names.out$duplicates), " :"),  "[", 1)
          #sum(!is.na(unlist(chooseData$taxa.names.out$duplicates)))
          taxon.inplot <- unlist(lapply(forestplot.data$all.text.tab, `[`, i =, j = 3))
          duplicate.texts <- sum(taxon.inplot %in% duplicate.taxa)
        } else {
          duplicate.texts <- 0
        }
        
        if (duplicate.texts >0) {
          output$taxa_displaylong = renderUI({
            do.call(tabsetPanel, lapply(1:nrow, function(i) {
              tabPanel(title = paste0("Page ", i), align = "center",
                       plotOutput(paste0("forestlong", i), height = 800, width = 750),
                       plotOutput(paste0("duplicateslong", i), height = 11.5*duplicate.texts+10, width = 750))
            }))
          })
        } else {
          output$taxa_displaylong = renderUI({
            do.call(tabsetPanel, lapply(1:nrow, function(i) {
              tabPanel(title = paste0("Page ", i), align = "center",
                       plotOutput(paste0("forestlong", i), height = 800, width = 750))
            }))
          })
        }
        
        
        lapply(1:nrow, function(j) {
          output[[paste0("forestlong", j)]] <- renderPlot({
            taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data, page = j)
          })
        })
        
        lapply(1:nrow, function(j) {
          output[[paste0("duplicateslong", j)]] <- renderPlot({
            duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
          })
        })
        
        taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutputlong, chooseData$NAadded$tax.tab, "twopi", include)
        
        flow.text <- taxa.sig$flow.text
        taxon.tab <- taxa.sig$taxon.tab
        ci.tab.all <- taxa.sig$ci.tab.all
        
        if ( include == FALSE){
          ci.tab.all <- ci.tab.all[-(grep("s_", taxon.tab$Taxon, fixed = TRUE)+1)]
          taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
        }
        
        if ( length(ci.tab.all) > 1 ){
          
          for( i in 1:nrow(taxon.tab)){
            if ( ci.tab.all[-1][i] < 0){
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
            else{
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
            
          }
        }
        
        N <- dim(taxon.tab)[1]
        itr <- ceiling(N/5)
        tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.three <- data.frame( matrix(ncol=2,nrow=0) )
        tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
        tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
        
        colnames(tab.one) <- c("ID", "Taxon")
        colnames(tab.two) <- c("ID", "Taxon")
        colnames(tab.three) <- c("ID", "Taxon")
        colnames(tab.four) <- c("ID", "Taxon")
        colnames(tab.five) <- c("ID", "Taxon")
        
        if ( dim(taxon.tab)[1] > 0 ) {
          
          i = 1
          j = 1
          N.rmd <- N %% 5
          N.fix <- N + 5 - N.rmd
          while ( i <= itr) {
            tab.one  [i,] <- taxon.tab[j,]
            tab.two  [i,] <- taxon.tab[j+1,]
            tab.three[i,] <- taxon.tab[j+2,]
            tab.four [i,] <- taxon.tab[j+3,]
            tab.five [i,] <- taxon.tab[j+4,]
            i <- i + 1  
            j <- j + 5
          }
          row.names(tab.one) <- NULL
          row.names(tab.two) <- NULL
          row.names(tab.three) <- NULL
          row.names(tab.four) <- NULL
          row.names(tab.five) <- NULL
          
          tab.one <- na.omit(tab.one)
          tab.two <- na.omit(tab.two)
          tab.three <- na.omit(tab.three)
          tab.four <- na.omit(tab.four)
          tab.five <- na.omit(tab.five)
        }
        
        output$Ltaxa_display_dend = renderUI({
          
          box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
              
              fluidRow(width = 12, align = "center",
                       div(style = "display: inline-block:vertical-align:top;", grVizOutput("Ldendrogram", height = 1000, width = 1000)) ),
              br(),
              fluidRow(width = 12, align = "center",
                       tagList(
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lbin1_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lbin2_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lbin3_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lbin4_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lbin5_taxonlist") )
                       )
              )
          )
        })
        
        output$Ldendrogram = renderGrViz({
          flow.text
        })
        
        output$Lbin1_taxonlist <- renderText({
          sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab1
        })
        output$Lbin2_taxonlist <- renderText({
          sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab2
        })
        output$Lbin3_taxonlist <- renderText({
          
          sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab3
        })
        output$Lbin4_taxonlist <- renderText({
          
          sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab4
        })
        output$Lbin5_taxonlist <- renderText({
          
          sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab5
        })
        
        incProgress(1/10, message = "Displaying Results in progress")
        output$downloadTable_taxalong = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Summary Statistics"),
                downloadButton("tdownloadTabl1long", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("tdownloadTabl2long", "Download", width = '50%', style = "color:black; background-color: red3"),
                h5("Dendrogram"),
                downloadButton("gdownloadlong", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        output$tdownloadTabl1long <- downloadHandler(
          filename = function() {
            paste("Taxa.Sum.Table.zip")
          },
          content = function(sum.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.bin.sum.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=sum.file, files=dataFiles)
          }
        )
        
        output$tdownloadTabl2long <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Output.zip")
          },
          content = function(DA.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=DA.file, files=dataFiles)
          }
        )
        
        output$gdownloadlong <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
          },
          content = function(file) {
            htmlwidgets::saveWidget(as_widget(taxa.sig.dend(taxa.outputs$DAoutputlong, chooseData$tax.tab, "twopi", include)), file)
          }
        )
        
        ref_string = REFERENCE_CHECK(data_transform = input$dataType_taxa.long, method_name = isolate(input$chooseMethod_taxa.long), FDR = "Yes")
        if (is.null(ref_string)) {
          shinyjs::hide("taxa_references_long")
        } else {
          shinyjs::show("taxa_references_long")
          output$taxa_references_long = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        
        shinyjs::enable("runbtn_bin_taxa.long")
        shinyjs::enable("dataType_taxa.long")
        shinyjs::enable("primvar_taxa.long")
        shinyjs::enable("chooseMethod_taxa.long")
        shinyjs::enable("covariates_taxa.long")
        #shinyjs::enable("chooseRef_taxa")
        shinyjs::enable("taxaCat1")
        shinyjs::enable("taxaCat2")
        shinyjs::enable("clustervar_taxa")
        shinyjs::enable("include_species.dend.long")
        
      })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  observeEvent(input$runbtn_cont_taxa.long,{
    validate(
      if (input$covariates_taxa.long == "Covariate(s)" & is.null(input$covariatesOptions_taxa.long)) {
        showNotification("Error: Please select covariate(s) before you click 'Run!' button.",
                         type = "error")
      } else {
        NULL
      }
    )
    shinyjs::disable("runbtn_cont_taxa.long")
    shinyjs::disable("dataType_taxa.long")
    shinyjs::disable("primvar_taxa.long")
    shinyjs::disable("chooseMethod_taxa.long")
    shinyjs::disable("covariates_taxa.long")
    shinyjs::disable("rename.con.var_taxa")
    shinyjs::disable("include_species.dend.long")
    shinyjs::disable("clustervar_taxa")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        
        if (input$chooseMethod_taxa.long == "GLMM (Negative Binomial)") {
          taxa.types$dataType = "count"
        }
        
        incProgress(1/10, message = "Examining Data in progress")
        if (input$covariates_taxa.long == "None") {
          taxa.con.id.out <- taxa.con.id.recode.united.func(chooseData$sam.dat, input$primvar_taxa.long, input$clustervar_taxa, input$rename.con.var_taxa,
                                                            taxa = chooseData$taxa.out[[taxa.types$dataType]])
          
          taxa.results$con.var <- taxa.con.id.out$con.var
          taxa.results$id.var <- taxa.con.id.out$id.var
          taxa.results$taxa <- taxa.con.id.out$taxa
          taxa.results$taxa.con.sum.out <- taxa.sum.apply(taxa.con.id.out,2,taxa.ind.sum.func)
        } else if (input$covariates_taxa.long == "Covariate(s)") {
          taxa.con.cov.id.out <- taxa.con.id.cov.recode.united.func(chooseData$sam.dat, input$primvar_taxa.long, input$rename.con.var_taxa, 
                                                                    input$clustervar_taxa, input$covariatesOptions_taxa.long, chooseData$taxa.out[[taxa.types$dataType]])
          
          taxa.results$con.var <- taxa.con.cov.id.out$con.var
          taxa.results$id.var <- taxa.con.cov.id.out$id.var
          taxa.results$cov.var <- taxa.con.cov.id.out$cov.var
          taxa.results$taxa <- taxa.con.cov.id.out$taxa
          taxa.results$taxa.con.sum.out <- taxa.sum.apply(taxa.con.cov.id.out,2,taxa.ind.sum.func)
        }
        
        taxa_dataConvar <- taxa.results$con.var
        taxa_dataTaxa <- taxa.results$taxa
        
        if (input$include_species.dend.long == "Phylum - Species") {
          include = TRUE
        } else {
          include = FALSE
        }
        
        #incProgress(5/10, message = "Fitting Regression Model in progress")
        if (input$chooseMethod_taxa.long == "LMM") {
          incProgress(5/10, message = "LMM")
          if (input$covariates_taxa.long == "Covariate(s)") {
            taxa.con.cov.lmer.q.out <- all.taxa.con.cov.lmer(con.var = taxa.results$con.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, multi.method = "BH") 
            #taxa.data.results$table.out <- taxa.con.cov.lmer.q.out
            taxa.outputs$DAoutputlong <- taxa.con.cov.lmer.q.out
            
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
            
          } else if (input$covariates_taxa.long == "None") {
            taxa.con.lmer.q.out <- all.taxa.con.lmer(con.var = taxa.results$con.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, multi.method = "BH") 
            #taxa.data.results$table.out <- taxa.con.lmer.q.out
            taxa.outputs$DAoutputlong <- taxa.con.lmer.q.out
            
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
          }
        } else if (input$chooseMethod_taxa.long == "GLMM (Beta)") {
          incProgress(5/10, message = "GLMM (Beta)")
          if (input$covariates_taxa.long == "None") {
            taxa.con.beta.q.out <- all.taxa.con.glmm.beta(con.var = taxa.results$con.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.con.beta.q.out
            
          } else if (input$covariates_taxa.long == "Covariate(s)") {
            taxa.con.cov.beta.q.out <- all.taxa.con.glmm.cov.beta(con.var = taxa.results$con.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.con.cov.beta.q.out
          }
          nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
        } else if (input$chooseMethod_taxa.long == "GLMM (Negative Binomial)") {
          
          incProgress(5/10, message = "GLMM (Negative Binomial)")
          if (input$covariates_taxa.long == "Covariate(s)") {
            taxa.con.cov.glmm.nb.q.out <- all.taxa.con.cov.glmm.nb(con.var = taxa.results$con.var, id.var = taxa.results$id.var, cov.var = taxa.results$cov.var, taxa = taxa.results$taxa, taxa.results$lib.size, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.con.cov.glmm.nb.q.out
            
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
            
          } else if (input$covariates_taxa.long == "None") {
            taxa.con.glmm.nb.q.out <- all.taxa.con.glmm.nb(con.var = taxa.results$con.var, id.var = taxa.results$id.var, taxa = taxa.results$taxa, taxa.results$lib.size, multi.method = "BH")
            taxa.outputs$DAoutputlong <- taxa.con.glmm.nb.q.out
            
            nrow <- taxa.forest.plot.pages(taxa.outputs$DAoutputlong, species.include = include)
          }
        }
        forestplot.data <- taxa.forest.plot.pages1(taxa.outputs$DAoutputlong, chooseData$taxa.names.out, species.include = include)
        
        if (any(!is.na(unlist(chooseData$taxa.names.out$duplicates)))) {
          duplicate.taxa <- sapply(strsplit(unlist(chooseData$taxa.names.out$duplicates), " :"),  "[", 1)
          taxon.inplot <- unlist(lapply(forestplot.data$all.text.tab, `[`, i =, j = 3))
          duplicate.texts <- sum(duplicate.taxa %in% taxon.inplot)
        } else {
          duplicate.texts <- 0
        }
        
        incProgress(2/10, message = "Displaying Results in progress")
        
        if (duplicate.texts >0) {
          output$taxa_displaylong = renderUI({
            do.call(tabsetPanel, lapply(1:nrow, function(i) {
              tabPanel(title = paste0("Page ", i), align = "center",
                       plotOutput(paste0("forestlong", i), height = 800, width = 750), 
                       plotOutput(paste0("duplicateslong", i), height = 11.5*duplicate.texts+10, width = 750))
            }))
          })
        } else {
          output$taxa_displaylong = renderUI({
            do.call(tabsetPanel, lapply(1:nrow, function(i) {
              tabPanel(title = paste0("Page ", i), align = "center",
                       plotOutput(paste0("forestlong", i), height = 800, width = 750))
            }))
          })
        }
        
        
        lapply(1:nrow, function(j) {
          output[[paste0("forestlong", j)]] <- renderPlot({
            taxa.forest.plot.pages2(page.taxa.q.out = forestplot.data, page = j)
          })
        })
        
        lapply(1:nrow, function(j) {
          output[[paste0("duplicateslong", j)]] <- renderPlot({
            duplicate.list(duplicate.taxa, taxon.inplot, chooseData$taxa.names.out$duplicates)
          })
        })
        
        taxa.sig <- taxa.sig.dend(taxa.outputs$DAoutputlong, chooseData$NAadded$tax.tab, "twopi", include)
        
        flow.text <- taxa.sig$flow.text
        taxon.tab <- taxa.sig$taxon.tab
        ci.tab.all <- taxa.sig$ci.tab.all
        
        if ( include == FALSE){
          ci.tab.all <- ci.tab.all[-(grep("s_", taxon.tab$Taxon, fixed = TRUE)+1)]
          taxon.tab <- taxon.tab[ !grepl("s_", taxon.tab$Taxon, fixed = TRUE), ]
        }
        
        if ( length(ci.tab.all) > 1 ){
          
          for( i in 1:nrow(taxon.tab)){
            if ( ci.tab.all[-1][i] < 0){
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "blue")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
            else{
              taxon.tab[i,1] <- cell_spec(taxon.tab[i,1], "html", color = "red")
              taxon.tab[i,2] <- cell_spec(taxon.tab[i,2], "html", color = "black")
            }
            
          }
        }
        
        N <- dim(taxon.tab)[1]
        itr <- ceiling(N/5)
        tab.one   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.two   <- data.frame( matrix(ncol=2,nrow=0) )
        tab.three <- data.frame( matrix(ncol=2,nrow=0) )
        tab.four  <- data.frame( matrix(ncol=2,nrow=0) )
        tab.five  <- data.frame( matrix(ncol=2,nrow=0) )
        
        colnames(tab.one) <- c("ID", "Taxon")
        colnames(tab.two) <- c("ID", "Taxon")
        colnames(tab.three) <- c("ID", "Taxon")
        colnames(tab.four) <- c("ID", "Taxon")
        colnames(tab.five) <- c("ID", "Taxon")
        
        if ( dim(taxon.tab)[1] > 0 ) {
          
          i = 1
          j = 1
          N.rmd <- N %% 5
          N.fix <- N + 5 - N.rmd
          while ( i <= itr) {
            tab.one  [i,] <- taxon.tab[j,]
            tab.two  [i,] <- taxon.tab[j+1,]
            tab.three[i,] <- taxon.tab[j+2,]
            tab.four [i,] <- taxon.tab[j+3,]
            tab.five [i,] <- taxon.tab[j+4,]
            i <- i + 1  
            j <- j + 5
          }
          row.names(tab.one) <- NULL
          row.names(tab.two) <- NULL
          row.names(tab.three) <- NULL
          row.names(tab.four) <- NULL
          row.names(tab.five) <- NULL
          
          tab.one <- na.omit(tab.one)
          tab.two <- na.omit(tab.two)
          tab.three <- na.omit(tab.three)
          tab.four <- na.omit(tab.four)
          tab.five <- na.omit(tab.five)
        }
        
        output$Ltaxa_display_dend = renderUI({
          
          box(title = strong("Dendrogram"), width = 12, status = "primary", solidHeader = TRUE,
              
              fluidRow(width = 12, align = "center",
                       div(style = "display: inline-block:vertical-align:top;", grVizOutput("Ldendrogram", height = 1000, width = 1000)) ),
              br(),
              fluidRow(width = 12, align = "center",
                       tagList(
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lcon1_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lcon2_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lcon3_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lcon4_taxonlist") ),
                         div(style="display: inline-block;vertical-align:top;", htmlOutput("Lcon5_taxonlist") )
                       )
              )
          )
        })
        
        output$Ldendrogram = renderGrViz({
          flow.text
        })
        output$Lcon1_taxonlist <- renderText({
          sig.tab1 <- kable(tab.one, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab1
        })
        output$Lcon2_taxonlist <- renderText({
          
          sig.tab2 <- kable(tab.two, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab2
        })
        output$Lcon3_taxonlist <- renderText({
          
          sig.tab3 <- kable(tab.three, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab3
        })
        output$Lcon4_taxonlist <- renderText({
          
          sig.tab4 <- kable(tab.four, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab4
        })
        output$Lcon5_taxonlistt <- renderText({
          
          sig.tab5 <- kable(tab.five, 'html', booktabs =TRUE, escape = FALSE) %>%
            kable_styling(latex_options = c('hold_position'))
          sig.tab5
        })
        
        incProgress(1/10, message = "Displaying Results in progress")
        output$downloadTable_taxalong = renderUI({
          tagList(
            box(title = strong("Download Output Table", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                p("You can download the summary statistics and data analysis outputs.",
                  style = "font-size:11pt"),
                h5("Summary Statistics"),
                downloadButton("tdownloadTabl1long", "Download", width = '50%', style = "color:black; background-color: red3"), br(),
                h5("Data Analysis Outputs"),
                downloadButton("tdownloadTabl2long", "Download", width = '50%', style = "color:black; background-color: red3"),
                h5("Dendrogram"),
                downloadButton("gdownloadlong", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        incProgress(1/10, message = "Displaying Results in progress")
        output$tdownloadTabl1long <- downloadHandler(
          filename = function() {
            paste("Taxa.Sum.Table.zip")
          },
          content = function(sum.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.results$taxa.con.sum.out$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=sum.file, files=dataFiles)
          }
        )
        
        output$tdownloadTabl2long <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Output.zip")
          },
          content = function(DA.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(taxa.outputs$DAoutputlong$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=DA.file, files=dataFiles)
          }
        )
        
        output$gdownloadlong <- downloadHandler(
          filename = function() {
            paste("Taxa.Analysis.Graphical.Output", ".html", sep="")
          },
          content = function(file) {
            htmlwidgets::saveWidget(as_widget(taxa.sig.dend(taxa.outputs$DAoutputlong, chooseData$tax.tab, "twopi", include)), file)
          }
        )
        
        ref_string = REFERENCE_CHECK(data_transform = input$dataType_taxa.long, method_name = isolate(input$chooseMethod_taxa.long), FDR = "Yes")
        if (is.null(ref_string)) {
          shinyjs::hide("taxa_references_long")
        } else {
          shinyjs::show("taxa_references_long")
          output$taxa_references_long = renderUI({
            tagList(
              box(title = strong("References", style = "color:black"), width = NULL, status = "primary", solidHeader = TRUE,
                  HTML(paste(ref_string, collapse="<br/>"))
              )
            )
          })
        }
        shinyjs::enable("runbtn_cont_taxa.long")
        shinyjs::enable("dataType_taxa.long")
        shinyjs::enable("primvar_taxa.long")
        shinyjs::enable("chooseMethod_taxa.long")
        shinyjs::enable("covariates_taxa.long")
        shinyjs::enable("rename.con.var_taxa")
        shinyjs::enable("include_species.dend.long")
        shinyjs::enable("clustervar_taxa")
      }
    )
  }, ignoreNULL = TRUE, ignoreInit = TRUE)  
}

shinyApp(ui = ui, server = server)
