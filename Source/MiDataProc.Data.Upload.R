library(biomformat)
library(phangorn)
library(bios2mds)

preprocess.tax.tab = function(tax.tab){
  trans.tax.tab <- matrix(NA, nrow(tax.tab), 7)
  tax.list <- strsplit(as.character(tax.tab$Taxon), ";")
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_0__", "", tax.list[[i]][grepl("D_0__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,1] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_1__", "", tax.list[[i]][grepl("D_1__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,2] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_2__", "", tax.list[[i]][grepl("D_2__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,3] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_3__", "", tax.list[[i]][grepl("D_3__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,4] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_4__", "", tax.list[[i]][grepl("D_4__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,5] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_5__", "", tax.list[[i]][grepl("D_5__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,6] <- taxon
    }
  }
  
  for (i in 1:length(tax.list)) {
    taxon <- gsub("D_6__", "", tax.list[[i]][grepl("D_6__", tax.list[[i]])])
    if (length(taxon) != 0) {
      trans.tax.tab[i,7] <- taxon
    }
  }
  rownames(trans.tax.tab) <- tax.tab$Feature.ID
  colnames(trans.tax.tab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  return(trans.tax.tab)
}

biom.check.samples <- function(otu.table, sample.data) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  sam.dat <- sample_data(sample.data)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  
  if (length(intersect(colnames(otu.tab), rownames(sam.dat))) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
#h4("Error: There is no common samples among OTU table and Sample Data")
#h4("Error: There is no common OTUs among OTU table, taxonomic table and tree tip labels")
biom.check.otu <- function(otu.table, tax.table, tre.data) {
  otu.tab <- otu_table(otu.table, taxa_are_rows = TRUE)
  tax.tab <- tax_table(tax.table)
  tree <- phy_tree(tre.data)
  
  if (length(intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))) == 0) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

#methods
REFERENCE_CHECK <- function(data_transform = "", method_name = "", FDR = ""){
  reference_lists <- NULL
  if(data_transform == "CLR (Default)"){
    reference_lists = c(reference_lists, "Aitchison J. The statistical analysis of compositional data. J R Statist Soc B. 1982;44(2):139-77.", "Martin-Fernandez JA, Hron K, Templ M, Filzmoser P, Palarea-Albaladejo J. Bayesian-multiplicative treatment of count zeros in compositional data sets. Statistical Modelling 2015:15(2):134-158.")
  }else if(data_transform == "Proportion"){
    reference_lists = c(reference_lists, "Martin-Fernandez JA, Hron K, Templ M, Filzmoser P, Palarea-Albaladejo J. Bayesian-multiplicative treatment of count zeros in compositional data sets. Statistical Modelling 2015:15(2):134-158.")
  }
  
  if(method_name == "Wilcoxon rank-sum test"){
    reference_lists = c(reference_lists, "Mann HB, Whitney DR. On a test of whether one of two random variables is stochastically larger than the other. Ann Math Statist. 1947;18(1):50-60.")
  }else if(method_name == "LMM"){
    reference_lists = c(reference_lists, "Laird NM, Ware JH. Random-effects models for longitudinal data. Biometrics. 1982:38:963-974.")
  }else if(method_name == "MiRKAT"){
    reference_lists = c(reference_lists, 
                        "Zhao N, Chen J, Carroll IM, Ringel-Kulka T, Epstein MP, Zhou H, Zhou JJ, Ringel Y, Li H, Wu MC. Testing in microbiome-profiling studies with MiRKAT, the microbiome regression-based kernel association test. Am J Hum Genet. 2015:96:5:797-807", 
                        "Wilson N, Zhao N, Zhan X, Koh H, Fu W, Chen J, Li H, Wu MC, Plantinga AM. MiRKAT: kernel machine regression-based global association tests for the microbiome. Bioinformatics. 2020:37(11):1-3.")
  }else if(method_name == "GLMM-MiRKAT"){
    reference_lists = c(reference_lists, 
                        "Koh H, Li Y, Zhan X, Chen J, Zhao N. A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. Front Genet. 2019;458(10).", 
                        "Wilson N, Zhao N, Zhan X, Koh H, Fu W, Chen J, Li H, Wu MC, Plantinga AM. MiRKAT: kernel machine regression-based global association tests for the microbiome. Bioinformatics. 2020;37(11):1-3.")
  }else if(method_name == "LMM"){
    reference_lists = c(reference_lists, "Laird NM, Ware JH. Random-effects models for longitudinal data. Biometrics. 1982:38:963-974.")
  }else if(method_name == "GLMM (Binomial)" | method_name == "GLMM (Negative Binomial)" | method_name == "GLMM (Beta)"){
    reference_lists = c(reference_lists, "Breslow NE, Clayton DG. Approximate inference in generalized linear mixed models. J Am Stat Assoc. 1993:88(421):9-25.")
  }else if(method_name == "GEE (Binomial)"){
    reference_lists = c(reference_lists, "Liang KY, Zeger SL. Longitudinal data analysis using generalized linear models. Biometrika. 1986:73(1):13-22")
  }
  
  if(FDR == "Yes"){
    reference_lists = c(reference_lists, "Benjamini Y, Hochberg Y. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J R Statist 
Soc B. 1995;57(1):289-300.")
  }
  
  if(length(reference_lists) == 0){
    return(NULL)
  }
  else{
    #referece_lists = referece_lists[order(referece_lists)]
    for (i in seq(1:length(reference_lists))){
      reference_lists[i] = paste(i, ". ", reference_lists[i], sep = "")
    }
    return(reference_lists)
  }
}