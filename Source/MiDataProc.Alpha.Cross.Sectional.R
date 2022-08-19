######################################
# Quality control and transformation #
######################################

library(phangorn)
library(phyloseq)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
library(quantreg)
library(fossil)
library(picante)
library(entropart)

rem.tax.d <- c("", "metagenome", "gut metagenome", "mouse gut metagenome")
rem.tax.str.d <- c("uncultured", "incertae", "Incertae", "unidentified", "unclassified", "unknown")

tax.tab.clean <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- na.code
    tax.tab.c[is.element(taxa, rem.tax), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.str, collapse="|"), taxa), i] <- na.code
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  ind <- which(tax.tab.c[,1] != na.code)
  tax.tab.c <- tax.tab.c[ind,]
  
  tax.tab.h <- tax.tab.c
  
  ind <- which(tax.tab.h[,1] != na.code)
  tax.tab.h[ind ,1] <- paste("k_", tax.tab.h[ind ,1], sep = "")
  ind <- which(tax.tab.h[,2] != na.code)
  ind_omit <- which(tax.tab.h[,2] != na.code & tax.tab.h[,1] == na.code)
  tax.tab.h[ind ,2] <- paste(tax.tab.h[ind,1],paste("p_",tax.tab.h[ind,2], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(2:7)] = na.code
  ind <- which(tax.tab.h[,3] != na.code)
  ind_omit <- which(tax.tab.h[,3] != na.code & tax.tab.h[,2] == na.code)
  tax.tab.h[ind ,3] <- paste(tax.tab.h[ind,2],paste("c_",tax.tab.h[ind,3], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(3:7)] = na.code
  ind <- which(tax.tab.h[,4] != na.code)
  ind_omit <- which(tax.tab.h[,4] != na.code & tax.tab.h[,3] == na.code)
  tax.tab.h[ind ,4] <- paste(tax.tab.h[ind,3],paste("o_",tax.tab.h[ind,4], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(4:7)] = na.code
  ind <- which(tax.tab.h[,5] != na.code)
  ind_omit <- which(tax.tab.h[,5] != na.code & tax.tab.h[,4] == na.code)
  tax.tab.h[ind ,5] <- paste(tax.tab.h[ind,4],paste("f_",tax.tab.h[ind,5], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(5:7)] = na.code
  ind <- which(tax.tab.h[,6] != na.code)
  ind_omit <- which(tax.tab.h[,6] != na.code & tax.tab.h[,5] == na.code)
  tax.tab.h[ind ,6] <- paste(tax.tab.h[ind,5],paste("g_",tax.tab.h[ind,6], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(6:7)] = na.code
  ind <- which(tax.tab.h[,7] != na.code)
  ind_omit <- which(tax.tab.h[,7] != na.code & tax.tab.h[,6] == na.code)
  tax.tab.h[ind ,7] <- paste(tax.tab.h[ind,6],paste("s_",tax.tab.h[ind,7], sep = ""), sep = ";")
  if(length(ind_omit)!=0) tax.tab.h[ind_omit,7] = na.code
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}

tax.tab.cleanNA <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.c <- tax.tab
  # print(head(tax.tab.c))
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- "NA "
    tax.tab.c[is.element(taxa, rem.tax), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.str, collapse="|"), taxa), i] <- na.code
    
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  
  ind <- which(tax.tab.c[,1] != na.code)
  tax.tab.c <- tax.tab.c[ind,]
  tax.tab.h <- tax.tab.c
  
  ind <- which(tax.tab.h[,1] != na.code)
  tax.tab.h[ind ,1] <- paste("k_", tax.tab.h[ind ,1], sep = "")
  
  ranks <- c("p_", "c_", "o_", "f_", "g_", "s_")
  for (i in 1:6) {
    na.num <- 1
    ind <- which(tax.tab.h[,i+1] != na.code)
    ind_omit <- which(tax.tab.h[,i+1] != na.code & tax.tab.h[,i] == na.code)
    tax.tab.h[ind ,i+1] <- paste(tax.tab.h[ind,i],paste(ranks[i],tax.tab.h[ind,i+1], sep = ""), sep = ";")
    # tax.tab.h[ind ,i+1] <- apply(as.matrix(tax.tab.h[ind ,i+1]), 1, function(x) 
    #   if(substring(x, nchar(x)-2) == "_NA") {
    #     x <- paste0(x, na.num); na.num <- na.num +1
    #     print(na.num)
    #   } else { 
    #     x <- x
    #   })
    # for(k in 1:length(tax.tab.h[ind, i+1])){
    #   if(substring(tax.tab.h[ind, i+1][k], nchar(tax.tab.h[ind, i+1][k])-2) == "_NA") {
    #     print(tax.tab.h[ind, i+1][k])
    #     tax.tab.h[ind, i+1] <- paste0(substring(tax.tab.h[ind, i+1], na.num))
    #     na.num <- na.num + 1          
    #   }
    # }
    # if(substring(tax.tab.h[ind ,i+1], nchar(tax.tab.h[ind ,i+1])-2))
    
    if(length(ind_omit)!=0) tax.tab.h[ind_omit,c((i+1):7)] = na.code
  }
  
  # ind <- which(tax.tab.h[,2] != na.code)
  # ind_omit <- which(tax.tab.h[,2] != na.code & tax.tab.h[,1] == na.code)
  # tax.tab.h[ind ,2] <- paste(tax.tab.h[ind,1],paste("p_",tax.tab.h[ind,2], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(2:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,3] != na.code)
  # ind_omit <- which(tax.tab.h[,3] != na.code & tax.tab.h[,2] == na.code)
  # tax.tab.h[ind ,3] <- paste(tax.tab.h[ind,2],paste("c_",tax.tab.h[ind,3], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(3:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,4] != na.code)
  # ind_omit <- which(tax.tab.h[,4] != na.code & tax.tab.h[,3] == na.code)
  # tax.tab.h[ind ,4] <- paste(tax.tab.h[ind,3],paste("o_",tax.tab.h[ind,4], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(4:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,5] != na.code)
  # ind_omit <- which(tax.tab.h[,5] != na.code & tax.tab.h[,4] == na.code)
  # tax.tab.h[ind ,5] <- paste(tax.tab.h[ind,4],paste("f_",tax.tab.h[ind,5], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(5:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,6] != na.code)
  # ind_omit <- which(tax.tab.h[,6] != na.code & tax.tab.h[,5] == na.code)
  # tax.tab.h[ind ,6] <- paste(tax.tab.h[ind,5],paste("g_",tax.tab.h[ind,6], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,c(6:7)] = na.code
  # 
  # ind <- which(tax.tab.h[,7] != na.code)
  # ind_omit <- which(tax.tab.h[,7] != na.code & tax.tab.h[,6] == na.code)
  # tax.tab.h[ind ,7] <- paste(tax.tab.h[ind,6],paste("s_",tax.tab.h[ind,7], sep = ""), sep = ";")
  # if(length(ind_omit)!=0) tax.tab.h[ind_omit,7] = na.code
  
  tax.tab.hh <- tax.tab.h
  
  return(tax.tab.hh)
}

otu.tab.clean <- function(biom, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  lib.size <- colSums(otu.tab)
  ind.low.lib.size <- which(lib.size > lib.size.cut.off)
  lib.size <- lib.size[ind.low.lib.size]
  otu.tab <- otu.tab[,ind.low.lib.size]
  sam.dat <- sam.dat[ind.low.lib.size,]
  
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  ind.low.mean.prop <- which(mean.prop > mean.prop.cut.off)
  taxa <- rownames(otu.tab)[ind.low.mean.prop]
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  biom <- prune_taxa(taxa, biom)
  
  return(biom)  
}

biom.clean <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  tax.tab <- tax_table(biom)
  
  if (kingdom != "all") {
    
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  ind.com.otu <- intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  tax.tab <- tax.tab[ind.com.2,]
  tree <- prune_taxa(ind.com.otu, tree)
  if(!is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  }
  if (tax.tab.c) {
    tax.tab <- tax.tab.clean(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

biom.cleanNA <- function(biom, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, kingdom = "Bacteria", tax.tab.c = TRUE, lib.size.cut.off = 3000, mean.prop.cut.off = 2e-05) {
  
  tax.tab <- tax_table(biom)
  
  if (kingdom != "all") {
    
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  ind.com.otu <- intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  tax.tab <- tax.tab[ind.com.2,]
  tree <- prune_taxa(ind.com.otu, tree)
  if(!is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  }
  if (tax.tab.c) {
    tax.tab <- tax.tab.cleanNA(tax.tab, rem.tax = rem.tax, rem.tax.str = rem.tax.str)
  }
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}

num.tax.rank <- function(tax.tab, rem.tax = rem.tax.d, rem.tax.str = rem.tax.str.d, na.code = "NANANA") {
  tax.tab.cleaned <- tax.tab.clean(tax.tab, rem.tax, rem.tax.str, na.code = na.code)
  num.taxa <- c()
  for (i in 1:6) {
    taxa <- unique(tax.tab.cleaned[,i+1])
    uni.taxa <- sum(taxa == na.code)
    num.taxa[i] <- nrow(taxa) - uni.taxa
  }
  return(num.taxa)
}

lib.size.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  lib.size.sum <- c(mean(lib.size), quantile(lib.size))
  names(lib.size.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(lib.size = lib.size, lib.size.sum = lib.size.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

mean.prop.func <- function(biom) {
  otu.tab <- otu_table(biom)
  lib.size <- colSums(otu.tab)
  prop.otu.tab <- otu.tab
  for (i in 1:length(lib.size)) {
    prop.otu.tab[,i] <- otu.tab[,i]/lib.size[i]
  }
  mean.prop <- rowMeans(prop.otu.tab)
  mean.prop.sum <- c(mean(mean.prop), quantile(mean.prop))
  names(mean.prop.sum) <- c("Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(list(mean.prop = mean.prop, mean.prop.sum = mean.prop.sum, num.sams = ncol(otu.tab), num.otus = nrow(otu.tab)))
}

rarefy.func <- function(biom, cut.off, multi.rarefy = FALSE) {
  
  if (!multi.rarefy | multi.rarefy == 1) {
    biom <- rarefy_even_depth(biom, cut.off, rngseed = 487)
  } else {
    otu.tab <- otu_table(biom)
    tax.tab <- tax_table(biom)
    tree <- phy_tree(biom)
    sam.dat <- sample_data(biom)
    
    otu.tab.list <- list()
    for (i in 1:multi.rarefy) {
      otu.tab.list[[i]] <- otu_table(rarefy_even_depth(biom, cut.off, rngseed = i), taxa_are_rows = TRUE)
    }
    
    sum.otu.tab <- otu.tab.list[[1]]
    for (i in 2:multi.rarefy) {
      sum.otu.tab <- sum.otu.tab + otu.tab.list[[i]]
    }
    otu.tab <- otu_table(round(sum.otu.tab/multi.rarefy), taxa_are_rows = TRUE)
    biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat) 
  }
  
  return(biom)
}

###################
# Alpha diversity #
###################
alpha.pe.pqe.func <- function(x, tree, norm = TRUE) {
  ind <- which(x != 0)
  s.tree <- prune_taxa(names(x[ind]), tree)
  pe <- AllenH(x[ind], 1, s.tree, Normalize = norm, CheckArguments = FALSE)
  pqe <- AllenH(x[ind], 2, s.tree, Normalize = norm, CheckArguments = FALSE)
  return(c(pe, pqe))
}

alpha.v1.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE")], alpha.ice, alpha.pd))
  colnames(alpha.div) <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
  
  return(alpha.div)
}

alpha.v2.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  alpha.pe.pqe <- t(apply(t.prop.otu.tab, 1, function(x) alpha.pe.pqe.func(x, tree)))
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "InvSimpson", "Chao1", "ACE")], alpha.ice, alpha.pd, alpha.pe.pqe))
  colnames(alpha.div) <- c("Observed", "Shannon", "InvSimpson", "Chao1", "ACE", "ICE", "PD", "PE", "PQE")
  
  return(alpha.div)
}

############################
# Alpha diversity analysis #
############################

##########################
# Binary - No covariates #
##########################

alpha.bin.t.test <- function(bin.var, alpha.div) {
  n.alpha <- ncol(alpha.div)
  out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    fit <- t.test(alpha.div[,i] ~ bin.var)
    out[i,] <- c(fit$statistic, fit$stderr, fit$parameter, fit$conf.int, fit$p.value)   
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha.div)
  colnames(out) <- c("t", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(out)
}

alpha.bin.wilcox.test <- function(bin.var, alpha.div) {
  n.alpha <- ncol(alpha.div)
  out <- matrix(NA, n.alpha, 2)
  for (i in 1:n.alpha) {
    fit <- wilcox.test(alpha.div[,i] ~ bin.var)
    out[i,] <- c(fit$statistic, fit$p.value)   
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha.div)
  colnames(out) <- c("W", "P.value")
  return(out)
}

alpha.bin.hist <- function(bin.var, alpha.div, test.out, mult.test.cor = TRUE) {
  
  n.alpha <- ncol(alpha.div)
  alpha.div.ind <- colnames(alpha.div)
  pvs <- test.out$P.value
  ind.p.sig <- which(pvs < 0.05)
  qvs <- test.out$Q.value
  ind.q.sig <- which(qvs < 0.05)
  
  par(mfrow = c(3, 3))
  for (i in 1:n.alpha) {
    if (!mult.test.cor) {
      if (is.element(i, ind.p.sig)) {
        xlab.v = paste("*p:", p.value.0.1(pvs[i]), sep="")
      } else {
        xlab.v = paste("p:", p.value.0.1(pvs[i]), sep="")
      }
    }
    if (mult.test.cor) {
      if (is.element(i, ind.q.sig)) {
        xlab.v = paste("*q:", p.value.0.1(qvs[i]), sep="")
      } else {
        xlab.v = paste("q:", p.value.0.1(qvs[i]), sep="")
      }
    }
    print(levels(bin.var))
    boxplot(alpha.div[,i] ~ bin.var, xlab=xlab.v, ylab=alpha.div.ind[i], names = substr(levels(bin.var), 1, 8), notch = TRUE, col=c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)))
  }
}

alpha.bin.var.func <- function(sam.dat) {
  var.names <- colnames(sam.dat)
  return(var.names)
}

###################
# Which variable? #
###################

# Which primary variable do you want to select?

is.mon.sin.rev.bin.con <- function(sam.dat) {
  
  n.var <- ncol(sam.dat)
  n.sam <- nrow(sam.dat)
  is.mon <- logical()
  is.rev <- logical()
  is.bin <- logical()
  is.con <- logical()
  is.sin <- logical()
  
  for (i in 1:n.var) {
    sam.var <- as.matrix(sam.dat[,i])
    if (length(table(sam.var)) == 1) {
      is.mon[i] <- TRUE
    }
    if (length(table(sam.var)) != 1) {
      is.mon[i] <- FALSE
    }
    if (length(table(sam.var)) == 2 & any(table(sam.var)==1)) {
      is.sin[i] <- TRUE
    }
    if (length(table(sam.var)) != 2 | !any(table(sam.var)==1)) {
      is.sin[i] <- FALSE
    }
    if (length(table(sam.var)) == n.sam & sum(is.na(as.numeric(sam.var))) == n.sam) {
      is.rev[i] <- TRUE
    }
    if (length(table(sam.var)) != n.sam | sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.rev[i] <- FALSE
    }
    if (length(table(sam.var)) == 2) {
      is.bin[i] <- TRUE
    }
    if (length(table(sam.var)) != 2) {
      is.bin[i] <- FALSE
    }
    if (length(table(sam.var)) != 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.con[i] <- TRUE
    }
    if (length(table(sam.var)) == 2 & sum(is.na(as.numeric(sam.var))) != n.sam) {
      is.con[i] <- FALSE
    }
    if (sum(is.na(as.numeric(sam.var))) == n.sam) {
      is.con[i] <- FALSE
    }
    
  }
  return(list(is.mon = is.mon, is.sin = is.sin, is.rev = is.rev, is.bin = is.bin, is.con = is.con))
}

pri.func <- function(sam.dat, mon.sin.rev.bin.con) {
  colnames(sam.dat)[(mon.sin.rev.bin.con$is.bin | mon.sin.rev.bin.con$is.con) & !mon.sin.rev.bin.con$is.mon & !mon.sin.rev.bin.con$is.sin]
}

is.bin.con.pri <- function(sam.dat, mon.sin.rev.bin.con, sel.pri.var) {
  ind <- which(colnames(sam.dat) == sel.pri.var)
  if(length(ind) != 0){
    if (mon.sin.rev.bin.con$is.bin[ind]) {
      out <- "Binary"
    } else {
      out <- "Continuous"
    }
  }else {
    out = "Neither"
  }
  return(out)
}

# What covariate(s) do you want to select?

cov.func <- function(sam.dat, mon.sin.rev.bin.con, sel.pri.var) {
  ind.pri <- colnames(sam.dat) == sel.pri.var
  ind.mon.sin.rev <- mon.sin.rev.bin.con$is.mon | mon.sin.rev.bin.con$is.rev | mon.sin.rev.bin.con$is.sin
  return(colnames(sam.dat)[!(ind.pri | ind.mon.sin.rev)])
}

alpha.bin.cat.func <- function(sam.dat, sel.bin.var) {
  # bin.var <- unlist(sam.dat[,sel.bin.var])
  # bin.var.no.na <- bin.var[!is.na(bin.var)]
  # bin.cat <- unique(bin.var.no.na)
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  
  return(bin.cat)
}

is.binary <- function(sam.dat, sel.bin.var) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  bin.var.no.na <- bin.var[!is.na(bin.var)]
  bin.cat <- unique(bin.var.no.na)
  if (length(bin.cat) != 2) {
    return(FALSE)
  }
  return(TRUE)
}

alpha.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

alpha.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  return(sam.dat)
}

alpha.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, alpha.div) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  alpha.div <- rbind(alpha.div[ind.ref,], alpha.div[ind.com,])
  return(list(bin.var = bin.var, alpha.div = alpha.div))
}

alpha.ind.sum.func <- function(x) {
  sum.out <- c(length(x), mean(x, na.rm = TRUE), quantile(x, na.rm = TRUE))
  names(sum.out) <-  c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(sum.out)
}

alpha.bin.sum.func <- function(bin.var, alpha.div) {
  n.alpha <- ncol(alpha.div)
  ref.sum <- matrix(NA, n.alpha, 7)
  com.sum <- matrix(NA, n.alpha, 7)
  print(1)
  for (i in 1:n.alpha) {
    ind.alpha <- alpha.div[,i]
    sum.out <- tapply(ind.alpha, bin.var, alpha.ind.sum.func)
    ref.sum[i,] <- sum.out[[1]]
    com.sum[i,] <- sum.out[[2]]
    print(2)
  }
  rownames(ref.sum) <- colnames(alpha.div)
  colnames(ref.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  rownames(com.sum) <- colnames(alpha.div)
  colnames(com.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  out <- list(as.data.frame(ref.sum), as.data.frame(com.sum))
  names(out) <- levels(bin.var)
  return(out)
}

#######################
# Binary - Covariates #
#######################

alpha.lm.bin.cov.func <- function(bin.var, cov.var, alpha.div, scale = TRUE) {
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, cov.var, alpha.div))
  
  lm.out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      lm.f <- formula(paste("scale(", alpha.ind[i], ")", "~", colnames(bin.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    if (!scale) {
      lm.f <- formula(paste(alpha.ind[i], "~", colnames(bin.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    fit.lm <- lm(lm.f, data = d)
    out.lm <- c(summary(fit.lm)$coefficients[2,c(1,2)], summary(fit.lm)$df[2], confint(fit.lm)[2,], summary(fit.lm)$coefficients[2,4])
    lm.out[i,] <- out.lm
  }
  lm.out <- as.data.frame(lm.out)
  rownames(lm.out) <- colnames(alpha.div)
  colnames(lm.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(lm.out)
}

alpha.bin.cov.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sel.cov.var, sam.dat, alpha.div) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  alpha.div <- rbind(alpha.div[ind.ref,], alpha.div[ind.com,])
  
  return(list(bin.var = bin.var, cov.var = cov.var, alpha.div = alpha.div))
}

alpha.forest.plot <- function(out, mult.test.cor = TRUE) {
  
  if (mult.test.cor) {
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "Est", "SE", "DF", "P-value", "Q-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), out[, 3], p.value.0.1(out[,6]), p.value.0.1(out[,7]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=0, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }else{
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "Est", "SE", "DF", "P-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), out[, 3], p.value.0.1(out[,6]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=0, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }
}

get.or.se <- function(model) {
  broom::tidy(model) %>%
    mutate(or = exp(estimate),
           var.diag = diag(vcov(model)),
           or.se = sqrt(or^2 * var.diag)) %>%
    dplyr::select(or.se) %>% unlist %>% unname
}

alpha.logit.reg.coef.bin.cov.func <- function(bin.var, cov.var, alpha.div, scale = TRUE) {
  
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, cov.var, alpha.div))
  
  logit.out <- matrix(NA, n.alpha, 6)
  
  for (i in 1:n.alpha) {
    if (scale) {
      logit.f <- formula(paste(colnames(bin.var), "~", "scale(", alpha.ind[i], ")", "+", paste(colnames(cov.var), collapse = "+")))
    }
    if (!scale) {
      logit.f <- formula(paste(colnames(bin.var), "~", alpha.ind[i], "+", paste(colnames(cov.var), collapse = "+")))
    }
    fit.logit <- glm(logit.f, data = d, family = "binomial")
    
    est <- summary(fit.logit)$coefficients[2,1]
    std.err <- summary(fit.logit)$coefficients[2,2]
    df <- summary(fit.logit)$df[2]
    ci <- c(est - qt(0.975, df)*std.err, est + qt(0.975, df)*std.err)
    out.logit <- c(est, std.err, df, ci, summary(fit.logit)$coefficients[2,4])
    logit.out[i,] <- out.logit
  }
  logit.out <- as.data.frame(logit.out)
  rownames(logit.out) <- colnames(alpha.div)
  colnames(logit.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  
  return(logit.out)
}

alpha.logit.bin.cov.func <- function(bin.var, cov.var, alpha.div, scale = TRUE) {
  
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, cov.var, alpha.div))
  
  logit.out <- matrix(NA, n.alpha, 6)
  
  for (i in 1:n.alpha) {
    
    if (scale) {
      logit.f <- formula(paste(colnames(bin.var), "~", "scale(", alpha.ind[i], ")", "+", paste(colnames(cov.var), collapse = "+")))
    }
    if (!scale) {
      logit.f <- formula(paste(colnames(bin.var), "~", alpha.ind[i], "+", paste(colnames(cov.var), collapse = "+")))
    }
    fit.logit <- glm(logit.f, data = d, family = "binomial")
    est <- summary(fit.logit)$coefficients[2,1]
    std.err <- summary(fit.logit)$coefficients[2,2]
    or.se <- sqrt(exp(est)^2*diag(vcov(fit.logit)))[2]
    
    df <- summary(fit.logit)$df[2]
    ci <- c(exp(est - qt(0.975, df)*std.err), exp(est + qt(0.975, df)*std.err))
    out.logit <- c(exp(est), or.se, df, ci, summary(fit.logit)$coefficients[2,4])
    
    logit.out[i,] <- out.logit
  }
  logit.out <- as.data.frame(logit.out)
  rownames(logit.out) <- colnames(alpha.div)
  colnames(logit.out) <- c("OR", "Std Err", "DF", "Lower", "Upper", "P.value")
  
  return(logit.out)
}

alpha.logit.forest.plot <- function(out, mult.test.cor = TRUE) {
  
  if (mult.test.cor) {
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "OR", "SE", "DF", "P-value", "Q-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), out[, 3], p.value.0.1(out[,6]), p.value.0.1(out[,7]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=1, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for OR",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }else{
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "OR", "SE", "DF", "P-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), out[, 3], p.value.0.1(out[,6]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=1, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for OR",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }
}

alpha.qr.bin.cov.func <- function(bin.var, cov.var, alpha.div, tau = 0.5, qr.method = "br", se.method = "nid", n.res = 1000, scale = TRUE) {
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, cov.var, alpha.div))
  
  qr.out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      qr.f <- formula(paste("scale(", alpha.ind[i], ")", "~", colnames(bin.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    if (!scale) {
      qr.f <- formula(paste(alpha.ind[i], "~", colnames(bin.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    fit.qr <- rq(qr.f, data = d, tau = tau, method = qr.method)
    if (se.method == "boot") {
      set.seed(i)
    }
    out.coef <- summary(fit.qr, se = se.method, R = n.res)$coefficients
    out.df <- summary(fit.qr, se = se.method, R = n.res)$rdf
    out.qr <- c(out.coef[2,c(1,2)], out.df, out.coef[2,1]-qt(0.975,out.df)*out.coef[2,2], out.coef[2,1]+qt(0.975,out.df)*out.coef[2,2], out.coef[2,4])
    qr.out[i,] <- out.qr
  }
  qr.out <- as.data.frame(qr.out)
  rownames(qr.out) <- colnames(alpha.div)
  colnames(qr.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(qr.out)
}

#############################
# Continuous - No covariate #
#############################

alpha.con.recode.func <- function(sam.dat, sel.con.var, rename.con.var, alpha.div) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  alpha.div = alpha.div[ind.nona,]
  return(list(con.var = con.var, alpha.div = alpha.div))
}

alpha.lm.con.func <- function(con.var, alpha.div, scale = TRUE) {
  n.alpha <- ncol(alpha.div)
  lm.out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      fit.lm <- lm(scale(alpha.div[,i]) ~ con.var[,1])
    }
    if (!scale) {
      fit.lm <- lm(alpha.div[,i] ~ con.var[,1])
    }
    out.lm <- c(summary(fit.lm)$coefficients[2,c(1,2)], summary(fit.lm)$df[2], confint(fit.lm)[2,], summary(fit.lm)$coefficients[2,4])
    lm.out[i,] <- out.lm
  }
  lm.out <- as.data.frame(lm.out)
  rownames(lm.out) <- colnames(alpha.div)
  colnames(lm.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(lm.out)
}

alpha.con.plot <- function(alpha.con.out, test.out, mult.test.cor = TRUE) {
  
  pvs <- test.out$P.value
  ind.p.sig <- which(pvs < 0.05)
  qvs <- test.out$Q.value
  ind.q.sig <- which(qvs < 0.05)
  n.alpha <- length(pvs)
  
  par(mfrow = c(3, 3))
  for (i in 1:n.alpha) {
    if (!mult.test.cor) {
      if (is.element(i, ind.p.sig)) {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\n*p:", p.value.0.1(pvs[i]), sep="")
      } else {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\np:", p.value.0.1(pvs[i]), sep="")
      }
    }
    if (mult.test.cor) {
      if (is.element(i, ind.q.sig)) {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\n*q:", p.value.0.1(qvs[i]), sep="")
      } else {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\nq:", p.value.0.1(qvs[i]), sep="")
      }
    }
    x <- alpha.con.out$con.var[,1]
    y <- alpha.con.out$alpha.div[,i]
    plot(x, y, xlab = xlab.v, ylab = colnames(alpha.con.out$alpha.div)[i], pch = 20, col = "gray")
    fit.lm <- lm(y ~ x)
    new.x = seq(min(x), max(x), by = (max(x)-min(x))/1000)
    conf.int <- predict(fit.lm, newdata = data.frame(x = new.x), interval = "confidence", level = 0.95)
    abline(fit.lm, col = "lightblue2", lwd = 2)
    lines(new.x, conf.int[,2], col = "blue", lwd = 2, lty="dotted")
    lines(new.x, conf.int[,3], col = "blue", lwd = 2, lty="dotted")
  }
}

alpha.qr.con.func <- function(con.var, alpha.div, tau = 0.5, qr.method = "br", se.method = "nid", n.res = 1000, scale = TRUE) {
  n.alpha <- ncol(alpha.div)
  qr.out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      fit.qr <- rq(scale(alpha.div[,i]) ~ con.var[,1], tau = tau, method = qr.method)
    }
    if (!scale) {
      fit.qr <- rq(alpha.div[,i] ~ con.var[,1], tau = tau, method = qr.method)
    }
    if (se.method == "boot") {
      set.seed(i)
    }
    out.coef <- summary(fit.qr, se = se.method, R = n.res)$coefficients
    out.df <- summary(fit.qr, se = se.method, R = n.res)$rdf
    out.qr <- c(out.coef[2,c(1,2)], out.df, out.coef[2,1]-qt(0.975,out.df)*out.coef[2,2], out.coef[2,1]+qt(0.975,out.df)*out.coef[2,2], out.coef[2,4])
    qr.out[i,] <- out.qr
  }
  qr.out <- as.data.frame(qr.out)
  rownames(qr.out) <- colnames(alpha.div)
  colnames(qr.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(qr.out)
}

alpha.qr.forest.plot <- function(alpha.con.out, test.out, mult.test.cor = TRUE, tau = 0.5, qr.method = "br", se.method = "nid", n.res = 1000) {
  pvs <- test.out$P.value
  ind.p.sig <- which(pvs < 0.05)
  qvs <- test.out$Q.value
  ind.q.sig <- which(qvs < 0.05)
  n.alpha <- length(pvs)
  par(mfrow = c(3, 3))
  for (i in 1:n.alpha) {
    if (!mult.test.cor) {
      if (is.element(i, ind.p.sig)) {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\n*p:", p.value.0.1(pvs[i]), sep="")
      } else {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\np:", p.value.0.1(pvs[i]), sep="")
      }
    }
    if (mult.test.cor) {
      if (is.element(i, ind.q.sig)) {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\n*q:", p.value.0.1(qvs[i]), sep="")
      } else {
        xlab.v = paste(colnames(alpha.con.out$con.var), "\nq:", p.value.0.1(qvs[i]), sep="")
      }
    }
    x <- alpha.con.out$con.var[,1]
    y <- alpha.con.out$alpha.div[,i]
    plot(x, y, xlab = xlab.v, ylab = colnames(alpha.con.out$alpha.div)[i], pch = 20, col = "gray")
    
    fit.qr <- rq(y ~ x, tau = tau, method = qr.method)
    new.x = seq(min(x), max(x), by = (max(x)-min(x))/length(x))
    conf.int <- predict(fit.qr, newdata = data.frame(x = new.x), interval = "confidence", level = 0.95)
    abline(fit.qr, col = "lightblue2", lwd = 2)
    lines(new.x, conf.int[,2], col = "blue", lwd = 2, lty="dotted")
    lines(new.x, conf.int[,3], col = "blue", lwd = 2, lty="dotted")
  }
}

###########################
# Continuous - Covariates #
###########################

alpha.con.cov.recode.func <- function(sam.dat, sel.con.var, sel.cov.var, rename.con.var, alpha.div) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona, ind.sel]))
  cov.var <- as.data.frame(sam.dat[ind.nona, sel.cov.var])
  colnames(cov.var) <- sel.cov.var
  alpha.div = alpha.div[ind.nona,]
  return(list(con.var = con.var, cov.var = cov.var, alpha.div = alpha.div))
}

alpha.lm.con.cov.func <- function(con.var, cov.var, alpha.div, scale = TRUE) {
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  d <- cbind(con.var, cov.var, alpha.div)
  
  lm.out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      lm.f <- formula(paste("scale(", colnames(alpha.div)[i], ")", "~", colnames(con.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    if (!scale) {
      lm.f <- formula(paste(colnames(alpha.div)[i], "~", colnames(con.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    fit.lm <- lm(lm.f, data = d)
    out.lm <- c(summary(fit.lm)$coefficients[2,c(1,2)], summary(fit.lm)$df[2], confint(fit.lm)[2,], summary(fit.lm)$coefficients[2,4])
    lm.out[i,] <- out.lm
  }
  lm.out <- as.data.frame(lm.out)
  rownames(lm.out) <- colnames(alpha.div)
  colnames(lm.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(lm.out)
}

alpha.con.cov.qr.func <- function(con.var, cov.var, alpha.div, tau = 0.5, qr.method = "br", se.method = "nid", n.res = 1000, scale = TRUE) {
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(con.var, cov.var, alpha.div))
  
  qr.out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      qr.f <- formula(paste("scale(", alpha.ind[i], ")", "~", colnames(con.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    if (!scale) {
      qr.f <- formula(paste(alpha.ind[i], "~", colnames(con.var), "+", paste(colnames(cov.var), collapse = "+")))
    }
    fit.qr <- rq(qr.f, data = d, tau = tau, method = qr.method)
    if (se.method == "boot") {
      set.seed(i)
    }
    out.coef <- summary(fit.qr, se = se.method, R = n.res)$coefficients
    out.df <- summary(fit.qr, se = se.method, R = n.res)$rdf
    out.qr <- c(out.coef[2,c(1,2)], out.df, out.coef[2,1]-qt(0.975,out.df)*out.coef[2,2], out.coef[2,1]+qt(0.975,out.df)*out.coef[2,2], out.coef[2,4])
    qr.out[i,] <- out.qr
  }
  qr.out <- as.data.frame(qr.out)
  rownames(qr.out) <- colnames(alpha.div)
  colnames(qr.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(qr.out)
}

###################
# Other functions #
###################

q.func <- function(out, method = c("BH", "BY")) {
  Q.value <- p.adjust(out$P.value, method = method)
  return(cbind(out, Q.value))
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000")
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000")
  x[ind.1] <- ">.999"
  return(x)
}
