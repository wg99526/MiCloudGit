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
library(gridGraphics)
library(gridExtra)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(erer)
library(DiagrammeR)
library(stringr)
library(devtools)
library(betareg)

tax.trans <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "NANANA") {
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  tax.arc.sin.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- numeric()
    try(tax.imp.prop <- zCompositions::cmultRepl(tax.count), silent = TRUE)
    if(length(tax.imp.prop) == 0) {
      tax.imp.prop <- matrix(NA, n, length(tax))
      for (i in 1:length(lib.size)) {
        tax.imp.prop[i,] <- (tax.count[i,]+0.1)/lib.size[i]
      }
      rownames(tax.imp.prop) <- colnames(otu.tab)
      colnames(tax.imp.prop) <- tax
    }
    
    tax.clr <- compositions::clr(tax.imp.prop)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    tax.arc.sin <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.arc.sin[i,] <- asin(sqrt(tax.prop[i,]))
    }
    rownames(tax.arc.sin) <- colnames(otu.tab)
    colnames(tax.arc.sin) <- tax
    
    ind <- which(tax == na.code)
    if (length(ind) != 0) {
      tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
      tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
      tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin[,-ind])
    }else{
      tax.count.out[[j]] <- as.data.frame(tax.count)
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count)
      tax.prop.out[[j]] <- as.data.frame(tax.prop)
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop)
      tax.clr.out[[j]] <- as.data.frame(tax.clr)
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin)
    }
  }
  
  names(tax.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.rare.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.imp.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.sub.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.arc.sin.out) <- c("phylum", "class", "order", "family", "genus", "species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.sub.clr.out, arcsin = tax.arc.sin.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.clr.out, arcsin = tax.arc.sin.out))
  }
  
}

tax.trans.na <- function(otu.tab, tax.tab, rare.otu.tab, rare.tax.tab, sub.com = TRUE, na.code = "") {
  
  # tax.tab <- matrix(nrow = nrow(tax.tab.ori), ncol = ncol(tax.tab.ori))
  
  n <- ncol(otu.tab)
  lib.size <- colSums(otu.tab)
  rare.n <- ncol(rare.otu.tab)
  
  tax.count.out <- list()
  tax.rare.count.out <- list()
  tax.prop.out <- list()
  tax.imp.prop.out <- list()
  tax.clr.out <- list()
  tax.sub.clr.out <- list()
  tax.arc.sin.out <- list()
  
  for (j in 1:6) {
    tax <- as.vector(unique(tax.tab[,j+1]))
    
    na.num <- 1
    for(i in 1:length(tax)) {
      na.taxa <- substring(tax[i], nchar(tax[i])-2)
      if(na.taxa == "_NA") {
        new.tax <- paste0(tax[i], na.num)
        ind.na.tab <- which(tax.tab[,j+1] == tax[i])
        tax.tab[ind.na.tab, j+1] <<- new.tax
        rare.tax.tab[ind.na.tab, j+1] <<- new.tax
        tax[i] <- new.tax
        na.num <- na.num + 1
      }
      
    }
    
    
    tax.count <- matrix(NA, n, length(tax))
    for (i in 1:length(tax)) {
      ind.tax <- which(tax.tab[,j+1] == tax[i])
      tax.count[,i] <- colSums(otu.tab[ind.tax,])
    }
    rownames(tax.count) <- colnames(otu.tab)
    colnames(tax.count) <- tax
    
    tax.prop <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.prop[i,] <- tax.count[i,]/lib.size[i]
    }
    rownames(tax.prop) <- colnames(otu.tab)
    colnames(tax.prop) <- tax
    
    tax.imp.prop <- numeric()
    try(tax.imp.prop <- zCompositions::cmultRepl(tax.count), silent = TRUE)
    if(length(tax.imp.prop) == 0) {
      tax.imp.prop <- matrix(NA, n, length(tax))
      for (i in 1:length(lib.size)) {
        tax.imp.prop[i,] <- (tax.count[i,]+0.1)/lib.size[i]
      }
      rownames(tax.imp.prop) <- colnames(otu.tab)
      colnames(tax.imp.prop) <- tax
    }
    
    tax.clr <- compositions::clr(tax.imp.prop)
    
    rare.tax <- as.vector(unique(rare.tax.tab[,j+1]))
    tax.rare.count <- matrix(NA, rare.n, length(rare.tax))
    for (i in 1:length(rare.tax)) {
      ind.tax <- which(rare.tax.tab[,j+1] == rare.tax[i])
      tax.rare.count[,i] <- colSums(rare.otu.tab[ind.tax,])
    }
    rownames(tax.rare.count) <- colnames(rare.otu.tab)
    colnames(tax.rare.count) <- rare.tax
    
    tax.arc.sin <- matrix(NA, n, length(tax))
    for (i in 1:length(lib.size)) {
      tax.arc.sin[i,] <- asin(sqrt(tax.prop[i,]))
    }
    rownames(tax.arc.sin) <- colnames(otu.tab)
    colnames(tax.arc.sin) <- tax
    
    ind <- which(tax == na.code)
    if (length(ind) != 0) {
      tax.count.out[[j]] <- as.data.frame(tax.count[,-ind])
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count[,-ind])
      tax.prop.out[[j]] <- as.data.frame(tax.prop[,-ind])
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop[,-ind])
      tax.clr.out[[j]] <- as.data.frame(tax.clr[,-ind])
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin[,-ind])
    }else{
      tax.count.out[[j]] <- as.data.frame(tax.count)
      tax.rare.count.out[[j]] <- as.data.frame(tax.rare.count)
      tax.prop.out[[j]] <- as.data.frame(tax.prop)
      tax.imp.prop.out[[j]] <- tax.sub.imp.prop <- as.data.frame(tax.imp.prop)
      tax.clr.out[[j]] <- as.data.frame(tax.clr)
      tax.sub.clr.out[[j]] <- as.data.frame(compositions::clr(tax.sub.imp.prop))
      tax.arc.sin.out[[j]] <- as.data.frame(tax.arc.sin)
    }
  }
  
  names(tax.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.rare.count.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.imp.prop.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.sub.clr.out) <- c("phylum", "class", "order", "family", "genus", "species")
  names(tax.arc.sin.out) <- c("phylum", "class", "order", "family", "genus", "species")
  
  if (sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.sub.clr.out, arcsin = tax.arc.sin.out))
  }
  if (!sub.com) {
    return(list(count = tax.count.out, rare.count = tax.rare.count.out, imp.prop = tax.imp.prop.out, prop = tax.prop.out, 
                clr = tax.clr.out, arcsin = tax.arc.sin.out))
  }
  
  
}

add_NA <- function(taxa.out, tax.tab) {
  taxa.out.ori <- taxa.out
  tax.tab <- as.data.frame(tax.tab)
  for(type in 1:length(taxa.out)) {
    na.num <- 1
    for (rank in 1:6) {
      for (i in 1:length(colnames(taxa.out[[type]][[rank]]))) {
        if(substring(colnames(taxa.out[[type]][[rank]])[i], nchar(colnames(taxa.out[[type]][[rank]])[i])-2) == "NA ") {
          colnames(taxa.out[[type]][[rank]])[i] <- paste0(colnames(taxa.out[[type]][[rank]])[i], na.num)
          na.num <- na.num + 1
        }
      }
    }
    for (rank in 1:5) {
      for (i in 1:length(colnames(taxa.out[[type]][[rank]]))){
        colnames(taxa.out[[type]][[rank+1]]) <- str_replace(colnames(taxa.out[[type]][[rank+1]]), colnames(taxa.out.ori[[type]][[rank]])[i], colnames(taxa.out[[type]][[rank]])[i])
      }
    }
    
    for(rank in 1:6) {
      for(i in 1:length(unique(colnames(taxa.out[[type]][[rank]])))) {
        ind <- tax.tab[[rank+1]] == unique(colnames(taxa.out.ori[[type]][[rank]]))[i]
        
        tax.tab[[rank+1]][ind] <- unique(colnames(taxa.out[[type]][[rank]]))[i]
      }
    }
  }
  return(list(taxa.out = taxa.out, tax.tab = tax.tab))
}

remove.NA.out <- function(bin.out, test.out, names) {
  remove.out <- lapply(test.out, function(x) {
    ind <- grep("_NA ", rownames(x))
    x <- x[-ind,]
  })
  
  remove.names <- names
  
  remove.names$names <- lapply(names$names, function(x) {
    ind <- grep("NA ", x)
    x <- x[-ind]
  })
  
  remove.bin <- bin.out
  
  remove.bin$taxa <- lapply(bin.out$taxa, function(x) {
    ind <- grep("_NA ", colnames(x))
    x <- x[,-ind]
  })
  
  return(list(bin.out = remove.bin, test.out = remove.out, names.out = remove.names))
}

########
# Taxa #
########

#################
# Taxa analysis #
#################

##########################
# Binary - No covariates #
##########################

taxa.ind.sum.func <- function(x) {
  sum.out <- c(length(x), mean(x, na.rm = TRUE), quantile(x, na.rm = TRUE))
  names(sum.out) <-  c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(sum.out)
}

taxa.bin.sum.func <- function(bin.var, taxa) {
  
  if(is.null(ncol(taxa))){
    out <- NULL
    
  }else{
    n.taxa <- ncol(taxa)
    ref.sum <- matrix(NA, n.taxa, 7)
    com.sum <- matrix(NA, n.taxa, 7)
    for (i in 1:n.taxa) {
      ind.taxa <- taxa[,i]
      sum.out <- tapply(ind.taxa, bin.var, taxa.ind.sum.func)
      ref.sum[i,] <- sum.out[[1]]
      com.sum[i,] <- sum.out[[2]]
    }
    rownames(ref.sum) <- colnames(taxa)
    colnames(ref.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
    rownames(com.sum) <- colnames(taxa)
    colnames(com.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
    out <- list(as.data.frame(ref.sum), as.data.frame(com.sum))
    names(out) <- levels(bin.var)
  }
  return(out)
}

taxa.bin.sum.united.func <- function(bin.var, taxa.out) {
  taxa.bin.sum <-list()
  for(i in 1:6) {
    taxa.bin.sum[[i]] <- taxa.bin.sum.func(bin.var, taxa.out[[i]])
  }
  names(taxa.bin.sum) <- names(taxa.out)
  return(taxa.bin.sum)
}

# Multiple testing correction

bin.q.func <- function(out, method = c("BH", "BY")) {
  if(is.null(out)){
    return(NULL)
  }else{
    Q.value <- p.adjust(out$P.value, method = method)
    return(cbind(out, Q.value))
  }
}

bin.q.united.func <- function(taxa.out, method = "BH") {
  q.out <- list()
  for(i in 1:6) {
    q.out[[i]] <- bin.q.func(taxa.out[[i]], method)
  }
  names(q.out) <- names(taxa.out)
  return(q.out)
}

bin.t.test.q.fcr.func <- function(bin.var, taxa, method = c("BH", "BY")) {
  n.taxa <- ncol(taxa)
  out <- matrix(NA, n.taxa, 6)
  for (i in 1:n.taxa) {
    fit <- t.test(taxa[,i] ~ bin.var)
    out[i,] <- c(fit$statistic, fit$stderr, fit$parameter, fit$conf.int, fit$p.value)   
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(taxa)
  colnames(out) <- c("t", "Std Err", "DF", "Lower", "Upper", "P.value")
  Q.value <- p.adjust(out$P.value, method = method)
  R <- sum(Q.value < 0.05)
  conf.level <- 1 - R*0.05/n.taxa
  fcr.out <- matrix(NA, n.taxa, 2)
  for (i in 1:n.taxa) {
    fit <- t.test(taxa[,i] ~ bin.var, conf.level = conf.level)
    fcr.out[i,] <- fit$conf.int   
  }
  fcr.out <- as.data.frame(fcr.out)
  rownames(fcr.out) <- colnames(taxa)
  colnames(fcr.out) <- c("FCR.Lower", "FCR.Upper")
  return(cbind(out, fcr.out, Q.value))
}

taxa.bin.cat.ref.logit.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  #bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa <- rbind(taxa[ind.ref,], taxa[ind.com,])
  
  return(list(bin.var = bin.var, taxa = taxa))
}

taxa.bin.cat.ref.logit.united.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  #bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa.out <- list()
  for(i in 1:6) {
    taxon <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
    taxa.out[[i]] <- taxon
  }
  names(taxa.out) <- names(taxa)
  
  return(list(bin.var = bin.var, taxa = taxa.out))
}

taxa.bin.cat.ref.beta.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa){
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  taxa <- rbind(taxa[ind.ref,], taxa[ind.com,])
  
  #taxa <- zCompositions::cmultRepl(taxa)
  
  return(list(bin.var = bin.var, taxa = taxa))
}

taxa.bin.var.func <- function(sam.dat) {
  var.names <- colnames(sam.dat)
  return(var.names)
}

taxa.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  if (length(bin.cat) != 2) {
    stop(paste(sel.bin.var, " is not binary", sep = ""))
  }
  return(bin.cat)
}

taxa.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}

taxa.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  return(sam.dat)
}

taxa.bin.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa <- rbind(taxa[ind.ref,], taxa[ind.com,])
  
  return(list(bin.var = bin.var, taxa = taxa))
}

taxa.bin.cat.ref.united.func <- function(sel.bin.var, sel.ref, sel.com, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa.out <- list()
  for(i in 1:6) {
    taxon <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
    taxa.out[[i]] <- taxon
  }
  names(taxa.out) <- names(taxa)
  
  return(list(bin.var = bin.var, taxa = taxa.out))
}

taxa.bin.cov.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sel.cov.var, sam.dat, taxa) { ## 0527
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  taxa <- rbind(taxa[ind.ref,], taxa[ind.com,])
  
  return(list(bin.var = bin.var, cov.var = cov.var, taxa = taxa))
}

taxa.bin.cov.cat.ref.united.func <- function(sel.bin.var, sel.ref, sel.com, sel.cov.var, sam.dat, taxa) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  taxa.list <- list()
  for (i in 1:6) {
    taxa.rank <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
    taxa.list[[i]] <- taxa.rank
  }
  names(taxa.list) <- names(taxa)
  return(list(bin.var = bin.var, cov.var = cov.var, taxa = taxa.list))
}

get.or.se <- function(model) {
  broom::tidy(model) %>%
    mutate(or = exp(estimate),
           var.diag = diag(vcov(model)),
           or.se = sqrt(or^2 * var.diag)) %>%
    select(or.se) %>% unlist %>% unname
}

get.or.se2 <- function(model){
  or = exp(broom::tidy(model)$estimate)
  var.diag = diag(vcov(model))
  or.se = sqrt(or^2 * var.diag)
  return(unname(unlist(or.se)))
}

taxa.bin.t.test <- function(bin.var, taxa) {
  if(is.null(taxa)){
    out <- NULL
  }else{
    n.taxa <- ncol(taxa)
    out <- matrix(NA, n.taxa, 6)
    for (i in 1:n.taxa) {
      fit <- t.test(taxa[,i] ~ bin.var)
      out[i,] <- c(-fit$statistic, fit$stderr, fit$parameter, -fit$conf.int[2], -fit$conf.int[1], fit$p.value)   
    }
    out <- as.data.frame(out)
    rownames(out) <- colnames(taxa)
    colnames(out) <- c("t", "Std Err", "DF", "Lower", "Upper", "P.value")
    
  }
  return(out)
}

taxa.bin.t.test.united <- function(bin.var, taxa) {
  t.test <- list()
  for(i in 1:6) {
    taxon <- taxa[[i]]
    t.test[[i]]<-taxa.bin.t.test(bin.var, taxon)
  }
  names(t.test) <- names(taxa)
  return(t.test)
}

taxa.bin.wilcox.test <- function(bin.var, taxa) {
  if(is.null(taxa)){
    out <- NULL
  }else{
    n.taxa <- ncol(taxa)
    out <- matrix(NA, n.taxa, 2)
    for (i in 1:n.taxa) {
      fit <- wilcox.test(taxa[,i] ~ bin.var)
      out[i,] <- c(fit$statistic, fit$p.value)   
    }
    out <- as.data.frame(out)
    rownames(out) <- colnames(taxa)
    colnames(out) <- c("W", "P.value")
  }
  return(out)
}

taxa.bin.wilcox.test.united <- function(bin.var, taxa) {
  wilcox.test <- list()
  for (i in 1:6) {
    taxon <- taxa[[i]]
    wilcox.test[[i]] <- taxa.bin.wilcox.test(bin.var, taxon)
  }
  names(wilcox.test) <- names(taxa)
  return(wilcox.test)
}

taxa.wilcox.test.est.func <- function(bin.var, taxa, sel.ref, sel.com, q.out) {
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  
  taxa.added <- list()
  for(i in 1:6) {
    if(!is.null(taxa[[i]])){
      n.tax <- length(taxa[[i]])
      med.diff <- numeric()
      for(j in 1:n.tax) {  
        taxon <- taxa[[i]][,j]
        med.ref <- median(taxon[ind.ref], na.rm = TRUE)
        med.com <- median(taxon[ind.com], na.rm = TRUE)
        med.diff[j] <- med.com - med.ref
      }
      q.out[[i]] <- cbind(Est = med.diff, q.out[[i]])
    }
  }
  return(q.out)
}

all.taxa.bin.t.test <- function(taxa.out, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    clr.taxa <- taxa.out$clr[[i]]
    clr.taxa.out <- taxa.bin.cat.ref.func(sel.bin.var = "ecig_status", sel.ref = "Non", sel.com = "E-cig", sam.dat = sam.dat, taxa = clr.taxa)
    clr.taxa.sum.out <- taxa.bin.sum.func(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.out <- taxa.bin.t.test(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.q.out <- bin.q.func(clr.taxa.test.out, method = multi.method)
    tax.out[[i]] <- list(clr.sum = clr.taxa.sum.out, clr.test.out = clr.taxa.test.q.out)  
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.bin.wilcox.test <- function(taxa.out, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    clr.taxa <- taxa.out$clr[[i]]
    clr.taxa.out <- taxa.bin.cat.ref.func(sel.bin.var = "ecig_status", sel.ref = "Non", sel.com = "E-cig", sam.dat = sam.dat, taxa = clr.taxa)
    clr.taxa.sum.out <- taxa.bin.sum.func(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.out <- taxa.bin.wilcox.test(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.q.out <- bin.q.func(clr.taxa.test.out, method = multi.method)
    tax.out[[i]] <- list(clr.sum = clr.taxa.sum.out, clr.test.out = clr.taxa.test.q.out)  
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.bin.t.test.united <- function(sel.bin.var, sel.ref, sel.com, taxa.out, sam.dat, multi.method = "BH") {
  tax.out <- list()
  
  for (i in 1:6) {
    clr.taxa <- taxa.out$clr[[i]]
    clr.taxa.out <- taxa.bin.cat.ref.func(sel.bin.var, sel.ref, sel.com, sam.dat = sam.dat, taxa = clr.taxa)
    
    clr.taxa.sum.out <- taxa.bin.sum.func(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.out <- taxa.bin.t.test(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.q.out <- bin.q.func(clr.taxa.test.out, method = multi.method)
    
    rank <- names(taxa.out[[1]])[i]
    rank <- paste(toupper(substr(rank,1,1)),substr(rank,2,nchar(rank)), sep = "")
    clr.taxa.sum.out[[1]]<-cbind(Rank = rank, clr.taxa.sum.out[[1]])
    clr.taxa.sum.out[[2]]<-cbind(Rank = rank, clr.taxa.sum.out[[2]])
    clr.taxa.test.q.out<-cbind(Rank = rank, clr.taxa.test.q.out)
    if (i==1) {
      tax.out[[1]] <- clr.taxa.sum.out
      tax.out[[2]] <- clr.taxa.test.q.out
    }
    else {
      tax.out[[1]] <- Map(rbind,tax.out[[1]],clr.taxa.sum.out)
      tax.out[[2]] <- rbind(tax.out[[2]],clr.taxa.test.q.out)
    }
  }
  names(tax.out) <- c("clr.sum","clr.test.out")
  return(tax.out)
}

all.taxa.bin.wilcox.test.united <- function(taxa.out, sam.dat, multi.method = "BH") {
  tax.out <- list()
  
  for (i in 1:6) {
    clr.taxa <- taxa.out$clr[[i]]
    clr.taxa.out <- taxa.bin.cat.ref.func(sel.bin.var = "ecig_status", sel.ref = "Non", sel.com = "E-cig", sam.dat = sam.dat, taxa = clr.taxa)

    clr.taxa.sum.out <- taxa.bin.sum.func(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.out <- taxa.bin.wilcox.test(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.test.q.out <- bin.q.func(clr.taxa.test.out, method = "BH")
    
    rank <- names(taxa.out[[1]])[i]
    rank <- paste(toupper(substr(rank,1,1)),substr(rank,2,nchar(rank)), sep = "")

    clr.taxa.sum.out[[1]]<-cbind(Rank = rank, clr.taxa.sum.out[[1]])
    clr.taxa.sum.out[[2]]<-cbind(Rank = rank, clr.taxa.sum.out[[2]])
    clr.taxa.test.q.out<-cbind(Rank = rank, clr.taxa.test.q.out)
    if (i==1) {
      tax.out[[1]] <- clr.taxa.sum.out
      tax.out[[2]] <- clr.taxa.test.q.out
    }
    else {
      tax.out[[1]] <- Map(rbind,tax.out[[1]],clr.taxa.sum.out)
      tax.out[[2]] <- rbind(tax.out[[2]],clr.taxa.test.q.out)
    }
  }
  names(tax.out) <- c("clr.sum","clr.test.out")
  return(tax.out)
}

taxa.bin.lm.func <- function(bin.var, taxa) {  # without covariate
  
  n.tax <- ncol(taxa)
  lm.out <- matrix(NA, n.tax, 6)
  for (i in 1:n.tax) {
    taxon <- taxa[,i]
    fit <- try(lm(taxon ~ bin.var), silent = TRUE)
    
    est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
    std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
    df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
    if(is.na(df)){
      ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) C(NA,NA))
    }else{
      ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) C(NA,NA))
    }
    
    p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
    
    out <- c(est, std.err, df, ci, p.val)
    
    #out.lm <- c(summary(fit.lm)$coefficients[2,c(1,2)], summary(fit.lm)$df[2], confint(fit.lm)[2,], summary(fit.lm)$coefficients[2,4])
    lm.out[i,] <- out
  }
  #here3
  lm.out <- as.data.frame(lm.out)
  rownames(lm.out) <- colnames(taxa)
  colnames(lm.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(lm.out)
}

taxa.bin.lm.united.func <- function(bin.var, taxa) {
  lm.test <- list()
  for(i in 1:6) {
    taxon <- taxa[[i]]
    lm.test[[i]] <- taxa.bin.lm.func(bin.var, taxon)
  }
  names(lm.test) <- names(taxa)
  return(lm.test)
}

all.taxa.bin.glm.nb <- function(bin.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    count.taxa.test.out <- taxa.bin.glm.nb.func(bin.var, taxa[[i]], library.size)
    count.taxa.test.q.out <- bin.q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.glm.nb.func <- function(bin.var, taxa, library.size) {
  n.tax <- ncol(taxa)
  lmer.out <- matrix(NA, n.tax, 6)
  for (i in 1:n.tax) {
    taxon <- taxa[,i]
    dat <- as.data.frame(cbind(bin.var, taxon))
    f <- formula(paste(colnames(dat)[2], " ~ ", colnames(dat)[1], "+ offset(log(library.size))", sep = ""))
    
    fit <- try(glm.nb(f, data = dat), silent = TRUE)
    est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
    std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
    DF <- NA
    ci <- c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err)
    p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
    out <- c(est, std.err, DF, ci, p.val)
    lmer.out[i,] <- out
  }
  lmer.out <- as.data.frame(lmer.out)
  rownames(lmer.out) <- colnames(taxa)
  colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(lmer.out)
}

# taxa.bin.logit.func <- function(bin.var, taxa){
#   n.tax <- ncol(taxa)
#   log.out <- matrix(NA, n.tax, 6)
#   for (i in 1:n.tax) {
#     taxon <- taxa[,i]
#     fit <- try(glm(unlist(bin.var) ~ taxon, family = binomial()), silent = TRUE)
#     #here1
#     
#     est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
#     std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
#     df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
#     ci <- c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err)
#     p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
#     
#     out <- c(est, std.err, df, ci, p.val)
#     #out.log <- c(summary(fit.log)$coefficients[2,c(1,2)], summary(fit.log)$df[2], confint(fit.log)[2,], summary(fit.log)$coefficients[2,4])
#     log.out[i,] <- out
#   }
#   log.out <- as.data.frame(log.out)
#   rownames(log.out) <- colnames(taxa)
#   colnames(log.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
#   return(log.out)
# }

all.taxa.logit.reg.coef.bin.func <- function(bin.var, taxa.out, scale = TRUE) {
  all.logit.out <- list()
  bin.var <- unlist(bin.var)
  for(i in 1:6) {
    if(!is.null(taxa.out[[i]])){
      taxa <- taxa.out[[i]]
      n.tax <- ncol(taxa)
      logit.out <- matrix(NA, n.tax, 6)
      
      for (j in 1:n.tax) {
        taxon <- taxa[,j]
        if(scale == TRUE){
          fit <- try(glm(bin.var ~ scale(taxon), family = "binomial"), silent = TRUE)
        } else {
          fit <- try(glm(bin.var ~ taxon, family = "binomial"), silent = TRUE)
        }
        #here2
        est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
        std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
        df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
        if(is.na(df)){
          ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) C(NA,NA))
        }else{
          ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) C(NA,NA))
        }
        
        p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
        
        out <- c(est, std.err, df, ci, p.val)
        
        logit.out[j,] <- out
      }
      logit.out <- as.data.frame(logit.out)
      rownames(logit.out) <- colnames(taxa)
      colnames(logit.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
      logit.q.out <- bin.q.func
      all.logit.out[[i]] <-logit.out
    }
    
  }
  names(all.logit.out) <- names(taxa.out)
  return(all.logit.out)
}

all.taxa.logit.bin.func <- function(bin.var, taxa.out, scale = TRUE) {
  all.logit.out <- list()
  for(i in 1:6) {
    if(!is.null(taxa.out[[i]])){
      taxa <- taxa.out[[i]]
      n.tax <- ncol(taxa)
      logit.out <- matrix(NA, n.tax, 6)
      for (j in 1:n.tax) {
        taxon <- taxa[,j]
        d <- as.data.frame(cbind(bin.var, taxon))
        #d[,1] <- as.numeric(d[,1])
        #d[,2] <- as.numeric(d[,2])
        if(scale){
          logit.f <- formula(paste(colnames(bin.var), "~" , "scale(taxon)"))
        } else {
          logit.f <- formula(paste(colnames(bin.var), "~" , "taxon"))
        }
        #here7
        fit <- try(glm(logit.f, data = d , family = "binomial"), silent = TRUE)

        or <- tryCatch(exp(summary(fit)$coefficients[2,1]), error = function(err) NA)
        or.se <- tryCatch(sqrt(or^2*diag(vcov(fit)))[2], error = function(err) NA)
        df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
        if(is.na(df)){
          ci <- tryCatch(c(exp(summary(fit)$coefficients[2,1] - qnorm(0.975)*summary(fit)$coefficients[2,2]), exp(summary(fit)$coefficients[2,1] + qnorm(0.975)*summary(fit)$coefficients[2,2])), error = function(err) C(NA,NA))
        }else{
          ci <- tryCatch(c(exp(summary(fit)$coefficients[2,1] - qt(0.975, df = df)*summary(fit)$coefficients[2,2]), exp(summary(fit)$coefficients[2,1] + qt(0.975, df = df)*summary(fit)$coefficients[2,2])), error = function(err) C(NA,NA))
        }
        p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
        out <- c(or, or.se, df, ci, p.val)
        logit.out[j,] <- out
      }
      logit.out <- as.data.frame(logit.out)
      rownames(logit.out) <- colnames(taxa)
      colnames(logit.out) <- c("OR", "Std Err", "DF", "Lower", "Upper", "P.value")
      logit.q.out <- bin.q.func
      all.logit.out[[i]] <-logit.out
    }
  }
  names(all.logit.out) <- names(taxa.out)
  return(all.logit.out)
}

all.taxa.bin.beta <- function(bin.var, taxa, multi.method = "BH"){
  
  tax.out <- list()
  for (i in 1:6) {
    #prop.taxa.out <- taxa.bin.cat.ref.beta.func(sel.bin.var = sel.bin.var, sel.ref = sel.ref, sel.com = sel.com, sam.dat = sam.dat, taxa[[i]])
    prop.taxa.test.out <- taxa.bin.beta.func(bin.var = bin.var, taxa = taxa[[i]])
    prop.taxa.test.q.out <- bin.q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
  
}

taxa.bin.beta.func <- function(bin.var, taxa) {
  n.tax <- ncol(taxa)
  beta.out <- matrix(NA, n.tax, 6)
  for (i in 1:n.tax) {
    taxon <- taxa[,i]
    dat <- as.data.frame(cbind(bin.var, taxon))
    dat[,2] <- as.numeric(dat[,2])
    dat[,1] <- as.numeric(dat[,1])
    f <- formula(paste(colnames(dat)[2], " ~ ", colnames(dat)[1], sep = ""))
    fit.beta <- try(betareg(f, data = dat, maxit = 50), silent = TRUE)
    
    if(class(fit.beta) != "try-error"){
      est <- summary(fit.beta)$coefficients$mean[2,1]
      std.err <- summary(fit.beta)$coefficients$mean[2,2]
      ci <- c(est + qnorm(c(0.025, 0.975))*std.err)
      out.beta <- c(summary(fit.beta)$coefficients$mean[2,c(1,2)], NA, ci, summary(fit.beta)$coefficients$mean[2,4])
    } else {
      out.beta <- rep(NA, 6)
    }
    beta.out[i,] <- out.beta
  }
  beta.out <- as.data.frame(beta.out)
  rownames(beta.out) <- colnames(taxa)
  colnames(beta.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(beta.out)
}

#######################
# Binary - Covariates #
#######################

taxa.bin.cov.lm.united.func <- function(bin.var, cov.var, taxa) {
  n.cov <- ncol(cov.var)
  lm.out.list <- list()
  for(i in 1:6) {
    if(!is.null(taxa[[i]])){
      n.tax <- ncol(taxa[[i]])
      d <- as.data.frame(cbind(bin.var, cov.var, taxa[[i]]))
      lm.out <- matrix(NA, n.tax, 6)
      for (j in 1:n.tax) {
        taxon <- taxa[[i]][,j]
        lm.f <- formula(paste("taxon", "~", colnames(bin.var), "+", paste(colnames(cov.var), collapse = "+")))
        fit <- try(lm(lm.f, data = d), silent = TRUE)
        est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
        std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
        df <- tryCatch(summary(fit)$df[2], error = function(err) NA)

        if(is.na(df)){
          ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) C(NA,NA))
        }else{
          ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) C(NA,NA))
        }
        p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
        
        out <- c(est, std.err, df, ci, p.val)
        lm.out[j,] <- out
      }
      lm.out <- as.data.frame(lm.out)
      rownames(lm.out) <- colnames(taxa[[i]])
      colnames(lm.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
      lm.out.list[[i]] <- lm.out
    }
  }
  names(lm.out.list) <- names(taxa)
  return(lm.out.list)
}

all.taxa.lm.bin.united <- function(taxa.out, sam.dat, multi.method = "BH") { # + cov
  
  tax.out <- list()
  
  for (i in 1:6) {
    clr.taxa <- taxa.out$clr[[i]]
    clr.taxa.out <- taxa.bin.cat.ref.func(sel.bin.var = "ecig_status", sel.ref = "Non", sel.com = "E-cig", sam.dat = sam.dat, taxa = clr.taxa)
    
    clr.taxa.sum.out <- taxa.bin.sum.func(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.lm.out <- taxa.bin.lm.func(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.lm.q.out <- bin.q.func(clr.taxa.lm.out, method = "BH")
    
    rank <- names(taxa.out[[1]])[i]
    rank <- paste(toupper(substr(rank,1,1)),substr(rank,2,nchar(rank)), sep = "")
    clr.taxa.sum.out[[1]]<-cbind(Rank = rank, clr.taxa.sum.out[[1]])
    clr.taxa.sum.out[[2]]<-cbind(Rank = rank, clr.taxa.sum.out[[2]])
    clr.taxa.lm.q.out<-cbind(Rank = rank, clr.taxa.lm.q.out)
    if (i==1) {
      tax.out[[1]] <- clr.taxa.sum.out
      tax.out[[2]] <- clr.taxa.lm.q.out
    }
    else {
      tax.out[[1]] <- Map(rbind,tax.out[[1]],clr.taxa.sum.out)
      tax.out[[2]] <- rbind(tax.out[[2]],clr.taxa.lm.q.out)
    }
  }
  names(tax.out) <- c("clr.sum","clr.test.out")
  return(tax.out)
}

all.taxa.lm.bin.cov.united <- function(sel.bin.var, sel.ref, sel.com, sel.cov.var,sam.dat, taxa.out, multi.method = "BH") { # + cov
  
  tax.out <- list()
  
  for (i in 1:6) {
    clr.taxa <- taxa.out[[i]]
    clr.taxa.out <- taxa.bin.cov.cat.ref.func(sel.bin.var, sel.ref, sel.com, sel.cov.var , sam.dat, clr.taxa)
    clr.taxa.sum.out <- taxa.bin.sum.func(bin.var = clr.taxa.out$bin.var, taxa = clr.taxa.out$taxa)
    clr.taxa.lm.out <- taxa.bin.cov.lm.func(bin.var = clr.taxa.out$bin.var, cov.var = clr.taxa.out$cov.var, taxa = clr.taxa.out$taxa)
    clr.taxa.lm.q.out <- bin.q.func(clr.taxa.lm.out, method = multi.method)
    
    rank <- names(taxa.out)[i]
    rank <- paste(toupper(substr(rank,1,1)),substr(rank,2,nchar(rank)), sep = "")
    clr.taxa.sum.out[[1]]<-cbind(Rank = rank, clr.taxa.sum.out[[1]])
    clr.taxa.sum.out[[2]]<-cbind(Rank = rank, clr.taxa.sum.out[[2]])
    clr.taxa.lm.q.out<-cbind(Rank = rank, clr.taxa.lm.q.out)
    if (i==1) {
      tax.out[[1]] <- clr.taxa.sum.out
      tax.out[[2]] <- clr.taxa.lm.q.out
    }
    else {
      tax.out[[1]] <- Map(rbind,tax.out[[1]],clr.taxa.sum.out)
      tax.out[[2]] <- rbind(tax.out[[2]],clr.taxa.lm.q.out)
    }
  }
  names(tax.out) <- c("clr.sum","clr.test.out")
  return(tax.out)
}

all.taxa.logit.reg.coef.bin.cov.func <- function(bin.var, cov.var, taxa.out, scale = TRUE) {
  all.logit.out <- list()
  for(i in 1:6) {
    if(!is.null(taxa.out[[i]])){
      taxa <- taxa.out[[i]]
      n.tax <- ncol(taxa)
      logit.out <- matrix(NA, n.tax, 6)
      for (j in 1:n.tax) {
        taxon <- taxa[,j]
        d <- as.data.frame(cbind(bin.var, cov.var, taxon))

        if(scale){
          logit.f <- formula(paste(colnames(bin.var), "~" , "scale(taxon)", "+", paste(colnames(cov.var), collapse = "+")))
        }else {
          logit.f <- formula(paste(colnames(bin.var), "~" , "taxon", "+", paste(colnames(cov.var), collapse = "+")))
        }

        fit <- try(glm(logit.f, data = d, family = "binomial"), silent = TRUE)
        
        est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
        std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
        df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
        if(is.na(df)){
          ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) C(NA,NA))
        }else{
          ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) C(NA,NA))
        }
        p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
        out <- c(est, std.err, df, ci, p.val)
        
        logit.out[j,] <- out
      }
      logit.out <- as.data.frame(logit.out)
      rownames(logit.out) <- colnames(taxa)
      colnames(logit.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
      logit.q.out <- bin.q.func
      all.logit.out[[i]] <-logit.out
    }
  }
  names(all.logit.out) <- names(taxa.out)
  return(all.logit.out)
}

all.taxa.logit.bin.cov.func <- function(bin.var, cov.var, taxa.out, scale = TRUE) {
  all.logit.out <- list()
  for(i in 1:6) {
    if(!is.null(taxa.out[[i]])){
      taxa <- taxa.out[[i]]
      n.tax <- ncol(taxa)
      logit.out <- matrix(NA, n.tax, 6)
      for (j in 1:n.tax) {
        taxon <- taxa[,j]
        d <- as.data.frame(cbind(bin.var, cov.var, taxon))
        if(scale){
          logit.f <- formula(paste(colnames(bin.var), "~" , "scale(taxon)", "+", paste(colnames(cov.var), collapse = "+")))
        } else {
          logit.f <- formula(paste(colnames(bin.var), "~" , "taxon", "+", paste(colnames(cov.var), collapse = "+")))
        }
        fit <- try(glm(logit.f, data = d , family = "binomial"), silent = TRUE)
        
        or <- tryCatch(exp(summary(fit)$coefficients[2,1]), error = function(err) NA)
        or.se <- tryCatch(sqrt(or^2*diag(vcov(fit)))[2], error = function(err) NA)
        df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
        if(is.na(df)){
          ci <- tryCatch(c(exp(summary(fit)$coefficients[2,1] - qnorm(0.975)*summary(fit)$coefficients[2,2]), exp(summary(fit)$coefficients[2,1] + qnorm(0.975)*summary(fit)$coefficients[2,2])), error = function(err) c(NA,NA))
        }else{
          ci <- tryCatch(c(exp(summary(fit)$coefficients[2,1] - qt(0.975, df = df)*summary(fit)$coefficients[2,2]), exp(summary(fit)$coefficients[2,1] + qt(0.975, df = df)*summary(fit)$coefficients[2,2])), error = function(err) c(NA,NA))
        }
        p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
        out <- c(or, or.se, df, ci, p.val)
        logit.out[j,] <- out
      }
      logit.out <- as.data.frame(logit.out)
      rownames(logit.out) <- colnames(taxa)
      colnames(logit.out) <- c("OR", "Std Err", "DF", "Lower", "Upper", "P.value")
      logit.q.out <- bin.q.func
      all.logit.out[[i]] <-logit.out
    }
  }
  names(all.logit.out) <- names(taxa.out)
  return(all.logit.out)
}

taxa.bin.cov.glm.nb.func <- function(bin.var, cov.var, taxa, library.size) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, cov.var, taxon))
      f <- formula(paste(colnames(dat)[ncol(dat)], " ~ ", colnames(dat)[1], "+", paste(names(cov.var), collapse = "+"), "+ offset(log(library.size))", sep = ""))
      fit <- try(glm.nb(f, data = dat), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
      df <- NA
      ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) C(NA,NA))
      p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
      out <- c(est, std.err, df, ci, p.val)
      lmer.out[i,] <- out
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.bin.cov.glm.nb <- function(bin.var, cov.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  
  for (i in 1:6) {
    count.taxa.test.out <- taxa.bin.cov.glm.nb.func(bin.var, cov.var, taxa[[i]], library.size)
    count.taxa.test.q.out <- bin.q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <- count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.bin.cov.beta <- function(bin.var, cov.var, taxa, multi.method = "BH"){
  tax.out <- list()
  for (i in 1:6) {
    prop.taxa.test.out <- taxa.bin.cov.beta.func(bin.var, cov.var, taxa[[i]])
    prop.taxa.test.q.out <- bin.q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.cov.beta.func <- function(bin.var, cov.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    beta.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, cov.var, taxon))
      f <- formula(paste(colnames(dat)[ncol(dat)], " ~ ", colnames(dat)[1], "+", paste(names(cov.var), collapse = "+"), sep = ""))
      fit.beta <- try(betareg(f, data = dat, maxit = 50), silent = TRUE)
      
      if(class(fit.beta) != "try-error"){
        est <- tryCatch(summary(fit.beta)$coefficients$mean[2,1], error = function(err) NA)
        std.err <- tryCatch(summary(fit.beta)$coefficients$mean[2,2], error = function(err) NA)
        ci <- tryCatch(c(est + qnorm(c(0.025, 0.975))*std.err), error = function(err) NA)
        p.val <- tryCatch(summary(fit.beta)$coefficients$mean[2,4], error = function(err) NA)
        out.beta <- c(est, std.err, NA, ci, p.val)
        
      } else {
        out.beta <- rep(NA, 6)
      }
      
      beta.out[i,] <- out.beta
    }
    beta.out <- as.data.frame(beta.out)
    rownames(beta.out) <- colnames(taxa)
    colnames(beta.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(beta.out)
  }
}

taxa.forest.plot.pages <- function(all.taxa.q.out, species.include, mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out$Q.value < 0.05)  
    } else {
      ind.sig <- which(out$P.value < 0.05)
    }
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    plot.per.page <- total
    num.pages <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
    }
  } else if(total<=80 & total>=0) {
    num.pages = 1
    plot.per.page <-total
  }
  
  return(num.pages)
}

taxa.forest.plot.pages1 <- function(all.taxa.q.out, taxa.names.out, species.include, report.type = "Est", mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  if(report.type == "Est"){
    report.txt <- "Est."
  } else if(report.type == "OR") {
    report.txt <- "OR"
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out$Q.value < 0.05)  
    } else {
      ind.sig <- which(out$P.value < 0.05)
    }
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- mod(total, num.pages)
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  all.text.tab <- list()
  all.ci.tab <- list()
  if(total >0) {
    plot.taxa <- list()
    rank <- 1
    tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
    current.page <- 1
    text.tab.all <- 0
    text.tab.save <- numeric()
    ci.tab.save <- numeric()
    info.taxa <- numeric()
    info.data <- numeric()
    info.rank <- numeric()
    info.ci <- numeric()
    sig.out <- numeric()
    actual.plot.per.page <- plot.per.page
    
    if (text.tab.all == 0) {
      if (mult.test.cor){
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value", "Q-value"), 1, 6)  
      } else {
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value"), 1, 5)
      }
      ci.tab.all <- matrix(c(NA, NA, NA), 1, 3)
    }
    
    for (rank in 1:range) {
      if(!is.null(all.taxa.q.out[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- all.taxa.q.out[[rank]][sig.by.rank[[rank]],]
        
        info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
        info.data <- rbind(info.data, sig.out)
        info.ci <- rbind(info.ci, cbind(sig.out[,report.type], sig.out[,"Lower"], sig.out[,"Upper"]))
        info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
      }
    }
    info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.type, "P.value", "Q.value")],3), nsmall = 3)))
    
    for(p in 1:num.pages) {
      
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      names(info) <- names(text.tab.all)
      
      all.text.tab[[p]] <- rbind(as.matrix(text.tab.all), info[c(initial:(initial+actual.plot.per.page-1)),])
      all.ci.tab[[p]] <- rbind(ci.tab.all, info.ci[c(initial:(initial+actual.plot.per.page-1)),])
    }
  } else {
    all.text.tab <- NULL
    all.ci.tab <- NULL
  }
  return(list(all.text.tab = all.text.tab, all.ci.tab = all.ci.tab))
}

taxa.forest.plot.pages2 <- function(page.taxa.q.out, page) {
  
  text.tab.all <- page.taxa.q.out$all.text.tab[[page]]
  ci.tab.all <- page.taxa.q.out$all.ci.tab[[page]]
  
  if(is.null(text.tab.all) & is.null(ci.tab.all)){
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  }else{
    str.max <- lapply(page.taxa.q.out$all.text.tab, nchar)
    for(i in 1:length(page.taxa.q.out$all.text.tab)){
      str.max[[i]] <- str.max[[i]][,3]
    }
    maxStr <- max(unlist(str.max))
    if(!is.numeric(maxStr)){
      maxStr <- 0
    }
    
    text.tab.all[,3] <- substr(text.tab.all[,3], 1, 55)
    #par(mar=c(0, 0.2, 0, 0.2))
    if(text.tab.all[1,4] == "Est."){
      if(nrow(ci.tab.all) <= 10) {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = unit(0.85, "cm"), #line.margin = unit(0.2, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      else if(nrow(ci.tab.all) <= 45){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), 
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      
    } else {
      if(nrow(ci.tab.all) <= 10) {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = unit(0.85, "cm"), #line.margin = unit(0.2, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      else if(nrow(ci.tab.all) <= 45){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), 
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      
      #forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
      #           zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = unit(0.4,"cm"),
      #           col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",# mar = unit(c(0.5,0,0.5,0), "cm"),
      #           txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
      #                          ticks=gpar(fontfamily="", cex=0.75),
      #                          xlab=gpar(fontfamily="", cex=0.75)))
    } 
  }
}

# duplicate.list <- function(taxa.names.rank.out, species.include = FALSE){
#   if(species.include){
#     display <- 5
#   }else{
#     display <- 6
#   }
#   if(sum(!is.na(unlist(taxa.names.rank.out$duplicates))) > 0){
#     duplicate.taxa <- unlist(taxa.names.rank.out$duplicates[1:display])[!is.na(unlist(taxa.names.rank.out$duplicates[1:display]))]
#     par(mar=c(0, 0.5, 0, 0.5))
#     text(x=0, y=0.5, paste(duplicate.taxa, collapse = "\n"), cex = 0.75, adj = c(0, NA))
#   } else {
#     text(x=0.5, y=0.5, "")
#   }
# }

duplicate.list <- function(duplicate.taxa, taxon.inplot, duplicate.full.list){
  if(length(duplicate.taxa %in% taxon.inplot)>0) {
    duplicate.taxa <- unlist(duplicate.full.list)[duplicate.taxa %in% taxon.inplot]
    par(mar=c(0, 0.5, 0, 0.5))
    text(x=0, y=0.5, paste(duplicate.taxa, collapse = "\n"), cex = 0.75, adj = c(0, NA))
  } else {
    print("?!")
    text(x=0.5, y=0.5, "")
  }
}

#############################
# Continuous - No covariate #
#############################

taxa.con.recode.func <- function(sam.dat, sel.con.var, rename.con.var, taxa) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  taxa.dat <- list()
  for (i in 1:6) {
    taxon <- taxa[[i]]
    taxon <- taxon[ind.nona,]
    taxa.dat[[i]] <- taxon
  }
  names(taxa.dat) <- names(taxa)
  
  return(list(con.var = con.var, taxa = taxa.dat))
}

taxa.con.beta.recode.func <- function(sam.dat, sel.con.var, rename.con.var, taxa) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  taxa.dat <- list()
  for (i in 1:6) {
    taxon <- taxa[[i]]
    taxon <- taxon[ind.nona,]
    taxa.dat[[i]] <- taxon
    #taxa.dat[[i]] <- zCompositions::cmultRepl(taxon)
  }
  names(taxa.dat) <- names(taxa)
  
  return(list(con.var = con.var, taxa = taxa.dat))
}

taxa.sum.apply <- function(taxa.out, margin, func) {
  taxa.sum.out <- list()
  for(i in 1:6){
    if(!is.null(taxa.out$taxa[[i]])){
      taxa.sum.out[[i]] <- apply(taxa.out$taxa[[i]], margin, func)
    }
  }
  names(taxa.sum.out) <- names(taxa.out$taxa)
  return(taxa.sum.out)
}

taxa.lm.con.func <- function(con.var, taxa) {
  taxa.lm.out <- list()
  for (j in 1:6) {
    if(!is.null(taxa[[j]])){
      taxa.rank <- taxa[[j]]
      n.tax <- ncol(taxa.rank)
      lm.out <- matrix(NA, n.tax, 6)
      for (i in 1:n.tax) {
        taxon <- taxa.rank[,i]
        fit <- try(lm(taxon ~ con.var[,1]), silent = TRUE)
        
        est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
        std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
        df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
        if(is.na(df)){
          ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) C(NA,NA))
        }else{
          ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) C(NA,NA))
        }
        p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
        out <- c(est, std.err, df, ci, p.val)

        lm.out[i,] <- out
      }
      lm.out <- as.data.frame(lm.out)
      rownames(lm.out) <- colnames(taxa.rank)
      colnames(lm.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
      taxa.lm.out[[j]] <- lm.out
    }
  }
  names(taxa.lm.out) <- names(taxa)
  return(taxa.lm.out)
}

taxa.con.glm.nb.func <- function(con.var, taxa, library.size) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(con.var, taxon))
      dat[,2] <- as.numeric(dat[,2])
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[2], " ~ ", colnames(dat)[1], "+ offset(log(library.size))", sep = ""))
      fit <- try(glm.nb(f, data = dat), silent = TRUE)
      
      est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
      DF <- NA
      ci <- c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err)
      p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
      out <- c(est, std.err, DF, ci, p.val)
      lmer.out[i,] <- out
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.con.glm.nb <- function(con.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    count.taxa.test.out <- taxa.con.glm.nb.func(con.var = con.var, taxa = taxa[[i]], library.size)  
    count.taxa.test.q.out <- bin.q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <- count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.con.beta <- function(con.var, taxa, multi.method = "BH"){
  tax.out <- list()
  for (i in 1:6) {
    prop.taxa.test.out <- taxa.con.beta.func(con.var, taxa[[i]])
    prop.taxa.test.q.out <- bin.q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.con.beta.func <- function(con.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    beta.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(con.var, taxon))
      dat[,2] <- as.numeric(dat[,2])
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[2], " ~ ", colnames(dat)[1], sep = ""))
      fit.beta <- try(betareg(f, data = dat, maxit = 50), silent = TRUE)
      if(class(fit.beta) != "try-error"){
        est <- summary(fit.beta)$coefficients$mean[2,1]
        std.err <- summary(fit.beta)$coefficients$mean[2,2]
        ci <- c(est + qnorm(c(0.025, 0.975))*std.err)
        out.beta <- c(summary(fit.beta)$coefficients$mean[2,c(1,2)], NA, ci, summary(fit.beta)$coefficients$mean[2,4])
      } else {
        out.beta <- rep(NA, 6)
      }
      beta.out[i,] <- out.beta
    }
    beta.out <- as.data.frame(beta.out)
    rownames(beta.out) <- colnames(taxa)
    colnames(beta.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(beta.out)
  }
}

taxa.bin.boxplot <- function(bin.var, taxa.out, t.test.q.out, taxa.names.out, page, mult.test.cor = TRUE) {  ## the taxa.out is actually taxa.out$clr
  
  sig.list <- list()
  
  if(mult.test.cor) {
    for(i in 1:6) {
      sig.list[[i]] <- which(t.test.q.out[[i]]$Q.value < 0.05)
    }
  } else {
    for(i in 1:6) {
      sig.list[[i]] <- which(t.test.q.out[[i]]$P.value < 0.05)
    }
  }
  
  j <- page
  
  nrow <- ceiling(length(sig.list[[j]])/4)
  if(nrow > 0){
    nrow <- ceiling(length(sig.list[[j]])/4)
    par(mfrow = c(nrow,4))
    id <- 0
    if(page > 1) {
      for(m in 1:(page-1)){
        id <- id+length(sig.list[[m]])
      }
    }
    for(k in 1:length(sig.list[[j]])) {
      if (grepl(":",taxa.names.out$names[[j]][sig.list[[j]][k]])) {
        name<- gsub(":","_", taxa.names.out$names[[j]][sig.list[[j]][k]])
      } else {
        name <- taxa.names.out$names[[j]][sig.list[[j]][k]]
      }
      if(mult.test.cor){
        qval.clr <- t.test.q.out[[j]][sig.list[[j]][k],"Q.value"]
        round.qval.clr <- format(round(qval.clr, digits = 3), nsmall = 3)
        xlab.v = paste("*q:", round.qval.clr, sep="")
      } else {
        if(!mult.test.cor){
          pval.clr <- t.test.q.out[[j]][sig.list[[j]][k],"P.value"]
          round.pval.clr <- format(round(pval.clr, digits = 3), nsmall = 3)
          xlab.v = paste("*p:", round.pval.clr, sep="")
        }
      }
      id <- id + 1
      taxon <- as.matrix(taxa.out[[j]][sig.list[[j]][k]])
      boxplot(taxon ~ bin.var, main = id, xlab=xlab.v, ylab=name, names = substr(levels(bin.var), 1, 8), notch = TRUE, col=c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)))
      
    }
  } else {
    ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    plot.new()
    text(x = 0.5, y = 0.5, paste("No significant taxa are found in ", ranks[page], sep = ""), 
         cex = 1.2, col = "black")
  }
  
}

taxa.q.func <- function(out, method = c("BH", "BY")) {
  q.out <- list()
  for(i in 1:6) {
    if(!is.null(out[[i]])){
      Q.value <- p.adjust(out[[i]]$P.value, method = method)
      q.out[[i]] <- cbind(out[[i]],Q.value)
    }
  }
  names(q.out) <- names(out)
  return(q.out)
}

taxa.qr.con.func <- function(con.var, taxa, tau, qr.method, se.method, n.res) {
  taxa.qr.out <- list()
  for (i in 1:6) {
    taxa.rank <- taxa[[i]]
    n.tax <- ncol(taxa.rank)
    qr.out <- matrix(NA, n.tax, 6)
    for(j in 1:n.tax) {
      taxon <- taxa.rank[,j]
      fit.qr <- rq(taxon ~ con.var[,1], tau = tau, method = qr.method)
      if (se.method == "boot") {
        set.seed(j)
      }
      out.coef <- summary(fit.qr, se = se.method, R = n.res)$coefficients
      out.df <- summary(fit.qr, se = se.method, R = n.res)$rdf
      out.qr <- c(out.coef[2,c(1,2)], out.df, out.coef[2,1]-qt(0.975,out.df)*out.coef[2,2], out.coef[2,1]+qt(0.975,out.df)*out.coef[2,2], out.coef[2,4])
      qr.out[j,] <- out.qr
    }
    qr.out <- as.data.frame(qr.out)
    rownames(qr.out) <- colnames(taxa)
    colnames(qr.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    taxa.qr.out[[i]] <- qr.out
  }
  names(taxa.qr.out) <- names(taxa)
  return(taxa.qr.out)
}

###########################
# Continuous - Covariates #
###########################

taxa.con.cov.recode.func <- function(sam.dat, sel.con.var, sel.cov.var, rename.con.var, taxa) {  #use rare.sam.dat for sam.dat
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  cov.var <- as.data.frame(as.matrix(sam.dat[ind.nona, sel.cov.var]))
  colnames(cov.var) <- sel.cov.var
  taxa.dat <- list()
  for (i in 1:6) {
    taxon <- taxa[[i]]
    taxon <- taxon[ind.nona,]
    taxa.dat[[i]] <- taxon
  }
  names(taxa.dat) <- names(taxa)
  return(list(con.var = con.var, cov.var = cov.var, taxa = taxa.dat))
}

taxa.con.cov.beta.recode.func <- function(sam.dat, sel.con.var, sel.cov.var, rename.con.var, taxa) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  cov.var <- as.data.frame(as.matrix(sam.dat[ind.nona, sel.cov.var]))
  colnames(cov.var) <- sel.cov.var
  taxa.dat <- list()
  
  for (i in 1:6) {
    taxon <- taxa[[i]]
    taxa.dat[[i]] <- taxon[ind.nona,]
    #taxa.dat[[i]] <- zCompositions::cmultRepl(taxon)
  }
  names(taxa.dat) <- names(taxa)
  
  return(list(con.var = con.var, cov.var = cov.var, taxa = taxa.dat))
}

taxa.lm.con.cov.func <- function(con.var, cov.var, taxa) {
  taxa.lm.out <- list()
  for (j in 1:6) {
    if(!is.null(taxa[[j]])){
      taxa.rank <- taxa[[j]]
      n.cov <- ncol(cov.var)
      n.tax <- ncol(taxa.rank)
      lm.out <- matrix(NA, n.tax, 6)
      d <- cbind(con.var, cov.var, taxa.rank)
      for (i in 1:n.tax) {
        taxon <- taxa.rank[,i]
        lm.f <- formula(paste("taxon", "~", colnames(con.var), "+", paste(colnames(cov.var), collapse = "+")))
        fit <- try(lm(lm.f, data = d), silent = TRUE)
        
        est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
        std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
        df <- tryCatch(summary(fit)$df[2], error = function(err) NA)
        if(is.na(df)){
          ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) C(NA,NA))
        }else{
          ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) C(NA,NA))
        }
        p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
        out <- c(est, std.err, df, ci, p.val)
        
        lm.out[i,] <- out
      }
      lm.out <- as.data.frame(lm.out)
      rownames(lm.out) <- colnames(taxa.rank)
      colnames(lm.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
      taxa.lm.out[[j]] <- lm.out
    }
  }
  names(taxa.lm.out) <- names(taxa)
  return(taxa.lm.out)
}

taxa.con.cov.glm.nb.func <- function(con.var, cov.var, taxa, library.size) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    # library.size <- apply(taxa,1, sum)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(con.var, cov.var, taxon))
      f <- formula(paste(colnames(dat)[ncol(dat)], " ~ ", colnames(dat)[1], "+", paste(names(cov.var), collapse = "+"), "+ offset(log(library.size))", sep = ""))
      fit <- try(glm.nb(f, data = dat), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
      df <- NA
      ci <- c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err)
      p.val <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
      
      out <- c(est, std.err, df, ci, p.val)
      lmer.out[i,] <- out
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.con.cov.glm.nb <- function(con.var, cov.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    count.taxa.test.out <- taxa.con.cov.glm.nb.func(con.var = con.var, cov.var = cov.var, taxa = taxa[[i]], library.size)
    count.taxa.test.q.out <- bin.q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <- count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.con.cov.beta <- function(con.var, cov.var, taxa, multi.method = "BH"){
  tax.out <- list()
  for (i in 1:6) {
    prop.taxa.test.out <- taxa.con.cov.beta.func(con.var, cov.var, taxa[[i]])
    prop.taxa.test.q.out <- bin.q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.con.cov.beta.func <- function(con.var, cov.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    beta.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(con.var, cov.var, taxon))
      dat[,2] <- as.numeric(dat[,2])
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[ncol(dat)], " ~ ", colnames(dat)[1], "+", paste(colnames(cov.var), collapse = "+"), sep = ""))
      fit.beta <- try(betareg(f, data = dat, maxit = 50), silent = TRUE)
      if(class(fit.beta) != "try-error"){
        est <- summary(fit.beta)$coefficients$mean[2,1]
        std.err <- summary(fit.beta)$coefficients$mean[2,2]
        ci <- c(est + qnorm(c(0.025, 0.975))*std.err)
        out.beta <- c(summary(fit.beta)$coefficients$mean[2,c(1,2)], NA, ci, summary(fit.beta)$coefficients$mean[2,4])
      } else {
        out.beta <- rep(NA, 6)
      }
      beta.out[i,] <- out.beta
    }
    beta.out <- as.data.frame(beta.out)
    rownames(beta.out) <- colnames(taxa)
    colnames(beta.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(beta.out)
  }
}
  
###################
# Other functions #
###################

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  return(x)
}

taxa.names.rank <- function(taxa.out){ 
  taxon.names <- list()
  taxon.names <- lapply(taxa.out, function(x) str_split(names(x), ";"))
  dup.list <- list(NA,NA,NA,NA,NA,NA)
  
  ranks <- c("K_", "P_", "C_", "O_", "F_", "G_", "S_")
  
  taxon.names.rank <- list()
  for(rank in 1:6){
    taxon <- lapply(taxon.names[[rank]], function(x) str_sub(x,start = 3))
    taxon.names.rank[[rank]] <- sapply(taxon, tail, 1)
    
    if(length(taxon.names.rank[[rank]]) != length(unique(taxon.names.rank[[rank]]))){
      duplicated.taxons <- unique(taxon.names.rank[[rank]][duplicated(taxon.names.rank[[rank]])])
      
      for(i in 1:length(duplicated.taxons)){
        duplicated.taxon <- duplicated.taxons[i]
        ind.dup <- which(taxon.names.rank[[rank]] %in% duplicated.taxon)
        
        if(duplicated.taxon == "NA") {
          for(k in 1:length(ind.dup)){
            #duplicated.taxon <- paste(duplicated.taxon, k,collapse = "")
            duplicated.taxon <- paste0("NA",k)
            taxon.names.rank[[rank]][ind.dup[k]] <- duplicated.taxon 
            dup.list[[rank]][k] <- paste(duplicated.taxon, " : ", paste(paste(ranks[1:(rank+1)], unlist(taxon[ind.dup[k]]), sep = ""), collapse = " | "), sep = "")
          }
        } else {
          for(j in 1:length(ind.dup)){
            duplicated.taxon <- paste(duplicated.taxon,"*",collapse = "")
            taxon.names.rank[[rank]][ind.dup[j]] <- duplicated.taxon 
            dup.list[[rank]][j] <- paste(duplicated.taxon, " : ", paste(paste(ranks[1:(rank+1)], unlist(taxon[ind.dup[j]]), sep = ""), collapse = " | "), sep = "")
          }
        }
        
      }
    }
  }
  names(taxon.names.rank) <- names(taxa.out)
  return(list(names = taxon.names.rank, duplicates = dup.list))
}

###################
# Which variable? #
###################

taxa.sig.dend <- function(out, tax.tab, layout.type = "twopi", species.include = "FALSE") {
  
  names(out) <-c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  ci.tab.all <- matrix(NA, 1, 1)
  for (i in 1:6) {
    ind.sig <- which(out[[i]]$Q.value < 0.05)
    if (length(ind.sig) >= 1) {
      sig.out <- out[[i]][ind.sig,]
      taxa <- rownames(sig.out)
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,1], 3), nsmall = 3), 
                        format(round(sig.out[,"P.value"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.value"], 3), nsmall = 3))
      ci.tab <- cbind(sig.out[,1])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    }
  }
  
  text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
  tax.tab.num <- as.matrix(as.data.frame(tax.tab))
  tax.tab.numNA <- as.matrix(as.data.frame(tax.tab))
  
  for (i in 1:6) {
    ind <- which(tax.tab.num[,i+1] == "NANANA")
    if (length(ind) >= 1) {
      tax.tab.num[ind,i+1] <- NA
    }
  }
  rank <- 1
  for (tax.rank in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    ind.sig <- which(text.tab.all[,2] == tax.rank)
    rank.list <- c("p_", "c_", "o_", "f_", "g_", "s_")
    if (length(ind.sig) >= 1) {
      sig.taxa <- text.tab.all[ind.sig,3]
      ind.non.sig <- which(tax.tab.num[,tax.rank] %in% sig.taxa)
      
      if(length(ind.non.sig) >0) {
        non.sig.taxa <- names(table(tax.tab.num[-ind.non.sig,tax.rank]))
        for (i in 1:length(non.sig.taxa)) {
          ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
          tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
        }
      }
      
      for (i in 1:length(sig.taxa)) {
        ind <- which(text.tab.all[,3] == sig.taxa[i])
        ind.rank <- which(text.tab.all[,2] == tax.rank)
        ind <- ind[ind %in% ind.rank]
        sig.taxa.num <- text.tab.all[ind ,1]
        ind.num <- which(tax.tab.num[,tax.rank] == sig.taxa[i])
        tax.tab.num[ind.num,tax.rank] <- sig.taxa.num
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        
      }
      rank <- rank+1
      
    } else {
      non.sig.taxa <- names(table(tax.tab.num[,tax.rank]))
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      for (i in 1:length(non.sig.taxa)) {
        ind <- which(tax.tab.num[,tax.rank] == non.sig.taxa[i])
        tax.tab.num[ind,tax.rank] <- paste(tax.rank, i, sep = "", collapse = "")
      }
      ind.na <- grepl("_NA ", tax.tab.numNA[,tax.rank])
      
      if(sum(ind.na) > 0){
        na.name.rank <- unlist(lapply(str_split(tax.tab.numNA[ind.na,tax.rank], ";"), function(x) {x <- x[length(x)]}))
        tax.tab.numNA[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
        tax.tab.num[ind.na,tax.rank] <- str_remove_all(na.name.rank, rank.list[rank])
      }
      rank <- rank+1
    }
  }
  all.phylum <- tax.tab.num[,"Phylum"]
  all.class <- tax.tab.num[,"Class"]
  all.order <- tax.tab.num[,"Order"]
  all.family <- tax.tab.num[,"Family"]
  all.genus <- tax.tab.num[,"Genus"]
  all.species <- tax.tab.num[,"Species"]
  
  missing <- list()
  connection.with.upper <- list()
  connection <- character()
  
  tax.tab.num[,"Phylum"] <- as.character(lapply(tax.tab.num[,"Phylum"], function(x){ str_replace_all(x, "[[:punct:]]| |=", "")}))
  tax.tab.num[,"Class"] <- as.character(lapply(tax.tab.num[,"Class"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Order"] <- as.character(lapply(tax.tab.num[,"Order"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Family"] <- as.character(lapply(tax.tab.num[,"Family"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Genus"] <- as.character(lapply(tax.tab.num[,"Genus"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  tax.tab.num[,"Species"] <- as.character(lapply(tax.tab.num[,"Species"], function(x){ str_replace_all(x,"[[:punct:]]| |=", "")}))
  
  Phylum <- names(table(tax.tab.num[,"Phylum"]))
  Phylum.desc <- c()
  Phylum.desc.dum <- c()
  
  sig.Phylum <- Phylum[!grepl("Phylum", Phylum)]
  ind <- as.numeric(sig.Phylum)+1
  
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)
  
  pos.sig.Phylum <- sig.Phylum[ind.pos.sig] 
  neg.sig.Phylum <- sig.Phylum[ind.neg.sig]
  na.Phylum <- Phylum[grepl("NA", Phylum)]
  non.sig.Phylum <- Phylum[!is.element(Phylum,c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum))]
  
  dummy <- c()
  for(i in 1:sum(!is.na(Phylum))){
    dummy <- c(dummy, paste0("dum",i))
  }
  
  phylum.wdum1 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  phylum.wdum2 <- c(pos.sig.Phylum, neg.sig.Phylum, na.Phylum, non.sig.Phylum)
  
  for(i in seq(1,length(Phylum),2)){
    phylum.wdum1[i] <- dummy[i]
  }
  for(i in seq(2,length(Phylum),2)){
    phylum.wdum2[i] <- dummy[i]
  }
  
  phylum.wdum <- rbind(phylum.wdum1, phylum.wdum2)
  k <- 1
  for (i in 1:length(Phylum)) {                                 # just include the following for loops inside this for loop and save each graph as pdf
    ind <- which(tax.tab.num[,"Phylum"] == Phylum[i])           # for (i in 1: length(phylum)) {~~~~, for (i in 1:length(Class))~~~~, for (i in 1: length(Order)~~....)}
    if (sum(!is.na(tax.tab.num[ind,"Class"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Class"]))
      desc <- str_replace_all(desc, " ", "")
      Phylum.desc[i] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")
      
      if(Phylum[i] %in% phylum.wdum2){
        Phylum.desc.dum[k] <- paste(Phylum[i], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }else{
        Phylum.desc.dum[k] <- paste(phylum.wdum2[which(phylum.wdum1==Phylum[i])], " -> {", paste(desc, collapse = " "), "} ", sep = "", collapse = "")  
        k <- k+1
      }
      connection <- c(connection, desc)
    }
    connection.with.upper[[1]] <- connection
  }
  
  Phylum.desc <- Phylum.desc[!is.na(Phylum.desc)]
  
  Bacteria.desc.dum <- paste("Bacteria -> {", paste(phylum.wdum1, collapse = " "), "}", sep = "", collapse = "")
  phy.to.phy <- paste(phylum.wdum1, " -> {", phylum.wdum2, "} ", sep = "")
  
  Phylum.desc <- paste(Phylum.desc, collapse = "; ")
  
  connection <- character()
  Class <- names(table(tax.tab.num[,"Class"]))
  Class.desc <- c()
  for (i in 1:length(Class)) {
    ind <- which(tax.tab.num[,"Class"] == Class[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Order"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Order"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Class.desc[i] <- paste(Class[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[2]] <- connection
  }
  
  Class.desc <- Class.desc[!is.na(Class.desc)]
  Class.desc <- paste(Class.desc, collapse = "; ")
  
  connection <- character()
  Order <- names(table(tax.tab.num[,"Order"]))
  Order <- Order[!Order %in% "NA"]
  Order.desc <- c()
  for (i in 1:length(Order)) {
    ind <- which(tax.tab.num[,"Order"] == Order[i])
    
    if (sum(!is.na(tax.tab.num[ind,"Family"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Family"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Order.desc[i] <- paste(Order[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[3]] <- connection
  }
  Order.desc <- Order.desc[!is.na(Order.desc)]
  Order.desc <- paste(Order.desc, collapse = "; ")
  
  connection <- character()
  Family <- names(table(tax.tab.num[,"Family"]))
  Family <- Family[!Family %in% "NA"]
  Family.desc <- c()
  for (i in 1:length(Family)) {
    ind <- which(tax.tab.num[,"Family"] == Family[i])
    if (sum(!is.na(tax.tab.num[ind,"Genus"])) >= 1) {
      desc <- names(table(tax.tab.num[ind,"Genus"]))
      desc <- str_replace_all(desc, " ", "")
      desc <- desc[!(desc %in% "NA")]
      Family.desc[i] <- paste(Family[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
      connection <- c(connection, desc)
    }
    connection.with.upper[[4]] <- connection
  }
  Family.desc <- Family.desc[!Family.desc %in% "NA"]
  Family.desc <- paste(Family.desc, collapse = "; ")
  
  connection <- character()
  Genus <- names(table(tax.tab.num[,"Genus"]))
  Genus <- Genus[!Genus %in% "NA"]
  Genus.desc <- c()
  
  if(length(setdiff(Genus,connection.with.upper[[4]])) > 0){
    for (i in 1:length(Genus)) {
      if(Genus[i] != setdiff(Genus, connection.with.upper[[4]])){
        ind <- which(tax.tab.num[,"Genus"] == Genus[i])
        
        if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
          desc <- names(table(tax.tab.num[ind,"Species"]))
          desc <- str_replace_all(desc, " ", "")
          desc <- desc[!(desc %in% "NA")]
          Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
          connection <- c(connection, desc)
        }
        connection.with.upper[[5]] <- connection 
      }
    }
  } else if(length(setdiff(Genus,connection.with.upper[[4]])) == 0){
    for (i in 1:length(Genus)) {
      ind <- which(tax.tab.num[,"Genus"] == Genus[i])
      
      if (sum(!is.na(tax.tab.num[ind,"Species"])) >= 1) {
        desc <- names(table(tax.tab.num[ind,"Species"]))
        desc <- str_replace_all(desc, " ", "")
        desc <- desc[!(desc %in% "NA")]
        Genus.desc[i] <- paste(Genus[i], " -> {", paste(desc, collapse = " "), "}", sep = "", collapse = "")    
        connection <- c(connection, desc)
      }
      connection.with.upper[[5]] <- connection 
    }
  }
  Genus.desc <- Genus.desc[!is.na(Genus.desc)]
  Genus.desc <- paste(Genus.desc, collapse = "; ")
  
  #Species <- names(table(tax.tab.num[,"Species"]))
  Class <- connection.with.upper[[1]]
  Order <- connection.with.upper[[2]]
  Family <- connection.with.upper[[3]]
  Genus <- connection.with.upper[[4]]
  Species <- connection.with.upper[[5]]
  
  missing[[1]] <- NA
  missing[[2]] <- setdiff(all.class, connection.with.upper[[1]])
  missing[[3]] <- setdiff(all.order, connection.with.upper[[2]])
  missing[[4]] <- setdiff(all.family, connection.with.upper[[3]])
  missing[[5]] <- setdiff(all.genus, connection.with.upper[[4]])
  missing[[6]] <- setdiff(all.species, connection.with.upper[[5]])
  
  Class <- Class[!(Class %in% "NA")]
  sig.Class <- Class[!grepl("Class", Class)]
  #ind <- text.tab.all[,1] %in% Class[!grepl("Class", Class)]
  ind <- as.numeric(sig.Class)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Class <- sig.Class[ind.pos.sig] 
  neg.sig.Class <- sig.Class[ind.neg.sig]
  
  Order <- Order[!(Order %in% "NA")]
  sig.Order <- Order[!grepl("Order", Order)]
  #ind <- text.tab.all[,1] %in% Order[!grepl("Order", Order)]
  ind <- as.numeric(sig.Order)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Order <- sig.Order[ind.pos.sig] 
  neg.sig.Order <- sig.Order[ind.neg.sig]
  
  Family <- Family[!(Family %in% "NA")]
  sig.Family <- Family[!grepl("Family", Family)]
  na.Family <- Family[grepl("NA", Family)]
  #ind <- text.tab.all[,1] %in% Family[!grepl("Family", Family)]
  ind <- as.numeric(sig.Family)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Family <- as.character(sig.Family[ind.pos.sig])
  neg.sig.Family <- as.character(sig.Family[ind.neg.sig])
  
  Genus <- Genus[!(Genus %in% "NA")]
  sig.Genus <- Genus[!grepl("Genus", Genus)]
  ind <- as.numeric(sig.Genus)+1
  ind.pos.sig <- which(ci.tab.all[ind,1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind,1] < 0)  
  pos.sig.Genus <- sig.Genus[ind.pos.sig] 
  neg.sig.Genus <- sig.Genus[ind.neg.sig]
  
  Species <- Species[!(Species %in% "NA")]
  sig.Species <- Species[!grepl("Species", Species)]
  ind <- as.numeric(sig.Species)+1
  ind.pos.sig <- which(ci.tab.all[ind[!is.na(ind)],1] >= 0)
  ind.neg.sig <- which(ci.tab.all[ind[!is.na(ind)],1] < 0)  
  
  pos.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  neg.sig.Species <- sig.Species[!grepl("NA", sig.Species)]
  
  pos.sig.Species <- pos.sig.Species[ind.pos.sig] 
  neg.sig.Species <- neg.sig.Species[ind.neg.sig]
  
  if (species.include) {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, Genus.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus, pos.sig.Species)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus, neg.sig.Species)
    total.group <- c(Class, Order, Family, Genus, Species)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
    
  } else {
    Bacteria.desc.dum <- paste(Bacteria.desc.dum, sep = "; ", collapse = "")
    phy.to.phy <- paste(phy.to.phy, sep = "; ", collapse = "")
    Phylum.desc.dum <- paste(Phylum.desc.dum, sep = "; ", collapse = "")
    
    Taxa.desc <- paste(Bacteria.desc.dum, phy.to.phy, Phylum.desc.dum, Class.desc, Order.desc, Family.desc, sep = "; ", collapse = "")
    pos.sig.taxa <- c(pos.sig.Class, pos.sig.Order, pos.sig.Family, pos.sig.Genus)
    neg.sig.taxa <- c(neg.sig.Class, neg.sig.Order, neg.sig.Family, neg.sig.Genus)
    total.group <- c(Class, Order, Family, Genus)
    total.group <- total.group[!total.group %in% "NA"]
    na.taxa <- total.group[grepl("NA", total.group)]
    non.sig.taxa <- total.group[!total.group %in% c(pos.sig.taxa, neg.sig.taxa, na.taxa)] 
  }
  Taxa.desc <- str_remove_all(Taxa.desc, " NA;")
  
  flow.text.dum <- list()
  length(flow.text.dum) <- 6
  
  neg <- 1
  non <- 1
  na <- 1
  for(i in 1:length(phylum.wdum1)){
    if(phylum.wdum1[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum2[i] %in% pos.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(pos.sig.Phylum[i], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = indianred, label = ",pos.sig.Phylum[i],"]", sep = " "))
    }else if(phylum.wdum1[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum2[i] %in% neg.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(neg.sig.Phylum[neg], paste("[fontname = Helvetica, shape = circle, fontsize = 11, fixedsize = true, width = 0.3, style = filled, fillcolor = steelblue, label = ",neg.sig.Phylum[neg],"]", sep = " "))
      neg <- neg+1
    }else if(phylum.wdum1[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }else if(phylum.wdum2[i] %in% non.sig.Phylum){
      flow.text.dum[[1]][i] <- paste(non.sig.Phylum[non], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
      non <- non+1
    }
    else if(phylum.wdum1[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }else if(phylum.wdum2[i] %in% na.Phylum){
      flow.text.dum[[1]][i] <- paste(na.Phylum[na], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
      na <- na+1
    }
  }
  
  if(length(pos.sig.taxa) >0){
    for(i in 1:length(pos.sig.taxa)){
      flow.text.dum[[2]][i] <- paste(pos.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = indianred, label = ",pos.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(neg.sig.taxa) >0){
    for(i in 1:length(neg.sig.taxa)){
      flow.text.dum[[3]][i] <- paste(neg.sig.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 10, fixedsize = true, width = 0.25, style = filled, fillcolor = steelblue, label = ",neg.sig.taxa[i],"]", sep = " "))
    }
  }
  
  if(length(non.sig.taxa) >0){
    for(i in 1:length(non.sig.taxa)){
      flow.text.dum[[4]][i] <- paste(non.sig.taxa[i], paste("[fontname = Helvetica, shape = point, fixedsize = true, width = 0.1, style = filled, fillcolor = grey, label = '']", sep = " "))
    }
  }
  
  if(length(na.taxa) >0){
    for(i in 1:length(na.taxa)){
      flow.text.dum[[5]][i] <- paste(na.taxa[i], paste("[fontname = Helvetica, shape = circle, fontsize = 5, fixedsize = true, width = 0.1, style = '', label = '']", sep = " "))
    }
  }
  
  if(length(dummy) >0){
    for(i in 1:length(dummy)){
      flow.text.dum[[6]][i] <- paste(dummy[i], paste("[fontname = Helvetica, shape = point, width = 0, style = invi, label = '']", sep = " "))
    }
  }
  
  flow.text.dum.complete <- c()
  for(i in 1:6){
    flow.text.dum.complete <- paste(flow.text.dum.complete, flow.text.dum[[i]], sep = "\n ")
  }
  
  flow.text <- paste("digraph ", layout.type, " {",
                     paste("graph [layout = ",layout.type,", root = Bacteria, ranksep = \"2 0.3 1 1 1 \"]", collapse = ""),
                     "Bacteria [fontname = Helvetica, fontsize = 12, fixedsize = true, width = 1, style = filled, fillcolor = lightblue, label = Bacteria]",
                     paste(flow.text.dum.complete, collapse = ""),
                     "edge [color = black, arrowhead = none]",
                     Taxa.desc,
                     "}", sep = "\n")
  
  taxon.tab <- data.frame("ID" = text.tab.all[-1,1],
                          #"Status" = ifelse( ci.tab.all[-1] < 0, "neg", "pos"), # "blue", "red"),
                          "Taxon" = sub(".*;", "", text.tab.all[-1,3]))
  
  return(list(flow.text = grViz(flow.text), taxon.tab = taxon.tab, ci.tab.all = ci.tab.all) )
}