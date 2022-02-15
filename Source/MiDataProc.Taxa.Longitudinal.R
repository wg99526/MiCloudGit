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
library(lme4)
library(lmerTest)
library(MiRKAT)
library(DiagrammeR)
library(stringr)
library(devtools)
library(reticulate)
#library(broom.mixed)
library(NBZIMM)
library(nlme)
#library(pscl)
#library(betareg)
library(gee)
library(geepack)
library(gridExtra)
library(glmmTMB)
library(glmm)

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
        
        for(j in 1:length(ind.dup)){
          duplicated.taxon <- paste(duplicated.taxon,"*",collapse = "")
          taxon.names.rank[[rank]][ind.dup[j]] <- duplicated.taxon 
          dup.list[[rank]][j] <- paste(duplicated.taxon, " : ", paste(paste(ranks[1:(rank+1)], unlist(taxon[ind.dup[j]]), sep = ""), collapse = " | "), sep = "")
        }
      }
    }
  }
  names(taxon.names.rank) <- names(taxa.out)
  return(list(names = taxon.names.rank, duplicates = dup.list))
}

#####################
# Data manipulation #
#####################

cluster.func <- function(sam.dat, mon.sin.rev.bin.con, selected.var) {
  ind.pri <- colnames(sam.dat) %in% selected.var
  ind.mon.sin.rev <- mon.sin.rev.bin.con$is.mon | mon.sin.rev.bin.con$is.rev | mon.sin.rev.bin.con$is.sin
  return(colnames(sam.dat)[!(ind.pri | ind.mon.sin.rev)])
}

taxa.bin.var.func <- function(sam.dat) {
  var.names <- colnames(sam.dat)
  return(var.names)
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

taxa.bin.id.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sel.id.var, sam.dat, taxa) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  id.var <- unlist(sam.dat[,sel.id.var])
  id.var <- id.var[c(ind.ref, ind.com)]
  #bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxon <- rbind(taxa[ind.ref,], taxa[ind.com,])
  taxa.out <- taxon
  
  names(taxa.out) <- names(taxa)
  
  return(list(bin.var = bin.var, id.var = id.var, taxa = taxa.out))
}

taxa.bin.id.cat.ref.united.func <- function(sel.bin.var, sel.id.var, sel.ref, sel.com, sam.dat, taxa) {
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  id.var <- unlist(sam.dat[,sel.id.var])
  id.var <- id.var[c(ind.ref, ind.com)]
  taxa.out <- list()
  for(i in 1:6){
    taxa.out[[i]] <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
  }
  names(taxa.out) <- names(taxa)
  
  return(list(bin.var = bin.var, id.var = id.var, taxa = taxa.out))
}

# taxa.bin.id.cov.cat.ref.func <- function(sel.bin.var, sel.id.var, sel.cov.var, sel.ref, sel.com, sam.dat, taxa) {
#   bin.var <- unlist(sam.dat[,sel.bin.var])
#   ind.ref <- which(bin.var == sel.ref)
#   ind.com <- which(bin.var == sel.com)
#   bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
#   id.var <- unlist(sam.dat[,sel.id.var])
#   id.var <- id.var[c(ind.ref, ind.com)]
#   cov.var <- sam.dat[,sel.cov.var]
#   cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
#   taxa <- rbind(taxa[ind.ref,], taxa[ind.com,])
#   
#   return(list(bin.var = bin.var, id.var = id.var, cov.var = cov.var, taxa = taxa))
# }

taxa.con.id.recode.func <- function(sam.dat, sel.con.var, sel.id.var, rename.con.var, taxa) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  id.var <- unlist(sam.dat[ind.nona,sel.id.var])
  taxa <- taxa[ind.nona,]
  return(list(con.var = con.var, id.var = id.var, taxa = taxa))
}

taxa.con.id.recode.united.func <- function(sam.dat, sel.con.var, sel.id.var, rename.con.var, taxa) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  id.var <- unlist(sam.dat[ind.nona,sel.id.var])
  taxa.out <- list()
  for(i in 1:6){
    taxa.out[[i]] <- taxa[[i]][ind.nona,]
  }
  names(taxa.out) <- names(taxa)
  return(list(con.var = con.var, id.var = id.var, taxa = taxa.out))
}

taxa.con.id.cov.recode.func <- function(sam.dat, sel.con.var, rename.con.var, sel.id.var, sel.cov.var, taxa) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  id.var <- unlist(sam.dat[ind.nona,sel.id.var])
  cov.var <- sam.dat[ind.nona,sel.cov.var]
  taxa <- taxa[ind.nona,]
  return(list(con.var = con.var, id.var = id.var, cov.var = cov.var, taxa = taxa))
}

taxa.con.id.cov.recode.united.func <- function(sam.dat, sel.con.var, rename.con.var, sel.id.var, sel.cov.var, taxa) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  id.var <- unlist(sam.dat[ind.nona,sel.id.var])
  cov.var <- sam.dat[ind.nona,sel.cov.var]
  taxa.out <- list()
  for(i in 1:6){
    taxa.out[[i]] <- taxa[[i]][ind.nona,]
  }
  names(taxa.out) <- names(taxa)
  return(list(con.var = con.var, id.var = id.var, cov.var = cov.var, taxa = taxa.out))
}

######################
# Summary statistics #
######################

taxa.ind.sum.func <- function(x) {
  sum.out <- c(length(x), mean(x), quantile(x))
  names(sum.out) <-  c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  return(sum.out)
}

taxa.bin.id.sum.func <- function(bin.var, taxa, sel.ref, sel.com) {
  if(!is.null(taxa)){
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
    names(out) <- c(sel.ref, sel.com)
    
    return(out)
  }
}

taxa.bin.id.sum.united.func <- function(bin.var, taxa.out, sel.ref, sel.com) {
  taxa.bin.sum <-list()
  for(i in 1:6) {
    taxa.bin.sum[[i]] <- taxa.bin.id.sum.func(bin.var, taxa.out[[i]], sel.ref, sel.com)
  }
  names(taxa.bin.sum) <- names(taxa.out)
  return(taxa.bin.sum)
}

#################
# Data analysis #
#################

# Random intercept model

taxa.bin.lmer.func <- function(bin.var, id.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(taxon, bin.var, id.var))
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[1], " ~ ", colnames(dat)[2], "+ (1 | ", colnames(dat)[3], ")", sep = ""))
      
      fit <- try(lmer(f, data = dat), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
      #ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      df <- tryCatch(summary(fit)$coefficients[2,3], error = function(err) NA)
      if(is.na(df)){
        ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      }else{
        ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      }
      pvs <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
      out <- c(est, std.err, df, ci, pvs)
      lmer.out[i,] <- out
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

taxa.bin.cov.lmer.func <- function(bin.var, id.var, cov.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(taxon, bin.var, id.var, cov.var))
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[1], " ~ ", colnames(dat)[2], "+ (1 | ", colnames(dat)[3], ") + ", paste(colnames(dat)[4:ncol(dat)], collapse = " + "), sep = ""))
      
      fit <- try(lmer(f, data = dat), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
      df <- tryCatch(summary(fit)$coefficients[2,3], error = function(err) NA)
      if(is.na(df)){
        ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      }else{
        ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      }
      pvs <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
      out <- c(est, std.err, df, ci, pvs)
      
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

taxa.con.lmer.func <- function(con.var, id.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(taxon, con.var, id.var))
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[1], " ~ ", colnames(dat)[2], "+ (1 | ", colnames(dat)[3], ")", sep = ""))
      
      fit <- try(lmer(f, data = dat), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
      df <- tryCatch(summary(fit)$coefficients[2,3], error = function(err) NA)
      if(is.na(df)){
        ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      }else{
        ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      }
      pvs <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
      out <- c(est, std.err, df, ci, pvs)
      lmer.out[i,] <- out
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

taxa.con.cov.lmer.func <- function(con.var, id.var, cov.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(taxon, con.var, id.var, cov.var))
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[1], " ~ ", colnames(dat)[2], "+ (1 | ", colnames(dat)[3], ") + ", paste(colnames(dat)[4:ncol(dat)], collapse = " + "), sep = ""))
      fit <- try(lmer(f, data = dat), silent = TRUE)
      
      est <- tryCatch(summary(fit)$coefficients[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,2], error = function(err) NA)
      df <- tryCatch(summary(fit)$coefficients[2,3], error = function(err) NA)
      if(is.na(df)){
        ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      }else{
        ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      }
      
      pvs <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
      out <- c(est, std.err, df, ci, pvs)
      lmer.out[i,] <- out
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

#####################################
# Random effects model for all taxa #
#####################################

all.taxa.bin.lmer <- function(bin.var, id.var, taxa, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    clr.taxa.test.out <- taxa.bin.lmer.func(bin.var = bin.var, id.var = id.var, taxa = taxa[[i]])
    clr.taxa.test.q.out <- q.func(clr.taxa.test.out, method = multi.method)
    tax.out[[i]] <- clr.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.lmer.sig.graph <- function(out) {
  
  text.tab.all <- matrix(c("Rank", "Taxon", "Est.", "P-value", "Q-value"), 1, 5)
  ci.tab.all <- matrix(c(NA, NA, NA), 1, 3)
  for (i in 1:6) {
    ind.sig <- which(out[[i]]$Q.value < 0.05)
    if (length(ind.sig) >= 1) {
      sig.out <- out[[i]][ind.sig,]
      taxa <- rownames(sig.out)
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,"Est"], 3), nsmall = 3), 
                        format(round(sig.out[,"P.value"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.value"], 3), nsmall = 3))
      ci.tab <- cbind(sig.out[,"Est"], sig.out[,"Lower"], sig.out[,"Upper"])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    } 
  }
  
  if(nrow(ci.tab.all) >=2) {
    text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.1, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
    plot.taxa <- grid.grab()
  } else {
    warning("No significant cases found")
  }
}

all.taxa.lmer.logit.sig.graph <- function(out, mult.test.cor = TRUE) {
  
  text.tab.all <- matrix(c("Rank", "Taxon", "OR", "P-value", "Q-value"), 1, 5)
  ci.tab.all <- matrix(c(NA, NA, NA), 1, 3)
  for (i in 1:6) {
    if(mult.test.cor) {
      ind.sig <- which(out[[i]]$Q.value < 0.05)
    } else {
      ind.sig <- which(out[[i]]$P.value < 0.05)
    }
    if (length(ind.sig) >= 1) {
      sig.out <- out[[i]][ind.sig,]
      taxa <- rownames(sig.out)
      text.tab <- cbind(rep(names(out)[i], length(ind.sig)), taxa, format(round(sig.out[,"OR"], 3), nsmall = 3), 
                        format(round(sig.out[,"P.value"], 3), nsmall = 3), 
                        format(round(sig.out[,"Q.value"], 3), nsmall = 3))
      ci.tab <- cbind(sig.out[,"OR"], sig.out[,"Lower"], sig.out[,"Upper"])
      text.tab.all <- rbind(text.tab.all, text.tab)
      ci.tab.all <- rbind(ci.tab.all, ci.tab)
    } 
  }
  if(nrow(ci.tab.all) >=2) {
    text.tab.all <- cbind(c("ID", 1:(nrow(text.tab.all)-1)), text.tab.all)
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.1, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
    plot.taxa <- grid.grab()
  } else {
    warning("No significant cases found")
  }
}

all.taxa.bin.cov.lmer <- function(bin.var, id.var, cov.var, taxa, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    clr.taxa.test.out <- taxa.bin.cov.lmer.func(bin.var = bin.var, id.var = id.var, cov.var = cov.var, taxa = taxa[[i]])
    clr.taxa.test.q.out <- q.func(clr.taxa.test.out, method = multi.method)
    tax.out[[i]] <- clr.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.con.lmer <- function(con.var, id.var, taxa, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    clr.taxa.test.out <- taxa.con.lmer.func(con.var = con.var, id.var = id.var, taxa = taxa[[i]])
    clr.taxa.test.q.out <- q.func(clr.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  clr.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

all.taxa.con.cov.lmer <- function(con.var, id.var, cov.var, taxa, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    clr.taxa.test.out <- taxa.con.cov.lmer.func(con.var = con.var, id.var = id.var, cov.var = cov.var, taxa = taxa[[i]])
    clr.taxa.test.q.out <- q.func(clr.taxa.test.out, method = multi.method)
    tax.out[[i]] <- clr.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

#######################
# Logistic Regression #
#######################

get.or.se.lmer <- function(model) {
  t <- broom::tidy(model)
  or = exp(t$estimate[c(1,2)])
  var.diag = diag(vcov(model))
  or.se = sqrt(or^2 * var.diag)
  
  return(or.se = unname(unlist(or.se)))
}

get.or.se <- function(model) {
  broom::tidy(model) %>%
    mutate(or = exp(estimate),
           var.diag = diag(vcov(model)),
           or.se = sqrt(or^2 * var.diag)) %>%
    select(or.se) %>% unlist %>% unname
}

taxa.bin.id.cov.cat.ref.func <- function(sel.bin.var, sel.id.var, sel.cov.var, sel.ref, sel.com, sam.dat, taxa) {
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  id.var <- unlist(sam.dat[,sel.id.var])
  id.var <- id.var[c(ind.ref, ind.com)]
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  taxa <- rbind(taxa[ind.ref,], taxa[ind.com,])
  
  return(list(bin.var = bin.var, id.var = id.var, cov.var = cov.var, taxa = taxa))
}

taxa.bin.id.cov.cat.ref.united.func <- function(sel.bin.var, sel.id.var, sel.cov.var, sel.ref, sel.com, sam.dat, taxa) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  id.var <- unlist(sam.dat[,sel.id.var])
  id.var <- id.var[c(ind.ref, ind.com)]
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  #bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa.out <- list()
  for(i in 1:6) {
    taxon <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
    taxa.out[[i]] <- taxon
  }
  names(taxa.out) <- names(taxa)
  
  return(list(bin.var = bin.var, id.var = id.var, cov.var = cov.var, taxa = taxa.out))
}

taxa.bin.id.cat.ref.logit.united.func <- function(sel.bin.var, sel.ref, sel.com, sel.id.var, sam.dat, taxa) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  id.var <- unlist(sam.dat[,sel.id.var])
  id.var <- id.var[c(ind.ref, ind.com)]
  #bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  taxa.out <- list()
  for(i in 1:6) {
    taxon <- rbind(taxa[[i]][ind.ref,], taxa[[i]][ind.com,])
    taxa.out[[i]] <- taxon
  }
  names(taxa.out) <- names(taxa)
  
  return(list(bin.var = bin.var, id.var = id.var, taxa = taxa.out))
}

all.taxa.bin.logit.glmm.b <- function(bin.var, id.var, taxa, rare.count = FALSE, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    taxa.test.out <- taxa.bin.logit.glmm.b.func(bin.var = bin.var, id.var = id.var, taxa = taxa[[i]], rare.count = rare.count)
    taxa.test.q.out <- q.func(taxa.test.out, method = multi.method)
    tax.out[[i]] <- taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.logit.glmm.b.func <- function(bin.var, id.var, taxa, rare.count = FALSE) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, taxon, id.var))
      if(rare.count){
        f <- formula(paste(colnames(dat)[1], " ~ ", "scale(", colnames(dat)[2], ") + (1|", colnames(dat)[3], ")", sep = ""))
      }else{
        f <- formula(paste(colnames(dat)[1], " ~ ", colnames(dat)[2], " + (1|", colnames(dat)[3], ")", sep = ""))
      }
      
      fit <- try(glmer(f, data = dat, family = "binomial"(link = "logit")), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,"Std. Error"], error = function(err) NA)
      df <- tryCatch(summary(fit)$coefficients[2,3], error = function(err) NA)
      if(is.na(df)){
        ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      }else{
        ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      }
      pvs <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
      
      out.logit <-c(est, std.err, df, ci, pvs)
      lmer.out[i,] <- out.logit
      
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.bin.logit.gee.reg.coef <- function(bin.var, id.var, taxa, rare.count = FALSE, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    taxa.test.out <- taxa.bin.logit.gee.reg.coef.func(bin.var = bin.var, id.var = id.var, taxa = taxa[[i]], rare.count = rare.count)
    taxa.test.q.out <- q.func(taxa.test.out, method = multi.method)
    tax.out[[i]] <- taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.logit.gee.reg.coef.func <- function(bin.var, id.var, taxa, rare.count = FALSE) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, taxon, id.var))
      if(rare.count){
        f <- formula(paste(colnames(dat)[1], " ~ ", "scale(", colnames(dat)[2], ")", sep = ""))
      }else{
        f <- formula(paste(colnames(dat)[1], " ~ ", colnames(dat)[2], sep = ""))
      }
      
      fit <- try(gee(f, id = as.factor(id.var), data = dat, corstr="exchangeable", family = "binomial"(link = "logit")), silent = TRUE)
      
      est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,"Robust S.E."], error = function(err) NA)
      ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      z <- tryCatch(summary(fit)$coefficients[2,"Robust z"], error = function(err) NA)
      pvs <- tryCatch(2*pnorm(-abs(z), mean = 0, sd = 1), error = function(err) NA)
      
      out.logit <- c(est, std.err, z, ci, pvs)
      lmer.out[i,] <- out.logit
    }
    
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "z value", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.bin.cov.logit.glmm.b <- function(bin.var, id.var, cov.var, taxa, rare.count = FALSE, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    taxa.test.out <- taxa.bin.cov.logit.glmm.b.func(bin.var = bin.var, id.var = id.var, cov.var = cov.var, taxa = taxa[[i]], rare.count = rare.count)
    taxa.test.q.out <- q.func(taxa.test.out, method = multi.method)
    tax.out[[i]] <- taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.cov.logit.glmm.b.func <- function(bin.var, id.var, cov.var, taxa, rare.count = FALSE) {
  
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, taxon, id.var, cov.var))
      if(rare.count){
        f <- formula(paste(colnames(dat)[1], " ~ ", "scale(", colnames(dat)[2], ") + ", paste(colnames(cov.var), collapse = "+"), "+ (1|", colnames(dat)[3], ")", sep = ""))
      }else{
        f <- formula(paste(colnames(dat)[1], " ~ ", colnames(dat)[2], paste(colnames(cov.var), collapse = "+"), "+ (1|", colnames(dat)[3], ")", sep = ""))
      }
      
      fit <- try(glmer(f, data = dat, family = "binomial"(link = "logit")), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,"Std. Error"], error = function(err) NA)
      df <- tryCatch(summary(fit)$coefficients[2,3], error = function(err) NA)
      if(is.na(df)){
        ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      }else{
        ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      }
      pvs <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
      
      out.logit <-c(est, std.err, df, ci, pvs)
      lmer.out[i,] <- out.logit
      
    }
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.bin.cov.logit.gee.reg.coef <- function(bin.var, id.var, cov.var, taxa, rare.count = FALSE, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    taxa.test.out <- taxa.bin.cov.logit.gee.reg.coef.func(bin.var = bin.var, id.var = id.var, cov.var = cov.var, taxa = taxa[[i]], rare.count = rare.count)
    taxa.test.q.out <- q.func(taxa.test.out, method = multi.method)
    tax.out[[i]] <- taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.cov.logit.gee.reg.coef.func <- function(bin.var, id.var, cov.var, taxa, rare.count = FALSE) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, cov.var, taxon, id.var))
      if(rare.count){
        f <- formula(paste(colnames(dat)[1], " ~ ", "scale(taxon)", "+", paste(colnames(cov.var), collapse = "+"), sep = ""))
      }else{
        f <- formula(paste(colnames(dat)[1], " ~ ", "taxon", "+", paste(colnames(cov.var), collapse = "+"), sep = ""))
      }
      fit <- try(gee(f, id = as.factor(id.var), data = dat, corstr="exchangeable", family = "binomial"(link = "logit"), maxiter = 50), silent = TRUE)
      est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
      std.err <- tryCatch(summary(fit)$coefficients[2,"Robust S.E."], error = function(err) NA)
      ci <- tryCatch(c(est - qnorm(0.975)*std.err, est + qnorm(0.975)*std.err), error = function(err) c(NA,NA))
      z <- tryCatch(summary(fit)$coefficients[2,"Robust z"], error = function(err) NA)
      pvs <- tryCatch(2*pnorm(-abs(z), mean = 0, sd = 1), error = function(err) NA)
      
      out.logit <-c(est, std.err, z, ci, pvs)
      lmer.out[i,] <- out.logit
    }
    
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "z value", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

################################
# Negative Binomial Regression #
################################

all.taxa.bin.glmm.nb <- function(bin.var, id.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    count.taxa.test.out <- taxa.bin.glmm.nb.func(bin.var = bin.var, id.var = id.var, taxa = taxa[[i]], library.size)
    count.taxa.test.q.out <- q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <- count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.glmm.nb.func <- function(bin.var, id.var, taxa, library.size) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      y <- as.numeric(taxon)
      x <- as.numeric(unlist(bin.var))
      cluster <- id.var
      dat <- cbind.data.frame(y,x,cluster)
      
      fit.nb <- try(glmm.nb(y ~ x + offset(log(library.size)), random = ~1|cluster), silent = TRUE)
      est <- tryCatch(summary(fit.nb)$tTable[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit.nb)$tTable[2,2], error = function(err) NA)
      df <- tryCatch(summary(fit.nb)$tTable[2,3], error = function(err) NA)
      ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      p.val <- tryCatch(summary(fit.nb)$tTable[2,5], error = function(err) NA)
      
      out.nb <- c(est, std.err, df, ci, p.val)
      lmer.out[i,] <- out.nb
    }
    rownames(lmer.out) <- colnames(taxa)
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.bin.cov.glmm.nb <- function(bin.var, id.var, cov.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    count.taxa.test.out <- taxa.bin.cov.glmm.nb.func(bin.var = bin.var, id.var = id.var, cov.var = cov.var,taxa = taxa[[i]], library.size)
    count.taxa.test.q.out <- q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <- count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.cov.glmm.nb.func <- function(bin.var, id.var, cov.var, taxa, library.size) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      y <- as.numeric(taxon)
      x <- as.numeric(unlist(bin.var))
      cluster <- id.var
      dat <- cbind.data.frame(y,x,cov.var,cluster)
      
      f <- formula(paste("y ~ x +", paste(colnames(cov.var), collapse = "+"),"+ offset(log(library.size))", sep = ""))
      
      fit.nb <- try(glmm.nb(f, data = dat, random = ~1|cluster), silent = TRUE)
      est <- tryCatch(summary(fit.nb)$tTable[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit.nb)$tTable[2,2], error = function(err) NA)
      df <- tryCatch(summary(fit.nb)$tTable[2,3], error = function(err) NA)
      ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      p.val <- tryCatch(summary(fit.nb)$tTable[2,5], error = function(err) NA)
      
      out.nb <- c(est, std.err, df, ci, p.val)
      lmer.out[i,] <- out.nb
    }
    rownames(lmer.out) <- colnames(taxa)
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.con.glmm.nb <- function(con.var, id.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    count.taxa.test.out <- taxa.con.glmm.nb.func(con.var = con.var, id.var = id.var, taxa = taxa[[i]], library.size)  
    count.taxa.test.q.out <- q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <- count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.con.glmm.nb.func <- function(con.var, id.var, taxa, library.size) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      y <- as.numeric(taxon)
      x <- as.numeric(unlist(con.var))
      cluster <- id.var
      dat <- cbind.data.frame(y,x,cluster)
      
        fit.nb <- try(glmm.nb(y ~ x + offset(log(library.size)), data = dat, random = ~1|cluster), silent = TRUE)
        
        est <- tryCatch(summary(fit.nb)$tTable[2,1], error = function(err) NA)
        std.err <- tryCatch(summary(fit.nb)$tTable[2,2], error = function(err) NA)
        df <- tryCatch(summary(fit.nb)$tTable[2,3], error = function(err) NA)
        ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
        p.val <- tryCatch(summary(fit.nb)$tTable[2,5], error = function(err) NA)
        
        out.nb <- c(est, std.err, df, ci, p.val)
      
      lmer.out[i,] <- out.nb
    }
    rownames(lmer.out) <- colnames(taxa)
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

all.taxa.con.cov.glmm.nb <- function(con.var, id.var, cov.var, taxa, library.size, multi.method = "BH") {
  tax.out <- list()
  for (i in 1:6) {
    count.taxa.test.out <- taxa.con.cov.glmm.nb.func(con.var = con.var, id.var = id.var, cov.var = cov.var,taxa = taxa[[i]], library.size)
    count.taxa.test.q.out <- q.func(count.taxa.test.out, method = multi.method)
    tax.out[[i]] <- count.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.con.cov.glmm.nb.func <- function(con.var, id.var, cov.var, taxa, library.size) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    lmer.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      y <- as.numeric(taxon)
      x <- as.numeric(unlist(con.var))
      cluster <- id.var
      dat <- cbind.data.frame(y,x,cov.var,cluster)
      
      f <- formula(paste("y ~ x +", paste(colnames(cov.var), collapse = "+"), "+ offset(log(library.size))", sep = ""))
      
      fit.nb <- try(glmm.nb(f, data = dat, random = ~1|cluster), silent = TRUE)
      est <- tryCatch(summary(fit.nb)$tTable[2,1], error = function(err) NA)
      std.err <- tryCatch(summary(fit.nb)$tTable[2,2], error = function(err) NA)
      df <- tryCatch(summary(fit.nb)$tTable[2,3], error = function(err) NA)
      ci <- tryCatch(c(est - qt(0.975, df = df)*std.err, est + qt(0.975, df = df)*std.err), error = function(err) c(NA,NA))
      p.val <- tryCatch(summary(fit.nb)$tTable[2,5], error = function(err) NA)
      
      out.nb <- c(est, std.err, df, ci, p.val)
      lmer.out[i,] <- out.nb
    }
    rownames(lmer.out) <- colnames(taxa)
    lmer.out <- as.data.frame(lmer.out)
    rownames(lmer.out) <- colnames(taxa)
    colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
    return(lmer.out)
  }
}

###################
# Beta Regression #
###################

all.taxa.bin.glmm.beta <- function(bin.var, id.var, taxa, multi.method = "BH"){
  tax.out <- list()
  for (i in 1:6) {
    prop.taxa.test.out <- taxa.bin.glmm.beta.func(bin.var, id.var, taxa[[i]])
    prop.taxa.test.q.out <- q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.bin.glmm.beta.func <- function(bin.var, id.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    beta.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, taxon, id.var))
      dat[,2] <- as.numeric(dat[,2])
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[2], " ~ ", colnames(dat)[1], "+ (1 | id.var)", sep = ""))
      fit.beta <- try(glmmTMB(f, data = dat, family = beta_family()), silent = TRUE)       
      if(class(fit.beta) != "try-error"){
        est <- summary(fit.beta)$coefficients$cond[2,1]
        std.err <- summary(fit.beta)$coefficients$cond[2,2]
        ci <- c(est + qnorm(c(0.025, 0.975))*std.err)
        out.beta <- c(summary(fit.beta)$coefficients$cond[2,c(1,2)], NA, ci, summary(fit.beta)$coefficients$cond[2,4])
        
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

all.taxa.bin.glmm.cov.beta <- function(bin.var, id.var, cov.var, taxa, multi.method = "BH"){
  tax.out <- list()
  for (i in 1:6) {
    prop.taxa.test.out <- taxa.bin.glmm.cov.beta.func(bin.var, id.var, cov.var, taxa[[i]])
    prop.taxa.test.q.out <- q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
  
}

taxa.bin.glmm.cov.beta.func <- function(bin.var, id.var, cov.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    beta.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(bin.var, id.var, cov.var, taxon))
      f <- formula(paste(colnames(dat)[ncol(dat)], " ~ ", colnames(dat)[1], "+", paste(names(cov.var), collapse = "+"), "+ (1|id.var)", sep = ""))
      fit.beta <- try(glmmTMB(f, data = dat, family = beta_family()), silent = TRUE)
      if(class(fit.beta) != "try-error"){
        est <- summary(fit.beta)$coefficients$cond[2,1]
        std.err <- summary(fit.beta)$coefficients$cond[2,2]
        ci <- c(est + qnorm(c(0.025, 0.975))*std.err)
        out.beta <- c(summary(fit.beta)$coefficients$cond[2,c(1,2)], NA, ci, summary(fit.beta)$coefficients$cond[2,4])
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

all.taxa.con.glmm.beta <- function(con.var, id.var, taxa, multi.method = "BH"){
  tax.out <- list()
  for (i in 1:6) {
    prop.taxa.test.out <- taxa.con.glmm.beta.func(con.var, id.var, taxa[[i]])
    prop.taxa.test.q.out <- q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.con.glmm.beta.func <- function(con.var, id.var, taxa) {
  if (!is.null(taxa)){
    n.tax <- ncol(taxa)
    beta.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(con.var, taxon, id.var))
      dat[,2] <- as.numeric(dat[,2])
      dat[,1] <- as.numeric(dat[,1])
      f <- formula(paste(colnames(dat)[2], " ~ ", colnames(dat)[1], "+ (1 | id.var)", sep = ""))
      fit.beta <- try(glmmTMB(f, data = dat, family = beta_family()), silent = TRUE)  
      if (class(fit.beta) != "try-error"){
        est <- summary(fit.beta)$coefficients$cond[2,1]
        std.err <- summary(fit.beta)$coefficients$cond[2,2]
        ci <- c(est + qnorm(c(0.025, 0.975))*std.err)
        out.beta <- c(summary(fit.beta)$coefficients$cond[2,c(1,2)], NA, ci, summary(fit.beta)$coefficients$cond[2,4])
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

all.taxa.con.glmm.cov.beta <- function(con.var, id.var, cov.var, taxa, multi.method = "BH"){
  tax.out <- list()
  for (i in 1:6) {
    prop.taxa.test.out <- taxa.con.glmm.cov.beta.func(con.var, id.var, cov.var, taxa[[i]])
    prop.taxa.test.q.out <- bin.q.func(prop.taxa.test.out, method = multi.method)
    tax.out[[i]] <-  prop.taxa.test.q.out
  }
  names(tax.out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(tax.out)
}

taxa.con.glmm.cov.beta.func <- function(con.var, id.var, cov.var, taxa) {
  if(!is.null(taxa)){
    n.tax <- ncol(taxa)
    beta.out <- matrix(NA, n.tax, 6)
    for (i in 1:n.tax) {
      taxon <- taxa[,i]
      dat <- as.data.frame(cbind(con.var, id.var, cov.var, taxon))
      f <- formula(paste(colnames(dat)[ncol(dat)], " ~ ", colnames(dat)[1], "+", paste(names(cov.var), collapse = "+"), "+ (1|id.var)", sep = ""))
      fit.beta <- try(glmmTMB(f, data = dat, family = beta_family()), silent = TRUE)
      if(class(fit.beta) != "try-error"){
        est <- summary(fit.beta)$coefficients$cond[2,1]
        std.err <- summary(fit.beta)$coefficients$cond[2,2]
        ci <- c(est + qnorm(c(0.025, 0.975))*std.err)
        out.beta <- c(summary(fit.beta)$coefficients$cond[2,c(1,2)], NA, ci, summary(fit.beta)$coefficients$cond[2,4])
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

q.func <- function(out, method = c("BH", "BY")) {
  if(is.null(out)){
    return(NULL)
  }else{
    Q.value <- p.adjust(out$P.value, method = method)
    return(cbind(out, Q.value))
  }
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  return(x)
}

clr.lmer.table <- function(all.taxa.bin.lmer.out) {
  test <- list()
  for (i in 1: length(all.taxa.bin.lmer.out)) {
    test[[i]] <- all.taxa.bin.lmer.out[[i]]["clr.lmer.out"] 
  }
  names(test) <- names(all.taxa.bin.lmer.out)
  return(test)
}

clr.sum.table <- function(all.taxa.bin.lmer.out) {
  test <- list()
  for (i in 1: length(all.taxa.bin.lmer.out)) {
    test[[i]] <- all.taxa.bin.lmer.out[[i]]["clr.sum"]
  }
  names(test) <- names(all.taxa.bin.lmer.out)
  return(test)
}
