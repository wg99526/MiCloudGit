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
library(broom.mixed)
library(gee)
library(geepack)

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

alpha.bin.var.func <- function(sam.dat) {
  var.names <- colnames(sam.dat)
  return(var.names)
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

alpha.bin.sum.func <- function(bin.var, alpha.div) {
  n.alpha <- ncol(alpha.div)
  ref.sum <- matrix(NA, n.alpha, 7)
  com.sum <- matrix(NA, n.alpha, 7)
  for (i in 1:n.alpha) {
    ind.alpha <- alpha.div[,i]
    sum.out <- tapply(ind.alpha, bin.var, alpha.ind.sum.func)
    ref.sum[i,] <- sum.out[[1]]
    com.sum[i,] <- sum.out[[2]]
  }
  rownames(ref.sum) <- colnames(alpha.div)
  colnames(ref.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  rownames(com.sum) <- colnames(alpha.div)
  colnames(com.sum) <- c("N", "Mean", "Minimum", "1st quartile", "Median", "3rd quartile", "Maximum")
  out <- list(as.data.frame(ref.sum), as.data.frame(com.sum))
  names(out) <- levels(bin.var)
  return(out)
}

alpha.bin.id.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sel.id.var, sam.dat, alpha.div) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  id.var <- sam.dat[,sel.id.var]
  id.var <- as.data.frame(rbind(id.var[ind.ref,], id.var[ind.com,]))
  alpha.div <- rbind(alpha.div[ind.ref,], alpha.div[ind.com,])
  
  return(list(bin.var = bin.var, id.var = id.var, alpha.div = alpha.div))
}

alpha.lmer.bin.id.func <- function(bin.var, id.var, alpha.div, scale = TRUE) {
  
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, id.var, alpha.div))
  
  out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste("scale(", alpha.ind[i], ")", "~", colnames(bin.var), "+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    if (!scale) {
      f <- formula(paste(alpha.ind[i], "~", colnames(bin.var), "+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    fit <- try(lmer(f, data = d), silent = TRUE)
    est.se.df <- tryCatch(summary(fit)$coefficients[2,c(1,2,3)], error = function(err) c(NA,NA,NA))
    ci <- tryCatch(confint(fit)[4,], error = function(err) c(NA,NA))
    pval <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
    out[i,] <- c(est.se.df, ci, pval)
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha.div)
  colnames(out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(out)
}

alpha.forest.lmer.plot <- function(out, mult.test.cor = TRUE) {
  
  if (mult.test.cor) {
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "Est", "SE", "P-value", "Q-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), p.value.0.1(out[,6]), p.value.0.1(out[,7]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=0, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }else{
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "Est", "SE", "P-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), p.value.0.1(out[,6]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=0, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }
}

alpha.forest.lmer.or.plot <- function(out, mult.test.cor = TRUE) {
  
  if (mult.test.cor) {
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "OR", "SE", "P-value", "Q-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), p.value.0.1(out[,6]), p.value.0.1(out[,7]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=1, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }
  else{
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "OR", "SE", "P-value"), 
                                    cbind(rownames(out), format(round(out[, c(1, 2)], digits = 3), nsmall = 3), p.value.0.1(out[,6]))))
    ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,1], out[,c(4,5)])))
    
    forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
               hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=1, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval",
               txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                              ticks=gpar(fontfamily="", cex=0.7),
                              xlab=gpar(fontfamily="", cex=0.7)))
  }
}

#######################
# Binary - Covariates #
#######################

alpha.bin.id.cov.cat.ref.func <- function(sel.bin.var, sel.ref, sel.com, sel.id.var, sel.cov.var, sam.dat, alpha.div) {  
  bin.var <- sam.dat[,sel.bin.var]
  ind.ref <- which(bin.var == sel.ref)
  ind.com <- which(bin.var == sel.com)
  bin.var <- as.data.frame(c(rep(0, length(ind.ref)), rep(1, length(ind.com))))
  colnames(bin.var) <- sel.bin.var
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  id.var <- sam.dat[,sel.id.var]
  id.var <- as.data.frame(rbind(id.var[ind.ref,], id.var[ind.com,]))
  alpha.div <- rbind(alpha.div[ind.ref,], alpha.div[ind.com,])
  
  return(list(bin.var = bin.var, id.var = id.var, cov.var = cov.var, alpha.div = alpha.div))
}

alpha.lmer.bin.id.cov.func <- function(bin.var, id.var, cov.var, alpha.div, scale = TRUE) {
  
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, id.var, cov.var, alpha.div))
  
  out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste("scale(", alpha.ind[i], ")", "~", colnames(bin.var), "+", paste(colnames(cov.var), collapse = "+"), "+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    if (!scale) {
      f <- formula(paste(alpha.ind[i], "~", colnames(bin.var), "+", paste(colnames(cov.var), collapse = "+"), "+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    fit <- try(lmer(f, data = d), silent = TRUE)
    est.se.df <- tryCatch(summary(fit)$coefficients[2,c(1,2,3)], error = function(err) c(NA,NA,NA))
    ci <- tryCatch(confint(fit)[4,], error = function(err) c(NA,NA))
    pval <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
    out[i,] <- c(est.se.df, ci, pval)
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha.div)
  colnames(out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(out)
}

get.or.se.lmer <- function(model) {
  t <- broom::tidy(model)
  or = exp(t$estimate[c(1,2)])
  var.diag = diag(vcov(model))
  or.se = sqrt(or^2 * var.diag)
  
  return(or.se = unname(unlist(or.se)))
}

alpha.logit.reg.coef.bin.gee.func <- function(bin.var, id.var, alpha.div, scale = TRUE) {
  
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, alpha.div, id.var))
  logit.out <- matrix(NA, n.alpha, 6)
  
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste(colnames(bin.var), "~", "scale(", alpha.ind[i], ")"))
    }
    if (!scale) {
      f <- formula(paste(colnames(bin.var), "~",  alpha.ind[i]))
    }
    
    fit <- try(gee(f, id = as.factor(id), data = d, corstr="exchangeable", family = "binomial"(link = "logit")), silent = TRUE)
    
    est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
    se <- tryCatch(summary(fit)$coefficients[2,"Robust S.E."], error = function(err) NA)      
    ci <- tryCatch(c(est - qnorm(0.975)*se, est + qnorm(0.975)*se), error = function(err) c(NA, NA))
    z <- tryCatch(summary(fit)$coefficients[2,"Robust z"], error = function(err) NA)
    pvs <- tryCatch(2*pnorm(-abs(z), mean = 0, sd = 1), error = function(err) NA)
    
    out <- c(est, se, z, ci, pvs)
    logit.out[i,] <- out
  }
  logit.out <- as.data.frame(logit.out)
  rownames(logit.out) <- colnames(alpha.div)
  colnames(logit.out) <- c("Est", "Std Err", "z value", "Lower", "Upper", "P.value")
  
  return(logit.out)
}

alpha.logit.reg.coef.bin.cov.gee.func <- function(bin.var, id.var, cov.var, alpha.div, scale = TRUE) {
  
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, cov.var, alpha.div, id.var))
  logit.out <- matrix(NA, n.alpha, 6)
  
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste(colnames(bin.var), "~", "scale(", alpha.ind[i], ")", "+", paste(colnames(cov.var), collapse = "+")))
    }
    if (!scale) {
      f <- formula(paste(colnames(bin.var), "~",  alpha.ind[i], "+", paste(colnames(cov.var), collapse = "+")))
    }
    
    fit <- try(gee(f, id = as.factor(id), data = d, corstr="exchangeable", family = "binomial"(link = "logit")), silent = TRUE)
    
    est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
    se <- tryCatch(summary(fit)$coefficients[2,"Robust S.E."], error = function(err) NA)      
    ci <- tryCatch(c(est - qnorm(0.975)*se, est + qnorm(0.975)*se), error = function(err) c(NA, NA))
    z <- tryCatch(summary(fit)$coefficients[2,"Robust z"], error = function(err) NA)
    pvs <- tryCatch(2*pnorm(-abs(z), mean = 0, sd = 1), error = function(err) NA)
    
    out <- c(est, se, z, ci, pvs)
    logit.out[i,] <- out
  }
  logit.out <- as.data.frame(logit.out)
  rownames(logit.out) <- colnames(alpha.div)
  colnames(logit.out) <- c("Est", "Std Err", "z value", "Lower", "Upper", "P.value")
  
  return(logit.out)
}

alpha.logit.reg.coef.bin.glmm.b.func <- function(bin.var, id.var, alpha.div, scale = TRUE) {
  
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(unlist(bin.var), alpha.div, unlist(id.var)))
  logit.out <- matrix(NA, n.alpha, 6)
  
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste(colnames(d)[1], "~", "scale(", alpha.ind[i], ")", "+ (1|", colnames(d)[3],")"))
    }
    if (!scale) {
      f <- formula(paste(colnames(d)[1], "~",  alpha.ind[i], "+ (1|", colnames(d)[3],")"))
    }
    
    fit <- try(glmer(f, data = d, family = "binomial"(link = "logit")), silent = TRUE)
    
    est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
    se <- tryCatch(summary(fit)$coefficients[2,"Std. Error"], error = function(err) NA)      
    ci <- tryCatch(c(est - qnorm(0.975)*se, est + qnorm(0.975)*se), error = function(err) c(NA,NA))
    z <- tryCatch(summary(fit)$coefficients[2,"z value"], error = function(err) NA)
    pvs <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
    
    out <- c(est, se, z, ci, pvs)
    logit.out[i,] <- out
  }
  logit.out <- as.data.frame(logit.out)
  rownames(logit.out) <- colnames(alpha.div)
  colnames(logit.out) <- c("Est", "Std Err", "z value", "Lower", "Upper", "P.value")
  
  return(logit.out)
}

alpha.logit.reg.coef.bin.cov.glmm.b.func <- function(bin.var, id.var, cov.var, alpha.div, scale = TRUE) {
  
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(unlist(bin.var), alpha.div, unlist(id.var), cov.var))
  logit.out <- matrix(NA, n.alpha, 6)
  
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste(colnames(d)[1], "~", "scale(", alpha.ind[i], ")", "+", paste(colnames(cov.var), collapse = "+"), "+ (1|", colnames(d)[3],")"))
    }
    if (!scale) {
      f <- formula(paste(colnames(d)[1], "~",  alpha.ind[i], "+", paste(colnames(cov.var), collapse = "+"), "+ (1|", colnames(d)[3],")"))
    }
    
    fit <- try(glmer(f, data = d, family = "binomial"(link = "logit")), silent = TRUE)
    
    est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
    se <- tryCatch(summary(fit)$coefficients[2,"Std. Error"], error = function(err) NA)   
    ci <- tryCatch(c(est - qnorm(0.975)*se, est + qnorm(0.975)*se), error = function(err) c(NA,NA))
    z <- tryCatch(summary(fit)$coefficients[2,"z value"], error = function(err) NA)
    pvs <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
    
    out <- c(est, se, z, ci, pvs)
    logit.out[i,] <- out
  }
  logit.out <- as.data.frame(logit.out)
  rownames(logit.out) <- colnames(alpha.div)
  colnames(logit.out) <- c("Est", "Std Err", "z value", "Lower", "Upper", "P.value")
  
  return(logit.out)
}

alpha.logit.bin.cov.glmm.b.func <- function(bin.var, id.var, cov.var, alpha.div, scale = TRUE) {
  
  n.cov <- ncol(cov.var)
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(bin.var, cov.var, alpha.div, id.var))
  logit.out <- matrix(NA, n.alpha, 6)
  
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste(colnames(bin.var), "~", "scale(", alpha.ind[i], ")", "+", paste(colnames(cov.var), collapse = "+"),"+ (1|id.var)", sep = ""))
    }
    if (!scale) {
      f <- formula(paste(colnames(bin.var), "~", alpha.ind[i], "+", paste(colnames(cov.var), collapse = "+"),"+ (1|id.var)", sep = ""))
    }
    fit <- try(glmer(f, id = as.factor(id), data = d, corstr="exchangeable", family = "binomial"(link = "logit")), silent = TRUE)
    
    est <- tryCatch(summary(fit)$coefficients[2,"Estimate"], error = function(err) NA)
    se <- tryCatch(summary(fit)$coefficients[2,"Std.err"], error = function(err) NA)
    ci <- tryCatch(exp(est + qnorm(c(0.025, 0.975))*se), error = function(err) c(NA,NA))
    wald <- tryCatch(round(summary(fit)$coefficients[2,"Wald"], digits = 3), error = function(err) NA)
    pvs <- tryCatch(summary(fit)$coefficients[2,4], error = function(err) NA)
    
    or.se <- tryCatch(sqrt(exp(summary(fit)$coefficient[,1])^2 *diag(vcov(fit)))[2], error = function(err) NA)
    out.logit <-c(exp(est), or.se, wald, ci, pvs)
    logit.out[i,] <- out.logit
  }

  lmer.out <- as.data.frame(lmer.out)
  rownames(lmer.out) <- colnames(taxa)
  colnames(lmer.out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(lmer.out)
}

#############################
# Continuous - No covariate #
#############################

alpha.con.id.recode.func <- function(sam.dat, sel.con.var, sel.id.var, rename.con.var, alpha.div) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  id.var <- as.data.frame(as.matrix(sam.dat[ind.nona, sel.id.var]))
  colnames(id.var) <- sel.id.var
  alpha.div = alpha.div[ind.nona,]
  return(list(con.var = con.var, id.var = id.var, alpha.div = alpha.div))
}

alpha.lmer.con.id.func <- function(con.var, id.var, alpha.div, scale = TRUE) {
  
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(con.var, id.var, alpha.div))
  
  out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste("scale(", alpha.ind[i], ")", "~", colnames(con.var), "+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    if (!scale) {
      f <- formula(paste(alpha.ind[i], "~", colnames(con.var), "+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    fit <- try(lmer(f, data = d), silent = TRUE)
    est.se.df <- tryCatch(summary(fit)$coefficients[2,c(1,2,3)], error = function(err) c(NA,NA,NA))
    ci <- tryCatch(confint(fit)[4,], error = function(err) c(NA,NA))
    pval <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
    out[i,] <- c(est.se.df, ci, pval)
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha.div)
  colnames(out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(out)
}

###############################
# Continuous - with covariate #
###############################

alpha.con.id.cov.recode.func <- function(sam.dat, sel.con.var, sel.id.var, sel.cov.var, rename.con.var, alpha.div) {
  ind.sel <- which(colnames(sam.dat) == sel.con.var)
  colnames(sam.dat)[ind.sel] <- rename.con.var
  ind.nona <- !is.na(sam.dat[,ind.sel])
  con.var <- as.data.frame(as.matrix(sam.dat[ind.nona,ind.sel]))
  cov.var <- as.data.frame(sam.dat[ind.nona, sel.cov.var])
  id.var <- as.data.frame(as.matrix(sam.dat[ind.nona, sel.id.var]))
  alpha.div = alpha.div[ind.nona,]
  return(list(con.var = con.var, id.var = id.var, cov.var = cov.var, alpha.div = alpha.div))
}

alpha.lmer.con.id.cov.func <- function(con.var, cov.var, id.var, alpha.div, scale = TRUE) {
  n.alpha <- ncol(alpha.div)
  alpha.ind <- colnames(alpha.div)
  d <- as.data.frame(cbind(con.var, cov.var, id.var, alpha.div))
  
  out <- matrix(NA, n.alpha, 6)
  for (i in 1:n.alpha) {
    if (scale) {
      f <- formula(paste("scale(", alpha.ind[i], ")", "~", colnames(con.var), "+", paste(colnames(cov.var), collapse = "+") ,"+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    if (!scale) {
      f <- formula(paste(alpha.ind[i], "~", colnames(con.var), "+", paste(colnames(cov.var), collapse = "+") , "+", "(1 | ", paste(colnames(id.var)), ")"))
    }
    fit <- try(lmer(f, data = d), silent = TRUE)
    est.se.df <- tryCatch(summary(fit)$coefficients[2,c(1,2,3)], error = function(err) c(NA,NA,NA))
    ci <- tryCatch(confint(fit)[4,], error = function(err) c(NA,NA))
    pval <- tryCatch(summary(fit)$coefficients[2,5], error = function(err) NA)
    out[i,] <- c(est.se.df, ci, pval)
    
  }
  out <- as.data.frame(out)
  rownames(out) <- colnames(alpha.div)
  colnames(out) <- c("Est", "Std Err", "DF", "Lower", "Upper", "P.value")
  return(out)
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
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  return(x)
}
