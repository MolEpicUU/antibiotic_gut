# Project - Antibiotic use and the gut microbiota


# This script will compare the betas from two strata  of sex or age

rm(list=ls())
library(data.table)


setwd('results')


  # Import sens analysis results (stratified by sex or age (above or below 55 yrs))
  res_sa <- fread("meta_species__sa_stratified.tsv")
  res_sex <- res_sa[model %in% c("female", "male"), ]
  res_age_above55 <- res_sa[model %in% c("age > 55 years"), ]
  
  res_sa <- fread("meta_species__sa_stratified_withoutsimpler.tsv")
  res_age_below55 <- res_sa[model %in% c("age <= 55 years"), ]
  
  
  res_age <- rbindlist(list(res_age_below55, res_age_above55), fill=T)
  
  atbs <- unique(res_sa$exposure)
  atbs <- atbs[grep("FQs|lincosam|BetaR|BetaS|Peni_Ext|TCL", atbs)] # removing antibiotics that are not common
  
  res_sex <- res_sex[exposure %in% atbs, ]
  res_age <- res_age[exposure %in% atbs, ]
  
  
  # Heterogeneity function (Z-test for the difference between two estimates) 
  heterog.test.fun <- function(res){
    
    sig_species <- unique(res[q.value<0.05, .(exposure, outcome)])
    res <- merge(sig_species, res, by = c("exposure", "outcome"))
    res <- res[, .(exposure, outcome, model, beta, SE, q.value)]
    res <- dcast(res, formula = exposure + outcome ~ model , value.var = c("beta","SE", "q.value"))
    oldnames <- colnames(res)
    colnames(res) <- c("exposure", "outcome", "beta_1", "beta_2", "se_1", "se_2", "q.value_1", "q.value_2")
    
    diff <-   res$beta_1 - res$beta_2
    se <- sqrt((res$se_1^2)+(res$se_2^2))
    z.score <- abs(diff/se)
    pv <- pnorm(z.score, lower.tail = F)*2
    
    colnames(res) <- oldnames
    res$heterog_p.value <- pv
    res$heterog_q.value <- p.adjust(pv, method="BH")
    
    return(res)
    
  }
  
  
  res_sex <- heterog.test.fun(res_sex)
  res_age <- heterog.test.fun(res_age)
  
  fwrite(res_sex , file = "hetero_beta_sex.tsv", sep="\t")
  fwrite(res_age , file = "hetero_beta_age.tsv", sep="\t")
  
  