# Project Antibiotic Use and the Gut Microbiota

# MOS ------------------------------------------------------------------------

rm(list=ls())

# Load packages
  library(data.table)
  library(lmerTest)
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(car)})

  setwd('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/')

  # This script will investigate the associations of number of antibiotics courses with alpha diversity, including
  # courses after the fecal sample collection

  # Import data set
  mos <- fread("work/mos_working_dataset.tsv", na.strings = c("NA", NA, ""))
  load('work/mos_model.Rdata')
  
  
  # Association between number of antibiotics and alpha diversity  
  
  mosatbclasses <-  grep("Class_.*yr", colnames(mos), value=T)
  
  mix.fun <- function(y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), model = main.model, data){
    
      form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exp, collapse="+"),"+(1|family)"))
      
      fit <- lmer(form,data)
      
      temp.res <- summary(fit)
      #gvif <- car::vif(fit)[exp,1]
      
      N = length(temp.res$residuals)
      res <- temp.res$coef[exp, ]
      colnames(res) <- c("beta", "SE","df", "t.value","p.value")
      ci <- suppressMessages(confint(fit,parm=exp))
      lci = data.frame(LCI = ci[,1])
      hci = data.frame(HCI = ci[,2])  
      
      data.frame(outcome = y, exposure = rownames(res), res, lci=lci,hci=hci, N = N, cohort = "MOS")
    
  }
  
  
  run.fun <- function(model = main.model, cohort = mos){
  resmosclass_shannon <- mix.fun(y="shannon", exp = mosatbclasses, data = cohort, model = model)
  resmosclass_rich <- mix.fun(y="richness", exp = mosatbclasses, data = cohort, model= model )
  resmosclass_invsimpson <- mix.fun(y="invsimpson", exp = mosatbclasses, data = cohort, model= model )
  
  resmos <- rbind(resmosclass_shannon, resmosclass_rich, resmosclass_invsimpson)
  
  res <- resmos %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  return(res)
  }
  
  resmos_dis <- run.fun(model = full.model)
  resmos_dis$model <- "full.model"
  

  res.alpha = rbind(resmos_dis)
  cc <- complete.cases(mos[, full.model , with = F])
  Nexposed <- mos[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = mosatbclasses]
  
  res.alpha <- merge(res.alpha, Nexposed , by = "exposure")
  
  #Save results 
  fwrite(res.alpha, file='results/res_mos_alpha_negexposure.tsv')
  
  message("End negative exposure 1")
  
  # Negative exposure no previous antibiotic  ----------------------------------------------------------
  
  mos2 <- mos[N1yr ==0 & N1_4yr == 0 & N4_8yr == 0, ]
  mosatbclasses = grep("Class_.*after1yr", colnames(mos2), value=T)
  
  mix.fun <- function(y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), model = main.model, data){
    
    form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exp, collapse="+"),"+(1|family)"))
    
    fit <- lmer(form,data)
    
    temp.res <- summary(fit)
    #gvif <- car::vif(fit)[exp,1]
    
    N = length(temp.res$residuals)
    res <- temp.res$coef[exp, ]
    colnames(res) <- c("beta", "SE","df", "t.value","p.value")
    ci <- suppressMessages(confint(fit,parm=exp))
    lci = data.frame(LCI = ci[,1])
    hci = data.frame(HCI = ci[,2])  
    
    data.frame(outcome = y, exposure = rownames(res), res, lci=lci,hci=hci, N = N, cohort = "MOS")
    
  }
  
  
  resmos_dis_atb0 <- run.fun(model = full.model, cohort = mos2)
  resmos_dis_atb0$model <- "full.model"
  
  
  cc <- complete.cases(mos2[, full.model , with = F])
  Nexposed <- mos2[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = mosatbclasses]
  
  resmos_dis_atb0 <- merge(resmos_dis_atb0, Nexposed , by = "exposure")
  
  fwrite(resmos_dis_atb0, file='results/res_mos_alpha_negexposure_preatbzero.tsv')
  
  message("End negative exposure 2")
