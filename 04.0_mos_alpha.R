# Project Antibiotic Use and the Gut Microbiota

rm(list=ls())

# Load packages
  library(data.table)
  library(lmerTest)
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(car)})

  setwd('atb_gut/')

# This script will investigate the associations of number of previous antibiotics with alpha diversity


  # Import data set
  mos <- fread("work/mos_working_dataset.tsv", na.strings = c("NA", NA, ""))
  load('work/mos_model.Rdata')
  
  # Association between number of antibiotics and alpha diversity  
  

  # MOS ------------------------------------------------------------------------
  
  mosatbclasses <-  grep("Class_.*yr", colnames(mos), value=T)
  
  
  mix.fun <- function(y, exp = c("N1yr","N1_4yr","N4_8yr"), model = main.model, data){
    
      form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exp, collapse="+"),"+(1|family)"))
      
      fit <- lmer(form,data)
      
      temp.res <- summary(fit)
      gvif <- car::vif(fit)[exp,1]
      
      N = length(temp.res$residuals)
      res <- temp.res$coef[exp, ]
      colnames(res) <- c("beta", "SE","df", "t.value","p.value")
      ci <- suppressMessages(confint(fit,parm=exp))
      lci = data.frame(LCI = ci[,1])
      hci = data.frame(HCI = ci[,2])  
      
      data.frame(outcome = y, exposure = rownames(res), res, lci=lci,hci=hci, N = N, cohort = "MOS", gvif = gvif)
    
  }
  
  
  run.fun <- function(model = main.model){
  resmosshannon <- mix.fun(y="shannon",  data = mos, model = model)
  resmosrich <- mix.fun(y="richness",  data = mos, model = model)
  resmosinvsimpson <- mix.fun(y="invsimpson", data = mos, model = model)
  resmosclass_shannon <- mix.fun(y="shannon", exp = mosatbclasses, data = mos, model = model)
  resmosclass_rich <- mix.fun(y="richness", exp = mosatbclasses, data = mos, model= model )
  resmosclass_invsimpson <- mix.fun(y="invsimpson", exp = mosatbclasses, data = mos, model= model )
  
  resmos <- rbind(resmosshannon, resmosrich, resmosinvsimpson, resmosclass_shannon, resmosclass_rich, resmosclass_invsimpson)
  
  res <- resmos %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  return(res)
  }
  
  resmos <- run.fun(model = basic.model)
  resmos$model <- "basic.model"
  resmos_dis <- run.fun(model = full.model)
  resmos_dis$model <- "full.model"
  

  res.alpha = rbind(resmos, resmos_dis)
  
  #Save results 
  fwrite(res.alpha, file='results/res_mos_alpha.tsv')
  
  message("End linear regression")
  
