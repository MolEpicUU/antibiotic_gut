# Project Antibiotic Use and the Gut Microbiota

# MOS ------------------------------------------------------------------------

rm(list=ls())

# Load packages
  library(data.table)
  library(lmerTest)
    library(dplyr)
    library(tidyr)
    library(car)

  setwd('nobackup/users/baldanzi/atb_gut/')

  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  registerdata  <-  fread('nobackup/users/baldanzi/atb_gut/Data/MOS/MOS_Halftime_n2644_Prescribed_Drug_register_until_2019.tsv')
  load('work/mos_model_revision.Rdata')
  
  # Filter register data
  registerdata <- registerdata[grep("^J01", ATC),]
  registerdata <- registerdata[lopnrMOS %in% mos$lopnrMOS, ]
  registerdata <- merge(mos[, .(lopnrMOS, Visit1, Visit2)], registerdata, by="lopnrMOS")
  registerdata <- registerdata[edatum>as.Date(Visit1) & as.Date(edatum)<(Visit2+365.25), ]
  
  registerdata_after <- registerdata[edatum>Visit1 & edatum<(Visit2+365.25), .(lopnrMOS, ATC, edatum, Visit1, Visit2)]
  
  registerdata_after[grep("^J0",ATC), class:="other"]
  registerdata_after[grep("^J01D[B,C,D,E]",ATC), class:="Class_cephalosporins"]
  registerdata_after[grep("^J01FA",ATC), class:="Class_macrolides"]
  registerdata_after[grep("^J01FF",ATC), class:="Class_lincosamides"]
  registerdata_after[grep("^J01MA",ATC), class:="Class_FQs"]
  registerdata_after[grep("^J01A",ATC), class:="Class_TCLs"]
  registerdata_after[grep("^J01E",ATC), class:="Class_SMZTMP"]
  registerdata_after[grep("^J01XE",ATC), class:="Class_NIT"]
  registerdata_after[grep("^J01CA",ATC), class:="Class_Peni_Ext"]
  registerdata_after[grep("^J01CE",ATC), class:="Class_Peni_BetaS"]
  registerdata_after[grep("^J01CF",ATC), class:="Class_Peni_BetaR"]
  registerdata_after[grep("^J01CR",ATC), class:="Class_Peni_Comb"]
  
  N_atb_class <- registerdata_after[, .N, by=.(lopnrMOS,class)]
  
  ind_zero_atb <- mos$lopnrMOS[!mos$lopnrMOS %in% N_atb_class$lopnrMOS]
  
  N_atb_class <- rbind(N_atb_class, data.frame(lopnrMOS = ind_zero_atb), fill=T)
  N_atb_class[is.na(N), N:=0]
  
  N_atb_class[ , Nafter1yr := sum(N) , by = lopnrMOS]
  N_atb_class[ , class := paste0(class,"_after1yr")]
  
  N_atb_class <-  N_atb_class %>% pivot_wider(id_cols = c(lopnrMOS,Nafter1yr), names_from = class,
                                              values_from = N, values_fill = 0)
  
  N_atb_class$`NA` <- NULL
  N_atb_class$NA_after1yr <- NULL
  
  mos <- merge(mos, N_atb_class, by ="lopnrMOS")
  
  
  # Association between number of antibiotics and alpha diversity  
  
  mosatbclasses <-  grep("Class_.*yr", colnames(mos), value=T)
  
  mix.fun <- function(y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), model = basic.model, data){
    
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
  
  
  run.fun <- function(model = basic.model, cohort = mos){
    resmosclass_shannon <- mix.fun(y="shannon", exp = mosatbclasses, data = cohort, model = model)
    resmosclass_rich <- mix.fun(y="richness", exp = mosatbclasses, data = cohort, model= model )
    resmosclass_invsimpson <- mix.fun(y="invsimpson", exp = mosatbclasses, data = cohort, model= model )
  
    resmos <- rbind(resmosclass_shannon, resmosclass_rich, resmosclass_invsimpson)
  
    res <- resmos %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  
  resmos_full <- run.fun(model = full.model)
  resmos_full$model <- "full.model"
  

  res.alpha = rbind(resmos_full)
  cc <- complete.cases(mos[, full.model , with = F])
  Nexposed <- mos[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = mosatbclasses]
  
  res.alpha <- merge(res.alpha, Nexposed , by = "exposure")
  
  #Save results 
  fwrite(res.alpha, file='results/res_mos_alpha_negexposure.tsv')
  
  message("End negative exposure 1")
  
  # Negative exposure no previous antibiotic  ----------------------------------------------------------
  
  mos2 <- mos[N1yr ==0 & N1_4yr == 0 & N4_8yr == 0, ]
  mosatbclasses = grep("Class_.*after1yr", colnames(mos2), value=T)
  
  mix.fun <- function(y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), model = basic.model, data){
    
    form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exp, collapse="+"),"+(1|family)"))
    
    fit <- lmer(form,data)
    
    temp.res <- summary(fit)
    
    N = length(temp.res$residuals)
    res <- temp.res$coef[exp, ]
    colnames(res) <- c("beta", "SE","df", "t.value","p.value")
    ci <- suppressMessages(confint(fit,parm=exp))
    lci = data.frame(LCI = ci[,1])
    hci = data.frame(HCI = ci[,2])  
    
    data.frame(outcome = y, exposure = rownames(res), res, lci=lci,hci=hci, N = N, cohort = "MOS")
    
  }
  
  
  resmos_full_atb0 <- run.fun(model = full.model, cohort = mos2)
  resmos_full_atb0$model <- "full.model"
  
  
  cc <- complete.cases(mos2[, full.model , with = F])
  Nexposed <- mos2[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = mosatbclasses]
  
  resmos_full_atb0 <- merge(resmos_full_atb0, Nexposed , by = "exposure")
  
  fwrite(resmos_full_atb0, file='results/res_mos_alpha_negexposure_preatbzero.tsv')
  
  message("End negative exposure 2")
