# Project Antibiotic Use and the Gut Microbiota

rm(list=ls())

# Load packages
  library(data.table)
    library(dplyr)
    library(tidyr)
    library(car)

  setwd('users/baldanzi/')

# This script will investigate the associations of number of previous antibiotics with alpha diversity


  # Import data set
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  registerdata   <- fread('data/phenotypes/prescibed_drugs.csv', na.strings = c("NA",NA,"")) 
  load('work/simpler_model_revision.Rdata')
  
  # Filter register data
  registerdata <- registerdata[grep("^J01", ATC),]
  registerdata <- merge(simpler[, .(SIMPKEY, Visit1)], registerdata, by="SIMPKEY")
  
  # Atb after 
  registerdata_after <- registerdata[EDATUM>Visit1 & EDATUM<(Visit1+365.25), .(SIMPKEY, ATC, EDATUM, Visit1)]
  
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
  
  N_atb_class <- registerdata_after[, .N, by=.(SIMPKEY,class)]
  ind_zero_atb <- simpler$SIMPKEY[!simpler$SIMPKEY %in% N_atb_class$SIMPKEY]
  
  N_atb_class <- rbind(N_atb_class, data.frame(SIMPKEY = ind_zero_atb), fill=T)
  N_atb_class[is.na(N), N:=0]
  
  N_atb_class[ , Nafter1yr := sum(N) , by = SIMPKEY]
  N_atb_class[ , class := paste0(class,"_after1yr")]
  
  N_atb_class <-  N_atb_class %>% pivot_wider(id_cols = c(SIMPKEY,Nafter1yr), names_from = class,
                                              values_from = N, values_fill = 0)
  
  N_atb_class$`NA` <- NULL
  N_atb_class$NA_after1yr <- NULL
  
  simpler <- merge(simpler, N_atb_class, by ="SIMPKEY")


  # Association between number of antibiotics and alpha diversity  
  
  linear.fun <- function(model,y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), data = simpler){
  
    form <- as.formula(paste0(y,"~", paste0(model, collapse = "+"), "+", paste0(exp,collapse="+")))
  
    fit <- lm(form,data = data)
    ci <- confint(fit)
    LCI <- data.frame(LCI=ci[exp,1])
    HCI <- data.frame(HCI=ci[exp,2])
    temp.res <- summary(fit)
    N = length(temp.res$residuals)
    res <- temp.res$coef[exp, ]
    colnames(res) <- c("beta", "SE", "t.value","p.value")

    res <- data.frame(outcome = y, exposure = rownames(res), res, N = N, LCI = LCI, HCI = HCI)
    return(res)
  }
  
  
  # SIMPLER ---------------------------------------------------------------------
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  run.fun <- function(model = basic.model, cohort = simpler){
    resclass_shannon <- linear.fun(y="shannon", exp = atbclasses, data = cohort, model = model)
    resclass_rich <- linear.fun(y="richness", exp = atbclasses, data = cohort, model= model )
    resclass_invsimpson <- linear.fun(y="invsimpson", exp = atbclasses, data = cohort, model= model )
    
    res <- rbind(resclass_shannon, resclass_rich, resclass_invsimpson)
    
    res <- res %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  
  
  ressimpler_full <- run.fun(model = full.model, cohort = simpler)
  ressimpler_full$model <- "full.model"
  ressimpler_full$cohort <- "SIMPLER"
  
  
  res.alpha = rbind(ressimpler_full)
  
  cc <- complete.cases(simpler[, full.model , with = F])
  Nexposed <- simpler[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = atbclasses]
  
  res.alpha <- merge(res.alpha, Nexposed , by = "exposure")
  
  fwrite(res.alpha, file='results/res_simpler_alpha_negexposure.tsv')
  
  message("End negative exposure 1")
  
  # Negative exposure no previous antibiotic  ----------------------------------------------------------
  
  simpler2 <- simpler[N1yr ==0 & N1_4yr == 0 & N4_8yr == 0, ] #1306
  atbclasses = grep("Class_.*after1yr", colnames(simpler2), value=T)
  
  linear.fun <- function(model,y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), data = simpler){
    
    form <- as.formula(paste0(y,"~", paste0(model, collapse = "+"), "+", paste0(exp,collapse="+")))
    
    fit <- lm(form,data = data)
    ci <- confint(fit)
    LCI <- data.frame(LCI=ci[exp,1])
    HCI <- data.frame(HCI=ci[exp,2])
    temp.res <- summary(fit)
    N = length(temp.res$residuals)
    res <- temp.res$coef[exp, ]
    colnames(res) <- c("beta", "SE", "t.value","p.value")
    
    res <- data.frame(outcome = y, exposure = rownames(res), res, N = N, LCI = LCI, HCI = HCI)
    return(res)
  }
  
  
  ressimpler_full_atb0 <- run.fun(model = full.model, cohort = simpler2)
  ressimpler_full_atb0$model <- "full.model"
  ressimpler_full_atb0$cohort <- "SIMPLER"
  
  
  cc <- complete.cases(simpler2[, full.model , with = F])
  Nexposed <- simpler2[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = atbclasses]
  
  ressimpler_full_atb0 <- merge(ressimpler_full_atb0, Nexposed , by = "exposure")
  
  fwrite(ressimpler_full_atb0, file='results/res_simpler_alpha_negexposure_preatbzero.tsv')
  
  message("End negative exposure 2")