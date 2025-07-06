# Project Antibiotic Use and the Gut Microbiota

rm(list=ls())

# Load packages
  library(data.table)
    library(dplyr)
    library(tidyr)
    library(car)

  setwd('nobackup/users/baldanzi/atb_gut/')

# This script will investigate the associations of number of previous antibiotics with alpha diversity

  # Import data set
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  registerdata <-   fread('/proj/sens2019512/SCAPIS/Gutsy/Phenotypes/Medication/Raw/SCAPIS-REGISTERSU-1715-LMED-20220511-T_R_LMED_27947_2021.txt', na.strings = c("NA",NA,""))
  load('work/scapis_model_revision.Rdata')
  
  # Filter scapis and register data
  # Restrict scapis to the individuals with at least one year of register data follow-up after visit 2
  # The latest date in the drug register data we have is 2018-12-31
  scapis <- scapis[Visit2<"2017-12-31",]
  setnames(registerdata, "ID","Subject")
  registerdata <- registerdata[grep("^J01", ATC),]
  registerdata <- merge(scapis[, .(Subject, Visit1, Visit2)], registerdata, by="Subject")
  
  # Atb after 
  registerdata_after <- registerdata[EDATUM>Visit1 & EDATUM<(Visit2+365.25), .(Subject, ATC, EDATUM, Visit1, Visit2)]
  
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
  
  N_atb_class <- registerdata_after[, .N, by=.(Subject,class)]
  ind_zero_atb <- scapis$Subject[!scapis$Subject %in% N_atb_class$Subject]
  
  N_atb_class <- rbind(N_atb_class, data.frame(Subject = ind_zero_atb), fill=T)
  N_atb_class[is.na(N), N:=0]
  
  N_atb_class[ , Nafter1yr := sum(N) , by = Subject]
  N_atb_class[ , class := paste0(class,"_after1yr")]
  
  N_atb_class <-  N_atb_class %>% pivot_wider(id_cols = c(Subject,Nafter1yr), names_from = class,
                                              values_from = N, values_fill = 0)
  
  N_atb_class$`NA` <- NULL
  N_atb_class$NA_after1yr <- NULL
  
  scapis <- merge(scapis, N_atb_class, by ="Subject")
 
  # Association between number of antibiotics and alpha diversity  
  
  linear.fun <- function(model,y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), data = scapis){
  
  form <- as.formula(paste0(y,"~", paste0(model, collapse = "+"), "+", paste0(exp,collapse="+")))
  
  fit <- lm(form,data = data)
  ci <- confint(fit)
  #gvif <- car::vif(fit)[exp,1]
  LCI <- data.frame(LCI=ci[exp,1])
  HCI <- data.frame(HCI=ci[exp,2])
  temp.res <- summary(fit)
  N = length(temp.res$residuals)
  res <- temp.res$coef[exp, ]
  colnames(res) <- c("beta", "SE", "t.value","p.value")

  res <- data.frame(outcome = y, exposure = rownames(res), res, N = N, LCI = LCI, HCI = HCI)
  return(res)
  }
  
  
  # SCAPIS ---------------------------------------------------------------------
  
  
  atbclasses <-  grep("Class_.*yr", colnames(scapis), value=T)
  
  run.fun <- function(model = basic.model, cohort = scapis){
    resclass_shannon <- linear.fun(y="shannon", exp = atbclasses, data = cohort, model = model)
    resclass_rich <- linear.fun(y="richness", exp = atbclasses, data = cohort, model= model )
    resclass_invsimpson <- linear.fun(y="invsimpson", exp = atbclasses, data = cohort, model= model )
    
    res <- rbind(resclass_shannon, resclass_rich, resclass_invsimpson)
    
    res <- res %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  

  resscapis_full <- run.fun(model = full.model)
  resscapis_full$model <- "full.model"
  resscapis_full$cohort <- "SCAPIS"
  
  res.alpha = rbind(resscapis_full)
  cc <- complete.cases(scapis[, full.model , with = F])
  Nexposed <- scapis[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = atbclasses]
  
  res.alpha <- merge(res.alpha, Nexposed , by = "exposure")
  
  fwrite(res.alpha, file='results/res_scapis_alpha_negexposure.tsv')
  
  
  message("End negative exposure 1")
  
  # Negative exposure no previous antibiotic  ----------------------------------------------------------
  
  scapis2 <- scapis[N1yr ==0 & N1_4yr == 0 & N4_8yr == 0, ]
  atbclasses = grep("Class_.*after1yr", colnames(scapis2), value=T)
  
  linear.fun <- function(model,y, exp = c("N1yr","N1_4yr","N4_8yr", "Nafter1yr"), data = scapis){
    
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
  
  
  resscapis_full_atb0 <- run.fun(model = full.model, cohort = scapis2)
  resscapis_full_atb0$model <- "full.model"
  resscapis_full_atb0$cohort <- "SCAPIS"
  
  
  cc <- complete.cases(scapis2[, full.model , with = F])
  Nexposed <- scapis2[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = atbclasses]
  
  resscapis_full_atb0 <- merge(resscapis_full_atb0, Nexposed , by = "exposure")
  
  fwrite(resscapis_full_atb0, file='results/res_scapis_alpha_negexposure_preatbzero.tsv')
  
  message("End negative exposure 2")
  
  
