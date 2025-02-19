# Project Antibiotic Use and the Gut Microbiota

# Load packages
  library(data.table)
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(car)})

  # This script will investigate the associations of number of antibiotics courses with alpha diversity, including
  # courses after the fecal sample collection

  # Import data set
  scapis <- fread("work/scapis_working_dataset.tsv", na.strings = c("NA", NA, ""))
  load('work/scapis_model.Rdata')
 
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
  
  run.fun <- function(model = main.model, cohort = scapis){
    resclass_shannon <- linear.fun(y="shannon", exp = atbclasses, data = cohort, model = model)
    resclass_rich <- linear.fun(y="richness", exp = atbclasses, data = cohort, model= model )
    resclass_invsimpson <- linear.fun(y="invsimpson", exp = atbclasses, data = cohort, model= model )
    
    res <- rbind(resclass_shannon, resclass_rich, resclass_invsimpson)
    
    res <- res %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  

  resscapis_dis <- run.fun(model = full.model)
  resscapis_dis$model <- "full.model"
  resscapis_dis$cohort <- "SCAPIS"
  
  res.alpha = rbind(resscapis_dis)
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
  
  
  resscapis_dis_atb0 <- run.fun(model = full.model, cohort = scapis2)
  resscapis_dis_atb0$model <- "full.model"
  resscapis_dis_atb0$cohort <- "SCAPIS"
  
  
  cc <- complete.cases(scapis2[, full.model , with = F])
  Nexposed <- scapis2[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = atbclasses]
  
  resscapis_dis_atb0 <- merge(resscapis_dis_atb0, Nexposed , by = "exposure")
  
  fwrite(resscapis_dis_atb0, file='results/res_scapis_alpha_negexposure_preatbzero.tsv')
  
  message("End negative exposure 2")
  
  
