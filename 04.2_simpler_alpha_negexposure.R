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
  simpler <- fread("work/simpler_working_dataset.csv", na.strings = c("NA", NA, ""))
  load('work/simpler_model.Rdata')
  
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
  
  run.fun <- function(model = main.model, cohort = simpler){
    resclass_shannon <- linear.fun(y="shannon", exp = atbclasses, data = cohort, model = model)
    resclass_rich <- linear.fun(y="richness", exp = atbclasses, data = cohort, model= model )
    resclass_invsimpson <- linear.fun(y="invsimpson", exp = atbclasses, data = cohort, model= model )
    
    res <- rbind(resclass_shannon, resclass_rich, resclass_invsimpson)
    
    res <- res %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  
  

  ressimpler_dis <- run.fun(model = full.model, cohort = simpler)
  ressimpler_dis$model <- "full.model"
  ressimpler_dis$cohort <- "SIMPLER"
  
  
  res.alpha = rbind(ressimpler_dis)
  
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
  
  
  ressimpler_dis_atb0 <- run.fun(model = full.model, cohort = simpler2)
  ressimpler_dis_atb0$model <- "full.model"
  ressimpler_dis_atb0$cohort <- "SIMPLER"
  
  
  cc <- complete.cases(simpler2[, full.model , with = F])
  Nexposed <- simpler2[cc, .(exposure = colnames(.SD), Nexposed = lapply(.SD, function(x) sum(x>=1))), .SDcols = atbclasses]
  
  ressimpler_dis_atb0 <- merge(ressimpler_dis_atb0, Nexposed , by = "exposure")
  
  fwrite(ressimpler_dis_atb0, file='results/res_simpler_alpha_negexposure_preatbzero.tsv')
  
  message("End negative exposure 2")