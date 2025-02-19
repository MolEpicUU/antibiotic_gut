# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(car)})

  setwd('atb_gut/')

# This script will investigate the associations of number of previous antibiotics with alpha diversity

  # Import data set
  scapis <- fread("work/scapis_working_dataset.tsv", na.strings = c("NA", NA, ""))
  
  load('work/scapis_model.Rdata')
 
  # Association between number of antibiotics and alpha diversity  
  
  linear.fun <- function(model,y, exp = c("N1yr","N1_4yr","N4_8yr"), data = scapis){
  
  form <- as.formula(paste0(y,"~", paste0(model, collapse = "+"), "+", paste0(exp,collapse="+")))
  
  fit <- lm(form,data = data)
  ci <- confint(fit)
  gvif <- car::vif(fit)[exp,1]
  LCI <- data.frame(LCI=ci[exp,1])
  HCI <- data.frame(HCI=ci[exp,2])
  temp.res <- summary(fit)
  N = length(temp.res$residuals)
  res <- temp.res$coef[exp, ]
  colnames(res) <- c("beta", "SE", "t.value","p.value")

  res <- data.frame(outcome = y, exposure = rownames(res), res, N = N, LCI = LCI, HCI = HCI, 
                    gvif = gvif)
  return(res)
  }
  
  
  # SCAPIS ---------------------------------------------------------------------
  
  
  atbclasses = grep("Class_.*yr", colnames(scapis), value=T)
  
  run.fun <- function(model = main.model, cohort = scapis){
    resshannon <- linear.fun(y="shannon", data = cohort, model = model)
    resrich <- linear.fun(y="richness",  data = cohort, model = model)
    res_invsimposon <- linear.fun(y="invsimpson", data = cohort, model = model)
    resclass_shannon <- linear.fun(y="shannon", exp = atbclasses, data = cohort, model = model)
    resclass_rich <- linear.fun(y="richness", exp = atbclasses, data = cohort, model= model )
    resclass_invsimpson <- linear.fun(y="invsimpson", exp = atbclasses, data = cohort, model= model )
    
    res <- rbind(resshannon, resrich,res_invsimposon, resclass_shannon, resclass_rich, resclass_invsimpson)
    
    res <- res %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  
  
  resscapis <- run.fun(model = basic.model)
  resscapis$model <- "basic.model"
  resscapis$cohort <- "SCAPIS"
  resscapis_dis <- run.fun(model = full.model)
  resscapis_dis$model <- "full.model"
  resscapis_dis$cohort <- "SCAPIS"
  
  res.alpha = rbind(resscapis, resscapis_dis)
  
  
  fwrite(res.alpha, file='results/res_scapis_alpha.tsv')
  
  
  message("End linear regression")
  
  # Sensitivity analysis exclusion of hospitalized participants _____

  resscapis_sa1 <- run.fun(model = full.model, cohort = scapis[hospinfect=="no",])
  resscapis_sa1$model <- "full.model_hospitalizedinfect"
  resscapis_sa1$cohort <- "SCAPIS"
  
  resscapis_sa2 <- run.fun(model = full.model, cohort = scapis[hospgeneral=="no",])
  resscapis_sa2$model <- "full.model_hospitalizedgeneral"
  resscapis_sa2$cohort <- "SCAPIS"
  
  res.alpha_sa = rbind(resscapis_sa1, resscapis_sa2)
  
  fwrite(res.alpha_sa, file='results/res_scapis_alpha_sensitivityanalyses.tsv')
  