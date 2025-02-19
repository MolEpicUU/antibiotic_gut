# Project Antibiotic Use and the Gut Microbiota

# Load packages
  library(data.table)
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(car)})

# This script will investigate the associations of number of previous antibiotics with alpha diversity


  # Import data set
  simpler <- fread("work/simpler_working_dataset.csv", na.strings = c("NA", NA, ""))

  # Association between number of antibiotics and alpha diversity  
  
  linear.fun <- function(model,y, exp = c("N1yr","N1_4yr","N4_8yr"), data = simpler){
  
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
  
  
  # SIMPLER ---------------------------------------------------------------------
  
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  run.fun <- function(model = basic.model, cohort = simpler){
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
  
  
  load('work/simpler_model.Rdata')
  ressimpler <- run.fun(model = basic.model, cohort = simpler)
  ressimpler$model <- "basic.model"
  ressimpler$cohort <- "SIMPLER"  
  ressimpler_dis <- run.fun(model = full.model, cohort = simpler)
  ressimpler_dis$model <- "full.model"
  ressimpler_dis$cohort <- "SIMPLER"
  
  
  res.alpha = rbind(ressimpler, ressimpler_dis)
  
  
  fwrite(res.alpha, file='results/res_simpler_alpha.tsv')

  
  # Sensitivity analysis exclusion of hospitalized participants _____
  
  ressimpler_sa1 <- run.fun(model = full.model, cohort = simpler[hospinfect=="no",])
  ressimpler_sa1$model <- "full.model_hospitalizedinfect"
  ressimpler_sa1$cohort <- "SIMPLER"
  
  ressimpler_sa2 <- run.fun(model = full.model, cohort = simpler[hospgeneral=="no",])
  ressimpler_sa2$model <- "full.model_hospitalizedgeneral"
  ressimpler_sa2$cohort <- "SIMPLER"
  
  res.alpha_sa = rbind(ressimpler_sa1, ressimpler_sa2)
  
  fwrite(res.alpha_sa, file='results/res_simpler_alpha_sensitivityanalyses.tsv')
  
  
  message("End")
  
 