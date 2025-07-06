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
  
  load('work/simpler_model_revision.Rdata')

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
    resclass_shannon <- linear.fun(y="shannon", exp = atbclasses, data = cohort, model = model)
    resclass_rich <- linear.fun(y="richness", exp = atbclasses, data = cohort, model= model )
    resclass_invsimpson <- linear.fun(y="invsimpson", exp = atbclasses, data = cohort, model= model )
    
    res <- rbind(resclass_shannon, resclass_rich, resclass_invsimpson)
    
    res <- res %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  
  

  ressimpler <- run.fun(model = basic.model, cohort = simpler)
  ressimpler$model <- "basic.model"
  ressimpler$cohort <- "SIMPLER"  
  ressimpler_full <- run.fun(model = full.model, cohort = simpler)
  ressimpler_full$model <- "full.model"
  ressimpler_full$cohort <- "SIMPLER"
  
  
  res.alpha = rbind(ressimpler, ressimpler_full)
  
  
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
  
  
  # Sensitivity analysis using different time windows of recent antibiotic use 
  # Import data set
  simpler <- fread("work/simpler_working_dataset_revision_noatbexclusion.tsv", na.strings = c("NA", NA, ""))
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  
  ressimpler_sa3_none <- run.fun(model = full.model,  cohort = simpler)
  ressimpler_sa3_none$time_windown = "none"
  
  ressimpler_sa3_30days <- run.fun(model = full.model,  cohort = simpler[last_EDATUM<Visit1-30 | is.na(last_EDATUM),, ])
  ressimpler_sa3_30days$time_windown = "30 days"
  
  ressimpler_sa3_180days <- run.fun(model = full.model,  cohort = simpler[last_EDATUM<Visit1-180 | is.na(last_EDATUM),, ])
  ressimpler_sa3_180days$time_windown = "180 days"
  
  atbclasses <- atbclasses[!grepl("_1yr", atbclasses)]
  ressimpler_sa3_1y <- run.fun(model = full.model,  cohort = simpler[N1yr == 0, ])
  ressimpler_sa3_1y$time_windown = "1 year"
  
  ressimpler_sa3 <- rbind(ressimpler_sa3_none, ressimpler_sa3_30days, ressimpler_sa3_180days, ressimpler_sa3_1y)
  ressimpler_sa3$cohort <- "SIMPLER"
  
  fwrite(ressimpler_sa3, file='results/res_simpler_alpha_sensitivityanalyses_timewindows.tsv', sep = "\t")
  
  # Sensitivity analysis stratified age and sex
  # Import data set
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  full.model_sex = full.model[-which(full.model=="sex")]
  ressimpler_male <- run.fun(model = full.model_sex,  cohort = simpler[sex == "M", ])
  ressimpler_male$model = "male"
  
  ressimpler_female <- run.fun(model = full.model_sex,  cohort = simpler[sex == "F", ])
  ressimpler_female$model = "female"
  
  # Minimal age in SIMPLER is higher than 55
  #ressimpler_agebelow55 <- run.fun(model = full.model,  cohort = simpler[age <= 55, ])
  #ressimpler_agebelow55$model = "age <= 55 years"
  
  ressimpler_ageabove55 <- run.fun(model = full.model,  cohort = simpler[age > 55, ])
  ressimpler_ageabove55$model = "age > 55 years"
  
  
  ressimpler_age_sex <- rbind(ressimpler_male, ressimpler_female, ressimpler_ageabove55)
  ressimpler_age_sex$cohort = "SIMPLER"
  
  fwrite(ressimpler_age_sex, file='results/res_simpler_alpha_sensitivityanalyses_age_sex.tsv', sep = "\t")
  
  message("End")
  
 