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

  load('work/scapis_model_revision.Rdata')
 
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
  
  run.fun <- function(model = basic.model, cohort = scapis){
    resclass_shannon <- linear.fun(y="shannon", exp = atbclasses, data = cohort, model = model)
    resclass_rich <- linear.fun(y="richness", exp = atbclasses, data = cohort, model= model )
    resclass_invsimpson <- linear.fun(y="invsimpson", exp = atbclasses, data = cohort, model= model )
    
    res <- rbind(resclass_shannon, resclass_rich, resclass_invsimpson)
    
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
  
  
  fwrite(res.alpha, file='results/res_scapis_alpha.tsv', , sep = "\t")
  
  
  message("End linear regression")
  
  # Sensitivity analysis exclusion of hospitalized participants _____

  resscapis_sa1 <- run.fun(model = full.model, cohort = scapis[hospinfect=="no",])
  resscapis_sa1$model <- "full.model_hospitalizedinfect"
  resscapis_sa1$cohort <- "SCAPIS"
  
  resscapis_sa2 <- run.fun(model = full.model, cohort = scapis[hospgeneral=="no",])
  resscapis_sa2$model <- "full.model_hospitalizedgeneral"
  resscapis_sa2$cohort <- "SCAPIS"
  
  res.alpha_sa = rbind(resscapis_sa1, resscapis_sa2)
  
  fwrite(res.alpha_sa, file='results/res_scapis_alpha_sensitivityanalyses.tsv', sep = "\t")
  
  # Sensitivity analysis using different time windows of recent antibiotic use 
  # Import data set
  scapis <- fread("work/scapis_working_dataset_revision_noatbexclusion.tsv", na.strings = c("NA", NA, ""))
  
  atbclasses = grep("Class_.*yr", colnames(scapis), value=T)
  
  
  resscapis_sa3_none <- run.fun(model = full.model,  cohort = scapis)
  resscapis_sa3_none$time_windown = "none"
  
  resscapis_sa3_30days <- run.fun(model = full.model,  cohort = scapis[last_EDATUM<Visit1-30 | is.na(last_EDATUM),, ])
  resscapis_sa3_30days$time_windown = "30 days"
  
  resscapis_sa3_180days <- run.fun(model = full.model,  cohort = scapis[last_EDATUM<Visit1-180 | is.na(last_EDATUM),, ])
  resscapis_sa3_180days$time_windown = "180 days"
  
  atbclasses <- atbclasses[!grepl("_1yr", atbclasses)]
  resscapis_sa3_1y <- run.fun(model = full.model,  cohort = scapis[N1yr == 0, ])
  resscapis_sa3_1y$time_windown = "1 year"
  
  resscapis_sa3 <- rbind(resscapis_sa3_none, resscapis_sa3_30days, resscapis_sa3_180days, resscapis_sa3_1y)
  resscapis_sa3$cohort <- "SCAPIS"
  
  fwrite(resscapis_sa3, file='results/res_scapis_alpha_sensitivityanalyses_timewindows.tsv', sep = "\t")
  
  # Sensitivity analysis stratified age and sex
  # Import data set
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))

  atbclasses = grep("Class_.*yr", colnames(scapis), value=T)
  
  full.model_sex = full.model[-which(full.model=="Sex")]
  resscapis_male <- run.fun(model = full.model_sex,  cohort = scapis[Sex == "male", ])
  resscapis_male$model = "male"
  
  resscapis_female <- run.fun(model = full.model_sex,  cohort = scapis[Sex == "female", ])
  resscapis_female$model = "female"
  
  resscapis_agebelow55 <- run.fun(model = full.model,  cohort = scapis[age <= 55, ])
  resscapis_agebelow55$model = "age <= 55 years"
  
  resscapis_ageabove55 <- run.fun(model = full.model,  cohort = scapis[age > 55, ])
  resscapis_ageabove55$model = "age > 55 years"
  
  
  resscapis_age_sex <- rbind(resscapis_male, resscapis_female, resscapis_agebelow55, resscapis_ageabove55)
  resscapis_age_sex$cohort = "SCAPIS"
  
  fwrite(resscapis_age_sex, file='results/res_scapis_alpha_sensitivityanalyses_age_sex.tsv', sep = "\t")

  