# Project Antibiotic Use and the Gut Microbiota

rm(list=ls())

# Load packages
  library(data.table)
  library(lmerTest)
    library(dplyr)
    library(tidyr)
    library(car)

  setwd('nobackup/users/baldanzi/atb_gut/')

# This script will investigate the associations of number of previous antibiotics with alpha diversity


  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  load('work/mos_model_revision.Rdata')
  
  # Association between number of antibiotics and alpha diversity  
  

  # MOS ------------------------------------------------------------------------
  
  mosatbclasses <-  grep("Class_.*yr", colnames(mos), value=T)

  mix.fun <- function(y, exp = c("N1yr","N1_4yr","N4_8yr"), model = basic.model, data){
    
    
    tryCatch({
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
      
      data.frame(outcome = y, exposure = rownames(res), res, LCI=lci,HCI=hci, N = N, cohort = "MOS", gvif = gvif, message=NA)
    }, error = function(e){
      
      data.frame(outcome = y, exposure = NA, beta=NA,SE=NA,df=NA,t.value=NA,p.value=NA, LCI=NA,HCI=NA, N = NA, cohort = "MOS", gvif = NA, message=e$message)
      
    })
    
  }
  
  
  
  run.fun <- function(model = basic.model, cohort = mos){
    
    resmosclass_shannon <- mix.fun(y="shannon", exp = mosatbclasses, data = cohort, model = model)
    resmosclass_rich <- mix.fun(y="richness", exp = mosatbclasses, data = cohort, model= model )
    resmosclass_invsimpson <- mix.fun(y="invsimpson", exp = mosatbclasses, data = cohort, model= model )
    
    res <- rbindlist(list(resmosclass_shannon, resmosclass_rich, resmosclass_invsimpson), fill=T)
    
    res <- res %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    return(res)
  }
  
  # Run analyses
  resmos <- run.fun(model = basic.model)
  resmos$model <- "basic.model"
  resmos_full <- run.fun(model = full.model)
  resmos_full$model <- "full.model"
  

  res.alpha = rbind(resmos, resmos_full)
  
  # Save results 
  fwrite(res.alpha, file='results/res_mos_alpha.tsv')
  
  message("End linear regression")
  
  
  
  # Sensitivity analysis using different time windows of recent antibiotic use ####
  # Import data set
  mos <- fread("work/mos_working_dataset_revision_noatbexclusion.tsv", na.strings = c("NA", NA, ""))
  
  mosatbclasses = grep("Class_.*yr", colnames(mos), value=T)
  
  
  resmos_sa3_none <- run.fun(model = full.model,  cohort = mos)
  resmos_sa3_none$time_windown = "none"
  
  resmos_sa3_30days <- run.fun(model = full.model,  cohort = mos[last_EDATUM<Visit1-30 | is.na(last_EDATUM), ])
  resmos_sa3_30days$time_windown = "30 days"
  
  mosatbclasses = grep("Class_.*yr", colnames(mos), value=T)
  mosatbclasses = mosatbclasses[-which(mosatbclasses == "Class_Peni_Comb_1yr")] # zero users 
  resmos_sa3_180days <- run.fun(model = full.model,  cohort = mos[last_EDATUM<Visit1-180 | is.na(last_EDATUM), ])
  resmos_sa3_180days$time_windown = "180 days"
  
  mosatbclasses <- mosatbclasses[!grepl("_1yr", mosatbclasses)]
  resmos_sa3_1y <- run.fun(model = full.model,  cohort = mos[N1yr == 0, ])
  resmos_sa3_1y$time_windown = "1 year"
  
  resmos_sa3 <- rbind(resmos_sa3_none, resmos_sa3_30days, resmos_sa3_180days, resmos_sa3_1y)
  resmos_sa3$cohort <- "mos"
  
  fwrite(resmos_sa3, file='results/res_mos_alpha_sensitivityanalyses_timewindows.tsv', sep = "\t")
  
  
  message("End sensitivity analysis of time windonws \n")
  
  
  # Sensitivity analysis stratified age and sex
  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  
  mosatbclasses <- grep("Class_.*yr", colnames(mos), value=T)
  mosatbclasses <- grep("NIT", mosatbclasses , value=T, invert=T) 
  mosatbclasses <- grep("Peni_Comb", mosatbclasses , value=T, invert=T) 
  mosatbclasses <- grep("SMZTMP", mosatbclasses , value=T, invert=T) 
  
  full.model_sex = full.model[-which(full.model=="sex")]
  resmos_male <- run.fun(model = full.model_sex,  cohort = mos[sex == 1, ])
  resmos_male$model = "male"
  
  mosatbclasses = grep("Class_.*yr", colnames(mos), value=T)
  resmos_female <- run.fun(model = full.model_sex,  cohort = mos[sex == 2, ])
  resmos_female$model = "female"
  
  message("\nSex stratified")
  
  mosatbclasses = grep("Class_.*yr", colnames(mos), value=T)
  mosatbclasses <- grep("Peni_Comb", mosatbclasses , value=T, invert=T) 
  
  resmos_agebelow55 <- run.fun(model = full.model,  cohort = mos[age <= 55, ])
  resmos_agebelow55$model = "age <= 55 years"
  
  
  mosatbclasses = grep("Class_.*yr", colnames(mos), value=T)
  resmos_ageabove55 <- run.fun(model = full.model,  cohort = mos[age > 55, ])
  resmos_ageabove55$model = "age > 55 years"
  
  
  resmos_age_sex <- rbindlist(list(resmos_male, resmos_female, resmos_agebelow55, resmos_ageabove55), fill=T)
  resmos_age_sex$cohort = "MOS"
  
  fwrite(resmos_age_sex, file='results/res_mos_alpha_sensitivityanalyses_age_sex.tsv', sep = "\t")
  
  message("End")
  
  