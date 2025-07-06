  # Project antibiotic use and the gut microbiota
  
  
  # This script perform a meta-analysis of the association between
  # antibiotic use and species abundance in three cohorts (SCAPIS, MOS and SIMPLER)
  # after exclusion of the most influential observation in each cohort 

  # Only the species with high heterogeneity in the full model meta-analysis were included
  
  
  # Meta-analysis
  
  library(data.table)
  library(metafor)
  library(dplyr)
  library(tidyr)
  
  rm(list=ls())
 
  setwd('results')
  
  
  # Meta-analysis function
  meta_fun <- function(out,res,m){
    
    exp <- unique(res$exposure)
    
    t <- lapply(exp, function(x){
      temp.data <- res[exposure == x & outcome == out & model == m, ]
      
      fit <- rma(yi = beta_infl, sei = se_infl, data = temp.data , method = "EE" )
      fixed.res <- data.frame(beta = fit$b, SE = fit$se, LCI = fit$ci.lb , HCI = fit$ci.ub, p.value = fit$pval)
      Q <- fit$QE
      Qpval <- fit$QEp
      I2 <- round(fit$I2,1)
      
      data.frame(exposure = x, outcome = out, model = m,
                 fixed.res, 
                 Q = Q, Qpval = Qpval, I2 = I2, cohort = "Meta-analysis")
    })
    do.call(rbind,t)
  }
  
    # Import results for all antibiotics  -----------------------------------------------------------
    # SCAPIS 
    scapis <-  fread('scapis_infobs_byclass.tsv')
    scapis[,cohort:="SCAPIS"]
    var <- c("exposure","outcome","model","cohort","beta_infl","se_infl")
    scapis <- scapis[,var,with=F]

    
    # MOS
    mos <-  fread('mos_infobs_byclass.tsv')
    mos[,cohort:="SCAPIS"]
    var <- c("exposure","outcome","model","cohort","beta_infl","se_infl")
    mos <- mos[,var,with=F]
    
    # SIMPLER
    simpler <-  fread('simpler_infobs_byclass.tsv')
    simpler[,cohort:="SCAPIS"]
    var <- c("exposure","outcome","model","cohort","beta_infl","se_infl")
    simpler <- simpler[,var,with=F]
    
   
  res_main <- rbind(scapis, simpler, mos)
  
  
  # Meta-analysis full model results
  species_all_cohorts <- unique(res_main$outcome)
  metainflobs_byclass <- do.call(rbind, lapply(species_all_cohorts, meta_fun, res = res_main, m = "full.model"))
  
  fullmeta <- fread('meta_species_class.tsv')
  fullmetaexpout <- fullmeta[q.value<.05 & Qpval<.05 & model == "full.model",.(exposure, outcome)]
  
  stopifnot("Species missing from the main analysis" = all(fullmetaexpout$outcome %in% metainflobs_byclass$outcome))
  metainflobs_byclass <- merge(metainflobs_byclass, fullmetaexpout, by=c("exposure","outcome"))
  
  fwrite(metainflobs_byclass, "metainflobs_byclass.tsv", sep="\t")
  
  message("End influential observation for total atb")
  
  setDT(metainflobs_byclass)
  expoutinflobs <- metainflobs_byclass[p.value>0.05,.(exposure, outcome)]
  
  expout_tokeep <- fsetdiff(fullmeta[model=="full.model", .(exposure,outcome)], expoutinflobs)
  
  fullmetaalltab_clean <- merge(fullmeta, expout_tokeep, by=c("exposure","outcome"))
  
  
  temp <- merge(metainflobs_byclass[p.value<0.05,.(exposure, outcome, beta_io = beta)], fullmeta[, .(exposure, outcome, beta)], by = c("exposure","outcome"))
  stopifnot(nrow(temp[beta_io>0  & beta<0,])==0)
  stopifnot(nrow(temp[beta_io<0  & beta>0,])==0)
  
  
  fwrite(fullmetaalltab_clean, "meta_species_class_clean.tsv", sep = "\t")
  
  message("Meta all antibiotics clean saved")
  



