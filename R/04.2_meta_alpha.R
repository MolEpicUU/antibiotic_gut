# Project Antibiotic Use and the Gut Microbiota


rm(list=ls())

# Load packages
  library(data.table)
    library(dplyr)
    library(tidyr)
    library(metafor)

  setwd('revision_NatMed/revision__results')

  # Import alpha diversity results 
  scapis  <- fread('res_scapis_alpha.tsv')
  simpler <- fread('res_simpler_alpha.tsv')
  mos     <- fread('res_mos_alpha.tsv')
  
  
  
  # Meta-analysis function ---------------------------------------------------------------------------------------
  
  meta_alpha_fun <- function(out,res,m){
    
    exp <- unique(res$exposure)
    
    t <- lapply(exp, function(x){
      temp.data <- res.alpha[exposure == x & outcome == out & model == m, ]
      if(nrow(temp.data)==0) stop(paste("temp.data has 0 rows.", "exposure = ", x, "outcome = ", out, "model = ", m))
      
      fit <- rma(yi = beta, sei = SE, data = temp.data , method = "EE" )
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
  
  
  # Basic model results - Meta-analysis
  res.alpha <- rbindlist(list(scapis, simpler , mos[,-"df"] ), fill=T)
  meta.basic <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha, m = "basic.model"))
  meta.basic <- meta.basic %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  # Full model results - Meta-analysis
  meta.full <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun,res = res.alpha, m = "full.model"))
  meta.full <- meta.full %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  
  meta.alpha.res <- rbind(meta.basic, meta.full, res.alpha)
  rm(meta.basic, meta.full, res.alpha, mos, scapis, simpler)
  
  # Save results 
  fwrite(meta.alpha.res,'meta_alpha.tsv')
  
  
  # Hospitalization ------------------------------------------------------------
  
  message("Hospitalization")
  
  # Sensitivity analyses 
    scapis_sa  <- fread('res_scapis_alpha_sensitivityanalyses.tsv')
    simpler_sa <- fread('res_simpler_alpha_sensitivityanalyses.tsv')
    
    # Meta-analysis 
    res.alpha <- rbindlist(list(scapis_sa, simpler_sa), fill=T)
    
    meta.alpha.res_sa <- lapply(unique(res.alpha$model), function(mod) {
    
    meta.sa <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha, m = mod))
    meta.sa <- meta.sa %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    
    return(meta.sa)
    
    })
  
    meta.alpha.res_sa <- rbindlist(meta.alpha.res_sa, fill=T)  
    meta.alpha.res_sa <- rbindlist(list(meta.alpha.res_sa, res.alpha), fill=T)
    
    # Save results 
    fwrite(meta.alpha.res_sa,'meta_alpha_sa.tsv')
    
    
    # TIme windowns -------------------------------------------------------------
    
    message("Time-windowns")
    
    # Sensitivity analyses 
    scapis_sa  <- fread('res_scapis_alpha_sensitivityanalyses_timewindows.tsv', na.strings=c("NA", NA, ""))
    mos_sa <- fread('res_mos_alpha_sensitivityanalyses_timewindows.tsv', na.strings=c("NA", NA, ""))
    simpler_sa <- fread('res_simpler_alpha_sensitivityanalyses_timewindows.tsv', na.strings=c("NA", NA, ""))
    
    # Meta-analysis 
    res.alpha <- rbindlist(list(scapis_sa, simpler_sa, mos_sa), fill=T)
    setnames(res.alpha, "time_windown", "model")
    
    res.alpha <- res.alpha[!is.na(exposure), ]
    

      meta.alpha.res_sa <- lapply(c("none","30 days", "180 days"), function(mod) {
      
      meta.sa <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha[model == mod,], m = mod))
      meta.sa <- meta.sa %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
      
      return(meta.sa)
      
    })
      
      meta.alpha.res_sa_1year <- lapply(c("1 year"), function(mod) {
        
        meta.sa <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha[model == mod,], m = mod))
        meta.sa <- meta.sa %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
        
        return(meta.sa)
        
      })
    
    meta.alpha.res_sa <- rbindlist(c(meta.alpha.res_sa, meta.alpha.res_sa_1year), fill=T)  
    
    # Save results 
    fwrite(meta.alpha.res_sa,'meta_alpha_sa_timewindows.tsv')
    
    # Age and Sex stratified -------------------------------------------------------------
    
    message("Age and Sex stratified")
    
    # Sensitivity analyses 
    scapis_sa  <- fread('res_scapis_alpha_sensitivityanalyses_age_sex.tsv', na.strings=c("NA", NA, ""))
    mos_sa <- fread('res_mos_alpha_sensitivityanalyses_age_sex.tsv', na.strings=c("NA", NA, ""))
    mos_sa <- mos_sa[ model != "age > 55 years", ]
    simpler_sa <- fread('res_simpler_alpha_sensitivityanalyses_age_sex.tsv', na.strings=c("NA", NA, ""))
    
    # Meta-analysis 
    res.alpha <- rbindlist(list(scapis_sa, simpler_sa, mos_sa), fill=T)

    res.alpha <- res.alpha[!is.na(exposure), ]
    
    
    meta.alpha.res_sa <- lapply(unique(res.alpha$model), function(mod) {
      
      meta.sa <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha[model == mod,], m = mod))
      meta.sa <- meta.sa %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
      
      return(meta.sa)
      
    })
    
    
    meta.alpha.res_sa <- rbindlist(meta.alpha.res_sa, fill=T)  
    
    # Save results 
    fwrite(meta.alpha.res_sa,'meta_alpha_sa__age_sex.tsv')
    
    
  message("End meta-analysis")
  