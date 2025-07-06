# Project Antibiotic Use and the Gut Microbiota


rm(list=ls())

# Load packages
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(metafor)

  setwd("revision_NatMed/revision__results")

  # Import alpha diversity results 
  scapis  <- fread('res_scapis_alpha_negexposure.tsv')
  simpler <- fread('res_simpler_alpha_negexposure.tsv')
  mos     <- fread('res_mos_alpha_negexposure.tsv')
  scapis_preatbzero   <- fread('res_scapis_alpha_negexposure_preatbzero.tsv')
  simpler_preatbzero  <- fread('res_simpler_alpha_negexposure_preatbzero.tsv')
  mos_preatbzero      <- fread('res_mos_alpha_negexposure_preatbzero.tsv')
  
  
  # Meta-analysis function ---------------------------------------------------------------------------------------
  
  meta_alpha_fun <- function(out,res,m){
    exp <- unique(res$exposure)
    t <- lapply(exp, function(x){
      temp.data <- res.alpha[exposure == x & outcome == out & model == m, ]
      
      fit <- rma(yi = beta, sei = SE, data = temp.data , method = "EE" )
      fixed.res <- data.frame(beta = fit$b, SE = fit$se, LCI = fit$ci.lb , HCI = fit$ci.ub, p.value = fit$pval)
      Q <- fit$QE
      Qpval <- fit$QEp
      I2 <- round(fit$I2,1)
      
      data.frame(exposure = x, outcome = out, model = m,
                 fixed.res, 
                 Q = Q, Qpval = Qpval, I2 = I2, cohort = "Meta-analysis")
    })
    do.call(rbind,t)}
  
  
  # Basic model results - Meta-analysis
  res.alpha <- rbind(scapis, simpler , mos[,-"df"] )
  #meta_basic <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha, m = "basic.model"))
  #meta_basic <- meta_basic %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  # Full model results - Meta-analysis
  meta_full <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun,res = res.alpha, m = "full.model"))
  meta_full <- meta_full  %>% filter(exposure %in% grep("after", meta_full$exposure, value=T)) %>% 
                                     group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  
  meta.alpha.res <- rbind(meta_full, res.alpha[grep("after", exposure),])
  
  # Save results 
  fwrite(meta.alpha.res,'meta_alpha_negexposure.tsv')
  
  message("End meta-analysis 1")
  
  # Pre atb == zero ---------------------------------------------------------------
  
  res.alpha <- rbind(scapis_preatbzero, simpler_preatbzero , mos_preatbzero[,-"df"] )
  meta_full_preatbzero <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun,res = res.alpha, m = "full.model"))
  meta_full_preatbzero <- meta_full_preatbzero  %>% filter(exposure %in% grep("after", meta_full_preatbzero$exposure, value=T)) %>% 
    group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  
  meta.alpha.res_preatbzero <- rbind(meta_full_preatbzero, res.alpha[grep("after", exposure),])
  
  # Save results 
  fwrite(meta.alpha.res_preatbzero,'meta_alpha_negexposure_preatbzero.tsv')
  