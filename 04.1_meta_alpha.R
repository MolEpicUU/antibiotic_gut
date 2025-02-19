# Project Antibiotic Use and the Gut Microbiota

# Load packages
  library(data.table)
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(metafor)})
  

  # Import alpha diversity results 
  scapis  <- fread('res_scapis_alpha.tsv')
  simpler <- fread('res_simpler_alpha.tsv')
  mos     <- fread('res_mos_alpha.tsv')
  
  
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
    
    do.call(rbind,t)
    }
  
  
  # Main model results - Meta-analysis
  res.alpha <- rbind(scapis, simpler , mos[,-"df"] )
  meta.main <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha, m = "basic.model"))
  meta.main <- meta.main %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  # Disease model results - Meta-analysis
  meta.dis <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun,res = res.alpha, m = "full.model"))
  meta.dis <- meta.dis %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  
  meta.alpha.res <- rbind(meta.main, meta.dis, res.alpha)
  rm(meta.main, meta.dis, res.alpha, mos, scapis, simpler)
  
  # Save results 
  fwrite(meta.alpha.res,'meta_alpha.tsv')
  
  
  # Sensitivity analyses 
    scapis_sa  <- fread('res_scapis_alpha_sensitivityanalyses.tsv')
    simpler_sa <- fread('res_simpler_alpha_sensitivityanalyses.tsv')
    
    # Meta-analysis 
    res.alpha <- rbind(scapis_sa, simpler_sa)
    meta.sa1 <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha, m = "full.model_hospitalizedinfect"))
    meta.sa1 <- meta.sa1 %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    
    meta.sa2 <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha, m = "full.model_hospitalizedgeneral"))
    meta.sa2 <- meta.sa2 %>% group_by(outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
    meta.alpha.res_sa <- rbind(meta.sa1, meta.sa2, res.alpha)  
    
    # Save results 
    fwrite(meta.alpha.res_sa,'meta_alpha_sa.tsv')
    
    
    
  message("End meta-analysis")
  