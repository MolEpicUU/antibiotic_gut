  # Project antibiotic use and the gut microbiota
  
  
  # This script perform a meta-analysis of the association between
  # antibiotic use and species abundance in three cohorts (SCAPIS, MOS and SIMPLER)
  
  
  # Meta-analysis of All antibiotics 
  
  suppressMessages({library(data.table)
  library(metafor)
  library(dplyr)
  library(tidyr)
  library(BiocParallel)})
  
 
  setwd('results')
  
  
  # Meta-analysis function
  meta_fun <- function(out,res,m){
    
    exp <- unique(res$exposure)
    
    t <- lapply(exp, function(x){
      temp.data <- res[exposure == x & outcome == out & model == m, ]
      
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
  
  
  # Import results associations by antibiotic class ------------------------------------------------
  scapis  <-  fread('scapis_species_class.tsv') 
  scapis[,cohort:="SCAPIS"]
  var <- c("exposure","outcome","model","cohort","beta","SE", "N")
  scapis  <- scapis[,var,with=F]
  
  # MOS
  mos <- fread('mos_species_class.tsv')
  mos[,cohort:="MOS"]
  mos <- mos[,var,with=F]
  
  # SIMPLER
  simpler <- fread('simpler_species_class.tsv')
  simpler[,cohort:="SIMPLER"]
  simpler <- simpler[,var,with=F]
  
  res <- rbind(scapis, mos, simpler)
  res_basic <- res[model == "basic.model"]
  res_full <- res[model == "full.model"]
  
  
  # Meta-analysis basic model results
  species_all_cohorts <- res_basic[,.N,by=.(outcome,exposure)] %>% filter(N==3) %>% pull(outcome) %>% unique(.)
  meta.basic <- do.call(rbind, bplapply(species_all_cohorts, meta_fun,  res = res_basic, m = "basic.model", BPPARAM = MulticoreParam(16)))
  meta.basic <- meta.basic %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  meta.basic$N_scapis  <- unique(scapis[model == "basic.model", N])
  meta.basic$N_simpler <- unique(simpler[model == "basic.model", N])
  meta.basic$N_mos <- unique(mos[model == "basic.model", N])
  
  # Meta-analysis full model results
  species_all_cohorts <- res_full[,.N,by=.(outcome,exposure)] %>% filter(N==3) %>% pull(outcome) %>% unique(.)
  meta.full <- do.call(rbind, bplapply(species_all_cohorts, meta_fun, res = res_full, m = "full.model", BPPARAM = MulticoreParam(16)))
  meta.full <- meta.full %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  meta.full$N_scapis   <- unique(scapis[model == "full.model", N])
  meta.full$N_simpler  <- unique(simpler[model == "full.model", N])
  meta.full$N_mos <- unique(mos[model == "full.model", N])

  
  # Save results 
  meta.res <- rbind(meta.basic, meta.full)
  fwrite(meta.res,'meta_species_class.tsv', sep = "\t")
  
  message("End meta-analysis antibiotics by class")
  
  
  
