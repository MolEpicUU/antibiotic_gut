  # Project antibiotic use and the gut microbiota
  
  
  # This script perform a meta-analysis of the association between
  # a single antibiotic dose and species in three cohorts (SCAPIS, MOS and SIMPLER)
  
  
  # Meta-analysis
  
  library(data.table)
    library(metafor)
    library(dplyr)
    library(tidyr)
  
 
  setwd('results')
  
  
  # Meta-analysis function
  meta_fun <- function(out,res,m){
    
    exp <- unique(res$exposure[res$model==m])
    
    t <- lapply(exp, function(x){
      temp.data <- res[exposure == x & outcome == out & model == m, ]
      
      tryCatch(
        expr = {
          fit <- rma(yi = beta, sei = SE, data = temp.data , method = "EE" )
          fixed.res <- data.frame(beta = fit$b, SE = fit$se, LCI = fit$ci.lb , HCI = fit$ci.ub, p.value = fit$pval)
          Q <- fit$QE
          Qpval <- fit$QEp
          I2 <- round(fit$I2,1)
          
          data.frame(exposure = x, outcome = out, model = m, fixed.res,  Q = Q, Qpval = Qpval, I2 = I2, cohort = "Meta-analysis", war = NA)
        }, 
        warning = function(w) { fit <- rma(yi = beta, sei = SE, data = temp.data , method = "EE" )
        fixed.res <- data.frame(beta = fit$b, SE = fit$se, LCI = fit$ci.lb , HCI = fit$ci.ub, p.value = fit$pval)
        Q <- fit$QE
        Qpval <- fit$QEp
        I2 <- round(fit$I2,1)
        data.frame(exposure = x, outcome = out, model = m, fixed.res,  Q = Q, Qpval = Qpval, I2 = I2, cohort = "Meta-analysis", war = conditionMessage(w))
        }, 
        error = function(w) { 
          fixed.res <- data.frame(beta = NA, SE = NA, LCI = NA , HCI = NA, p.value = NA)
        data.frame(exposure = x, outcome = out, model = m, fixed.res,  Q = NA, Qpval = NA, I2 = NA, cohort = "Meta-analysis", war = conditionMessage(w))
        }
        )
      
    })
    do.call(rbind,t)
  }
  
    # Import results singledose species -----------------------------------------------------------
    # SCAPIS 
    scapis <-  fread('scapis_species_singledose.tsv')
    var <- c("exposure","outcome","model","cohort","beta","SE")
    scapis  <- scapis[,var,with=F]

    
    # MOS
    mos <- fread('mos_species_singledose.tsv')
    mos <- mos[,var,with=F]

    
    # SIMPLER
    simpler <- fread('simpler_species_singledose.tsv')
    simpler <- simpler[,var,with=F]
    
  res_allcohorts  <- rbind(scapis, mos, simpler)
  species_all_cohorts <- res_allcohorts[model=="basic.model",.N,by=.(outcome,exposure)] %>% filter(N==3) %>% pull(outcome) %>% unique(.)
  
  
  # Meta-analysis basic model results ####
  
  meta.basic <- do.call(rbind, lapply(species_all_cohorts, meta_fun, res =  res_allcohorts , m = "basic.model"))
  meta.basic <- meta.basic %>%  mutate(q.value = p.adjust(p.value, method = "BH"))
  
  # Meta-analysis full model results ####
  meta.full <- do.call(rbind, lapply(species_all_cohorts, meta_fun, res = res_allcohorts, m = "full.model"))
  meta.full <- meta.full %>%  mutate(q.value = p.adjust(p.value, method = "BH"))

  # Save results
  meta.res <- rbind(meta.basic, meta.full)
  fwrite(meta.res,'meta_species_singledose.tsv', sep = '\t')
  
  
  # ALPHA --------------------------------------------------------------------------------------------------------
  
  # SCAPIS 
  scapis <-  fread('scapis_alpha_singledose.tsv')
  var <- c("exposure","outcome","model","cohort","beta","SE")
  scapis  <- scapis[,var,with=F]
  
  
  # MOS
  mos <- fread('mos_alpha_singledose.tsv')
  mos <- mos[,var,with=F]
  
  
  # SIMPLER
  simpler <- fread('simpler_alpha_singledose.tsv')
  simpler <- simpler[,var,with=F]
  
  res_allcohorts  <- rbind(scapis, mos, simpler)
  
  # Meta-analysis basic model results ####
  
  meta.basic <- do.call(rbind, lapply(c("shannon","richness","invsimpson"), meta_fun, res =  res_allcohorts , m = "basic.model"))
  meta.basic <- meta.basic %>%  mutate(q.value = p.adjust(p.value, method = "BH"))
  
  # Meta-analysis full model results ####
  meta.full <- do.call(rbind, lapply(c("shannon","richness","invsimpson"), meta_fun, res = res_allcohorts, m = "full.model"))
  meta.full <- meta.full %>%  mutate(q.value = p.adjust(p.value, method = "BH"))
  
  # Save results
  meta.res <- rbind(meta.basic, meta.full)
  fwrite(meta.res,'meta_alpha_singledose.tsv', sep = '\t')
  
  message("End alpha meta")
  
 
  
  