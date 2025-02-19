  # Project antibiotic use and the gut microbiota
  
  
  # This script perform a meta-analysis of the associations 
  # from sensitivity analysis 
  
  
  # Meta-analysis of All antibiotics 
  
  suppressMessages({library(data.table)
  library(metafor)
  library(dplyr)
  library(tidyr)})
  
  
  # Meta-analysis function
  meta_fun <- function(out,res){
    
    exp <- unique(res$exposure)
    
    t <- lapply(exp, function(x){
      temp.data <- res[exposure == x & outcome == out, ]
      
      fit <- rma(yi = beta, sei = SE, data = temp.data , method = "EE" )
      fixed.res <- data.frame(beta = fit$b, SE = fit$se, LCI = fit$ci.lb , HCI = fit$ci.ub, p.value = fit$pval)
      Q <- fit$QE
      Qpval <- fit$QEp
      I2 <- round(fit$I2,1)
      
      data.frame(exposure = x, outcome = out,
                 fixed.res, 
                 Q = Q, Qpval = Qpval, I2 = I2, cohort = "Meta-analysis")
    })
    do.call(rbind,t)
  }
  
    # Import results for all antibiotics  -----------------------------------------------------------
    # SCAPIS 
    scapissa1_class <-  fread('scapis_class_sa1.tsv')
    scapissa2_class <- fread('scapis_class_sa2.tsv')
    scapissa1_class[,cohort:="SCAPIS"]
    scapissa2_class[,cohort:="SCAPIS"]
    var <- c("exposure","outcome","model","cohort","beta","SE")
    scapissa1_class <- scapissa1_class[,var,with=F]
    scapissa2_class <- scapissa2_class[,var,with=F]
    
    # SIMPLER
    simplersa1_class <-  fread('simpler_class_sa1.tsv')
    simplersa2_class <- fread('simpler_class_sa2.tsv')
    simplersa1_class[,cohort:="SIMPLER"]
    simplersa2_class[,cohort:="SIMPLER"]
    var <- c("exposure","outcome","model","cohort","beta","SE")
    simplersa1_class <- simplersa1_class[,var,with=F]
    simplersa2_class <- simplersa2_class[,var,with=F]
    
   
  res_sa1_class <- rbind(scapissa1_class, simplersa1_class)
  res_sa2_class <- rbind(scapissa2_class, simplersa2_class)
  

  # Meta-analysis sa1 (Hosp infect) model results
  species_all_cohorts <- res_sa1_class[,.N,by=.(outcome,exposure)] %>% filter(N==2) %>% pull(outcome) %>% unique(.)
  metaclasssa1 <- do.call(rbind, lapply(species_all_cohorts, meta_fun,  res = res_sa1_class))
  metaclasssa1 <- metaclasssa1 %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  metaclasssa1$N_scapis  <- unique(scapissa1_class$N)
  metaclasssa1$N_simpler <- unique(simplersa1_class$N)
  metaclasssa1$model <- "class_sa1"
  
  # Meta-analysis sa2 (Hosp any cause) model results
  species_all_cohorts <- res_sa2_class[,.N,by=.(outcome,exposure)] %>% filter(N==2) %>% pull(outcome) %>% unique(.)
  metaclasssa2 <- do.call(rbind, lapply(species_all_cohorts, meta_fun,  res = res_sa2_class))
  metaclasssa2 <- metaclasssa2 %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  metaclasssa2$N_scapis  <- unique(scapissa2_class$N)
  metaclasssa2$N_simpler <- unique(simplersa2_class$N)
  metaclasssa2$model <- "class_sa2"

  
  # Save results 
  meta.res <- rbind(metaclasssa1, metaclasssa2)
  fwrite(meta.res,'meta_byclass_sa12.tsv')
  
  message("End antibiotics by class")
  
