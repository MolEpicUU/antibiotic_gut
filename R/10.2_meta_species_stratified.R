  # Project antibiotic use and the gut microbiota
  
  # This script perform a meta-analysis of the association between
  # antibiotic use and species abundance in three cohorts (SCAPIS, MOS and SIMPLER)
  # stratified by age and sex
  
  # Meta-analysis of All antibiotics 
  
  library(data.table)
  library(metafor)
  library(dplyr)
  library(tidyr)
  library(BiocParallel)
  
 
  setwd('nobackup/users/baldanzi/atb_gut/results')
  
  
  # Meta-analysis function
  meta_fun <- function(out,res,m){
    
    exp <- unique(res$exposure)
    message(paste("Length exp ==", length(exp)))
    
    print(m)
    
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
  var <- c("exposure","outcome","model","cohort","beta","SE", "N")
  
  # SCAPIS
  scapis  <-  fread('scapis_species_sa_stratified.tsv') 
  scapis[,cohort:="SCAPIS"]
  scapis  <- scapis[,var,with=F]
  
  # MOS
  mos <- fread('mos_species_class_sa_stratified.tsv', na.strings=c("NA", "", NA))
  mos[,cohort:="MOS"]
  mos <- mos[,var,with=F]
  
  # SIMPLER
  simpler <- fread('simpler_species_sa_stratified.tsv', na.strings=c("NA","",NA))
  simpler[,cohort:="SIMPLER"]
  simpler <- simpler[,var,with=F]
  
  res <- rbind(scapis, mos, simpler)
  
  atbclasses <- unique(res$exposure)
  atbclasses <- atbclasses[grep("FQs|lincosam|BetaR|BetaS|Peni_Ext|TCL", atbclasses)]
  
  res <- res[exposure %in% atbclasses, ]

  
  # Meta-analysis by model 
  species_all_cohorts <- res[model == "female",.N, by=.(outcome,exposure)] %>% filter(N==3) %>% pull(outcome) %>% unique(.)
  
  res <- res[outcome %in% species_all_cohorts, ]
  
  unique_models <-  unique(res$model)
  message(paste("Number of models =", length(unique_models)))
  
  meta_res <- lapply(unique_models, function(model_var){
    
    meta_res_temp <- do.call(rbind, bplapply(species_all_cohorts, meta_fun,  res = res, m = model_var, BPPARAM = MulticoreParam(16)))
    meta_res_temp <- meta_res_temp %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    
    meta_res_temp$N_scapis  <- unique(scapis[model == model_var, N])
    
    if(model_var == "age <= 55 years") {
      meta_res_temp$N_simpler <- 0
    } else {
      meta_res_temp$N_simpler <- unique(simpler[model == model_var, N])
    }
      
    if(model_var == "age > 55 years") {
      meta_res_temp$N_mos <- 0
    } else {
      meta_res_temp$N_mos <- unique(mos[model == model_var & !is.na(N), N])
    }
    
      
      return(meta_res_temp)
      
  })
  
  meta_res <- rbindlist(meta_res, fill=T)
  

  
  # Save results 
  fwrite(meta_res,'meta_species__sa_stratified.tsv', sep = "\t")
  
  
  # ------------------------------------------------------------------------------
  
  res <- rbind(scapis, mos)
  res <- res[exposure %in% atbclasses, ]
  
  # Meta-analysis by model 
  
  unique_models <-  c("age > 55 years", "age <= 55 years")
  message(paste("Number of models =", length(unique_models)))
  
  meta_res <- lapply(unique_models, function(model_var){
    
    meta_res_temp <- do.call(rbind, bplapply(species_all_cohorts, meta_fun,  res = res, m = model_var, BPPARAM = MulticoreParam(16)))
    meta_res_temp <- meta_res_temp %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    
    meta_res_temp$N_scapis  <- unique(scapis[model == model_var, N])
    
    if(model_var == "age <= 55 years") {
      meta_res_temp$N_simpler <- 0
    } else {
      meta_res_temp$N_simpler <- unique(simpler[model == model_var, N])
    }
    
    if(model_var == "age > 55 years") {
      meta_res_temp$N_mos <- 0
    } else {
      meta_res_temp$N_mos <- unique(mos[model == model_var & !is.na(N), N])
    }
    
    
    return(meta_res_temp)
    
  })
  
  meta_res <- rbindlist(meta_res, fill=T)
  
  
  
  # Save results 
  fwrite(meta_res,'meta_species__sa_stratified_withoutsimpler.tsv', sep = "\t")
  
  message("End antibiotics - species stratified")
  
  
