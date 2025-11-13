  # Project antibiotic use and the gut microbiota
  
  # This script perform a meta-analysis of the association between
  # antibiotic use and species abundance in three cohorts (SCAPIS, MOS and SIMPLER)
  # stratified by age and sex
  
  # Meta-analysis of All antibiotics 
  
  suppressMessages({library(data.table)
  library(metafor)
  library(dplyr)
  library(tidyr)
  library(BiocParallel)})
  
  
  stratum <- commandArgs(trailingOnly = TRUE)
  #  stratum <- "male"

  print(stratum)
  
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
  scapis  <-  fread(paste0('scapis_species_class_sa_stratified_allantibiotics__', stratum, '.tsv')) 
  scapis[,cohort:="SCAPIS"]
  scapis  <- scapis[,var,with=F]
  
  # MOS
  mos <- fread(paste0('mos_species_class_sa_stratified_allantibiotics__', stratum, '.tsv')) 
  mos[,cohort:="MOS"]
  mos <- mos[,var,with=F]
  
  # SIMPLER
  if(stratum!="agebelow55"){
  simpler <- fread(paste0('simpler_species_class_sa_stratified_allantibiotics__', stratum, '.tsv')) 
  simpler[,cohort:="SIMPLER"]
  simpler <- simpler[,var,with=F]
  
    res <- rbind(scapis, mos, simpler)
    
    species_all_cohorts <- Reduce(intersect, list(scapis$outcome, mos$outcome, simpler$outcome))
    
  } else{
    
    res <- rbind(scapis, mos)
    
    simpler <- fread('simpler_species_class_sa_stratified_allantibiotics__ageabove55.tsv')
    species_all_cohorts <- Reduce(intersect, list(scapis$outcome, mos$outcome, simpler$outcome))
    
  }
  
  atbclasses <- unique(res$exposure)

  
  # Meta-analysis by model 
  
  
  res <- res[outcome %in% species_all_cohorts, ]
  
  unique_models <-  unique(res$model)
  message(paste("Number of models =", length(unique_models)))
  

    meta_res <- do.call(rbind, bplapply(species_all_cohorts, meta_fun,  res = res, m = stratum, BPPARAM = MulticoreParam(16)))
    meta_res <- meta_res %>% mutate(q.value = p.adjust(p.value, method = "BH"))
    setDT(meta_res)
    
    meta_res  <- merge(meta_res, scapis[, .(exposure, outcome, N_scapis=N)], by= c("exposure", "outcome"), all.x=T)
    meta_res  <- merge(meta_res, mos[, .(exposure, outcome, N_mos=N)], by= c("exposure", "outcome"), all.x=T)

    
    if(stratum == "agebelow55") {
      meta_res$N_simpler <- 0
    } else {
      meta_res  <- merge(meta_res, simpler[, .(exposure, outcome, N_simpler=N)], by= c("exposure", "outcome"), all.x=T)
    }
      
  
  # Save results 
  fwrite(meta_res, paste0('meta_species__sa_stratified_allantibiotics_', stratum ,'.tsv'), sep = "\t")

  
  head(meta_res)
  
message(stratum)
message("END")
  
  
