# Project Antibiotic Use and the Gut Microbiota

rm(list=ls())

# Load packages
  library(data.table)
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(metap)
    library(metafor)})

  setwd('nobackup/users/baldanzi/atb_gut/results')

  # Import alpha diversity results 
  scapis  <- fread('scapis_alpha__interaction.tsv')
  simpler <- fread('simpler_alpha__interaction.tsv')
  mos     <- fread('mos_alpha__interaction.tsv')
  mos[, cohort := "MOS"]
  
  
  
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
  
  
  metapvalue_fun <- function(out,res,m){
    
    
    t <- lapply(classes, function(cl){
      
      temp.data <- res[grepl(cl, exposure) & outcome == out & model == m, .(outcome, model, p_lrt, N)]
      temp.data <- unique(temp.data)
      
      p_lrt <- ifelse(nrow(temp.data)>2, sumlog(temp.data$p_lrt)$p, temp.data$p_lrt)
      
      
      data.frame(class = cl, outcome = out, model = m, p_omnibus = p_lrt)
    })
    do.call(rbind,t)
  }
  
  # ===========================================================================
  classes <- c("Class_Peni_Ext",  "Class_Peni_BetaS", "Class_lincosamides", "Class_SMZTMP",  "Class_FQs" ,  "Class_Peni_Comb", 
               "Class_NIT", "Class_Peni_BetaR",  "Class_TCLs", "Class_macrolides",  "Class_cephalosporins")
  
  
  # Sex
  res.alpha <- rbindlist(list(scapis, simpler , mos[,-"df"] ), fill=T)
  meta_sex <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun, res = res.alpha[model == "sex"], m = "sex"))
  meta_sex_p <- do.call(rbind, lapply(unique(res.alpha$outcome), metapvalue_fun, res = res.alpha[model == "sex"], m = "sex"))
  setDT(meta_sex_p)
  meta_sex_p[, q_omnibus := p.adjust(p_omnibus, "BH"), by = outcome]
  
  setDT(meta_sex)
  meta_sex[, class := gsub(".*\\:", "", gsub("_\\d.*", "", exposure))]
  meta_sex <- merge(meta_sex, meta_sex_p, by = c("class", "outcome", "model"))
  
  # Age
  meta_age <- do.call(rbind, lapply(unique(res.alpha$outcome), meta_alpha_fun,res = res.alpha[model == "age"], m = "age"))
  meta_age_p <- do.call(rbind, lapply(unique(res.alpha$outcome), metapvalue_fun,res = res.alpha[model == "age"], m = "age"))
  setDT(meta_age_p)
  meta_age_p[, q_omnibus := p.adjust(p_omnibus, "BH"), by = outcome]
  
  setDT(meta_age)
  meta_age[, class := gsub(".*\\:", "", gsub("_\\d.*", "", exposure))]
  meta_age <- merge(meta_age, meta_age_p, by = c("class", "outcome", "model"))
  
  meta_res <- rbind(meta_sex, meta_age)
  
  # Save results 
  fwrite(meta_res,'meta_alpha__interaction.tsv', sep = "\t")
  
  