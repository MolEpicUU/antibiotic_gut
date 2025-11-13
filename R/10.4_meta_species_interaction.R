  # Project antibiotic use and the gut microbiome
  
  
  # This script perform a meta-analysis of the interaction terms and LRT p-values  
  # antibiotic use and species abundance in three cohorts (SCAPIS, MOS and SIMPLER)
  
  
  # Meta-analysis of All antibiotics 
  
  suppressMessages({library(data.table)
  library(metafor)
  library(dplyr)
  library(tidyr)
  library(metap)
  library(BiocParallel)})
  
 
  setwd('nobackup/users/baldanzi/atb_gut/results')
  
  
  # Meta-analysis function
  meta_fun <- function(out,res,m){
    
    exp <- unique(res$exposure[res$model == m])
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
  
  metapvalue_fun <- function(out,res,m){
    
    
    t <- lapply(classes, function(cl){
      
      temp.data <- res[grepl(cl, exposure) & outcome == out & model == m, .(outcome, model, p_omnibus, N)]
      temp.data <- unique(temp.data)
      
      p_om <- ifelse(nrow(temp.data)>2, sumlog(temp.data$p_lrt)$p, temp.data$p_lrt) # Fisher's method meta-analysis p-value
      

      data.frame(class = cl, outcome = out, model = m, p_omnibus = p_om)
    })
    do.call(rbind,t)
  }
  
  
  # Import results associations  ------------------------------------------------
  var <- c("exposure","outcome","model","cohort","beta","SE", "N", "p_omnibus", "p_lrt")
  
  # SCAPIS
  scapis <- list.files(pattern = "scapis_species__interaction_", path = "interaction")
  scapis <- file.path("interaction", scapis)
  scapis  <-  rbindlist(lapply(scapis, fread)) 
  scapis[,cohort:="SCAPIS"]
  scapis  <- scapis[,var,with=F]
  
  # MOS
  mos <- list.files(pattern = "mos_species__interaction_", path = "interaction")
  mos <- file.path("interaction", mos)
  mos  <-  rbindlist(lapply(mos, fread), fill = T) 
  mos[,cohort:="MOS"]
  mos <- mos[,var,with=F]
  
  # SIMPLER
  path_simpler_res <- '/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/atbgut/results/simpler_interaction'
  simpler <- list.files(pattern = "simpler_species__interaction_", path = path_simpler_res)
  simpler <- file.path(path_simpler_res, simpler)
  simpler  <-  rbindlist(lapply(simpler, fread)) 
  simpler[,cohort:="SIMPLER"]
  simpler <- simpler[,var,with=F]
  
  res <- rbind(scapis, mos, simpler)
  
  atbclasses <- unique(res$exposure)

  
  # Meta-analysis by model 
  species_all_cohorts <- Reduce(intersect, 
                                list(scapis = scapis$outcome, 
                                     simpler = simpler$outcome, 
                                     mos = mos$outcome))
  
  res <- res[outcome %in% species_all_cohorts & !is.na(exposure), ]
  
  unique_models <-  unique(res$model)
  message(paste("Number of models =", length(unique_models)))
  
  
  meta_res <- lapply(unique_models, function(model_var){
    
    rbindlist(bplapply(species_all_cohorts, meta_fun,  res = res, m = model_var, BPPARAM = MulticoreParam(16)))

    }
  )

  
  meta_res <- rbindlist(meta_res, fill=T)
  
  # Meta-analysis LRT test
  classes <- c("Class_Peni_Ext",  "Class_Peni_BetaS", "Class_lincosamides", "Class_SMZTMP",  "Class_FQs" ,  "Class_Peni_Comb", 
               "Class_NIT", "Class_Peni_BetaR",  "Class_TCLs", "Class_macrolides",  "Class_cephalosporins")
  
  
  meta_p <- lapply(unique_models, function(model_var){
    
    rbindlist(bplapply(species_all_cohorts, metapvalue_fun,  res = res, m = model_var, BPPARAM = MulticoreParam(16)))
    
  }
  )
  
  meta_p <- rbindlist(meta_p, fill=T)
  meta_p[, q_omnibus := p.adjust(p_omnibus, method = "BH"), by = model ]
  
  
  meta_res[, class := gsub(".*\\:", "", gsub("_\\d.*", "", exposure))]
  meta_res <- merge(meta_res, meta_p, by = c("class", "outcome", "model"))
  
  # Save results 
  fwrite(meta_res,'meta_species__sa_interaction.tsv', sep = "\t")
  
  message("End antibiotics - species interaction EE")
