  # Project antibiotic use and the gut microbiota
  
  # This script perform a meta-analysis of the association between
  # antibiotic use and species abundance in SCAPIS, MOS and SIMPLER
  
  
  # Meta-analysis of All antibiotics 
  
  suppressMessages({library(data.table)
  library(metafor)
  library(dplyr)
  library(tidyr)
  library(BiocParallel)})
  
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
  scapisbasic  <-  fread('scapis_species_linear_class_N1N2N3.tsv') 
  scapisfull   <- fread('scapis_species_linear_class_N1N2N3_Dis.tsv')
  scapisbasic[,cohort:="SCAPIS"]
  scapisfull[ ,cohort:="SCAPIS"]
  var <- c("exposure","outcome","model","cohort","beta","SE")
  scapisbasic  <- scapisbasic[,var,with=F]
  scapisfull   <- scapisfull[,var,with=F]
  
  # MOS
  mosbasic <- fread('mos_species_linear_class_N1N2N3.tsv')
  mosfull  <- fread('mos_species_linear_class_N1N2N3_Dis.tsv')
  mosbasic[,cohort:="MOS"]
  mosfull[, cohort:="MOS"]
  mosbasic <- mosbasic[,var,with=F]
  mosfull  <- mosfull[,var,with=F]
  
  # SIMPLER
  simplerbasic <- fread('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/atbgut/results/simpler_species_linear_class_N1N2N3.tsv')
  simplerfull  <- fread('/proj/nobackup/sens2019512/wharf/baldanzi/baldanzi-sens2019512/atbgut/results/simpler_species_linear_class_N1N2N3_Dis.tsv')
  simplerbasic[,cohort:="SIMPLER"]
  simplerfull[, cohort:="SIMPLER"]
  simplerbasic <- simplerbasic[,var,with=F]
  simplerfull  <- simplerfull[,var,with=F]
  
  res_basic <- rbind(scapisbasic, mosbasic, simplerbasic)
  res_full <- rbind(scapisfull, mosfull, simplerfull)
  
  
  # Meta-analysis basic model results
  species_all_cohorts <- res_basic[,.N,by=.(outcome,exposure)] %>% filter(N==3) %>% pull(outcome) %>% unique(.)
  meta.basic <- do.call(rbind, bplapply(species_all_cohorts, meta_fun,  res = res_basic, m = "basic.model", BPPARAM = MulticoreParam(16)))
  meta.basic <- meta.basic %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  meta.basic$N_scapis  <- unique(scapisbasic$N)
  meta.basic$N_simpler <- unique(simplerbasic$N)
  meta.basic$N_mos <- unique(mosbasic$N)
  
  # Meta-analysis Disease model results
  species_all_cohorts <- res_full[,.N,by=.(outcome,exposure)] %>% filter(N==3) %>% pull(outcome) %>% unique(.)
  meta.Dis <- do.call(rbind, bplapply(species_all_cohorts, meta_fun, res = res_full, m = "full.model", BPPARAM = MulticoreParam(16)))
  meta.Dis <- meta.Dis %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  meta.Dis$N_scapis   <- unique(scapisbasic$N)
  meta.Dis$N_simpler  <- unique(simplerbasic$N)
  meta.Dis$N_mos <- unique(mosbasic$N)

  
  # Save results 
  meta.res <- rbind(meta.basic, meta.Dis)
  fwrite(meta.res,'meta_byclass.tsv')

  message("End antibiotics by class")
  
