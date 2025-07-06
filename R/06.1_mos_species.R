# Project Antibiotic Use and the Gut Microbiota

# Script created to use the MOS data 

# Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  library(lmerTest)
  
  t0 = Sys.time()

# This script will investigate the associations of number of previous antibiotics with species abundance
  
# MIXED MODEL REGRESSION ####

  setwd('nobackup/users/baldanzi/atb_gut/')

  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  clr.species <- fread('work/mos_clr_species.tsv', na.strings = c("NA", NA, ""),  data.table = F)

  # Import Basic model and species by abundance
  load('work/mos_model_revision.Rdata')


  # Create the mix-model regression function 

  linear.fun <- function(species, exposure, model = basic.model){
  
    temp.data <- merge(mos[, c("lopnrMOS",exposure, model,"family") , with=F], clr.species, by = "lopnrMOS")
  
    if(!"data.table" %in% class(temp.data)) setDF(temp.data)
  
    rownames(temp.data) <- temp.data$lopnrMOS
  
  
  res <- bplapply(species, function(y){
    
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"),"+(1|family)"))
          
          fit <- lmer(form,temp.data)
          
          gvif <- car::vif(fit)[exposure, "GVIF", drop = F]
    
          temp.res <- summary(fit)
    
          N = length(temp.res$residuals)
          res <- temp.res$coef[exposure, ]
          colnames(res) <- c("beta", "SE","df", "t.value","p.value")
          suppressMessages(ci <- confint(fit,parm=exposure))
          lci = data.frame(LCI = ci[,1])
          hci = data.frame(HCI = ci[,2])  
    
          df <- data.frame(outcome = y, exposure = rownames(res), res, lci=lci,hci=hci, N = N)
          cbind(df, gvif)
  
    }, BPPARAM = MulticoreParam(16))
  
  res <- do.call(rbind, res)
  
  setDT(res)
  
  res[,q.value := p.adjust(p.value, method = "BH")]
  res[,model := ""]
  
}

  # Basic model  ####----------------------------------------------
  
  message("Basic model by Class ")
  
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
 

  res <- linear.fun(prevalent.species, exposure = atbclasses )
  res$model <- "basic.model"

  fwrite(res, file='results/mos_species_class.tsv', sep = "\t")
  
  
  # Full model #### ----------------------------------------------------
  
  res_full <- linear.fun(prevalent.species, exposure = atbclasses,  model = full.model)
  res_full$model <-  "full.model"
  
  res <- rbind(res, res_full)
  
  fwrite(res, file='results/mos_species_class.tsv', sep = "\t")

  t1 = Sys.time()
  print(t1-t0)
  message("End")
  
