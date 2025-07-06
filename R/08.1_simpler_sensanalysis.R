# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  library(BiocParallel)
  
  # SIMPLER

  # This script will investigate the associations of number of previous antibiotics with species abundance
  # after removing participants that have been hospitalized with a diagnosis treated with antibiotics 
  

  setwd('users/baldanzi/')
  
  # Import species data
  clr.species <- fread('work/simpler_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  
  # Import model
  load('work/simpler_model_revision.Rdata')
  

# Create the linear regression function 

  linear.fun <- function(species, exposure, model = full.model){
    
    temp.data <- merge(simpler[, c("SIMPKEY",exposure, model) , with=F], clr.species, by = "SIMPKEY")
    
    if(!"data.table" %in% class(temp.data)) setDF(temp.data)
    
    
    
    res <- bplapply(species, function(y){
      
      form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
      
      fit <- lm(form,temp.data)
      ci <- data.frame(confint(fit, parm = exposure))
      colnames(ci) <- c("LCI","HCI")
      
      temp.res <- summary(fit)
      
      N = length(temp.res$residuals)
      res <- temp.res$coef[exposure, ]
      colnames(res) <- c("beta", "SE", "t.value","p.value")
      
      df <- data.frame(outcome = y, exposure = rownames(res), res, N = N)
      cbind(df,ci)
      
    }, BPPARAM = MulticoreParam(16))
    
    res <- do.call(rbind, res)
    
    setDT(res)
    
    res[,q.value := p.adjust(p.value, method = "BH")]
    res[,model := "full.model"]
    
    
  }
  
  ###########################################################
  #### Exclude those with hospitalization with infection ####
  ###########################################################
  
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  simpler <- simpler[hospinfect== "no",]
 
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value  = T )

  
  res <- linear.fun(prevalent.species,exposure = atbclasses, model = full.model)
  
  fwrite(res, file='results/simpler_class_sa1.tsv')
  
  
  message("End  Hosp Infect Model")
  
  #####################################################
  #### Exclude all hospitalizations with infection ####
  #####################################################
  
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, "")) 
  simpler <- simpler[hospgeneral== "no",]
  

  res <- linear.fun(prevalent.species, exposure = atbclasses, model = full.model)
  
  fwrite(res, file='results/simpler_class_sa2.tsv')
  
  
  message("End All Hospitalization Model")
