# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  library(BiocParallel)
  
  # SCAPIS

  # This script will investigate the associations of number of previous antibiotics with species abundance
  # after removing participants that have been hospitalized with a diagnosis treated with antibiotics 
  

  setwd('nobackup/users/baldanzi/atb_gut/')
  
  # Import species data
  clr.species <- fread('work/scapis_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  
  # Import model
  load('work/scapis_model_revision.Rdata')
  

# Create the linear regression function 

  linear.fun <- function(species, exposure, model = full.model){
    
    temp.data <- merge(scapis[, c("Subject",exposure, model) , with=F], clr.species, by = "Subject")
    
    if(!"data.table" %in% class(temp.data)) setDF(temp.data)
    
    rownames(temp.data) <- temp.data$Subject
    
    
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
  
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  scapis <- scapis[hospinfect== "no",]
 
  
  atbclasses = grep("Class_.*yr", colnames(scapis), value  = T )

  
  res <- linear.fun(prevalent.species,exposure = atbclasses, model = full.model)
  
  fwrite(res, file='results/scapis_class_sa1.tsv')
  
  
  message("End  Hosp Infect Model")
  
  #####################################################
  #### Exclude all hospitalizations with infection ####
  #####################################################
  
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, "")) 
  scapis <- scapis[hospgeneral== "no",]
  
  res <- linear.fun(prevalent.species, exposure = atbclasses, model = full.model)
  
  fwrite(res, file='results/scapis_class_sa2.tsv')
  
  
  message("End All Hospitalization Model")
