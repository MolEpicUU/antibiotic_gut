# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  
  t0 = Sys.time()

# This script will investigate the associations of number of previous antibiotics with species abundance
  
# LINEAR REGRESSION ####

  setwd('nobackup/users/baldanzi/atb_gut/')
  
  # Import data set
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))

  clr.species <- fread('work/scapis_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  

  # Import Basic model
  load('work/scapis_model_revision.Rdata')


# Create the linear regression function 

linear.fun <- function(species, exposure, model = basic.model){
  
  temp.data <- merge(scapis[, c("Subject",exposure, model) , with=F], clr.species, by = "Subject")
  
  if(!"data.table" %in% class(temp.data)) setDF(temp.data)
  
  rownames(temp.data) <- temp.data$Subject
  
  
  res <- bplapply(species, function(y){
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
          
          fit <- lm(form,temp.data)
          
          ci <- data.frame(confint(fit, parm = exposure))
          colnames(ci) <- c("LCI", "HCI")
          
          gvif <- car::vif(fit)[exposure,"GVIF" , drop = F]
          
    
          temp.res <- summary(fit)
    
          N = length(temp.res$residuals)
          res <- temp.res$coef[exposure, ]
          colnames(res) <- c("beta", "SE", "t.value","p.value")
    
          df <- data.frame(outcome = y, exposure = rownames(res), res, N = N)
          cbind(df,ci, gvif)
  
    }, BPPARAM = MulticoreParam(16))
  
  res <- do.call(rbind, res)
  
  setDT(res)
  
  res[,q.value := p.adjust(p.value, method = "BH")]
  res[,model := ""]
  
  
}

  # Basic model ####

  
  atbclasses = grep("Class_.*yr", colnames(scapis), value=T)
  
  print(atbclasses)
  
  
  res <- linear.fun(prevalent.species,  exposure = atbclasses )
  res$model <- "basic.model"
  
  fwrite(res, file='results/scapis_species_class.tsv')
  
  
  message("End Basic model ")
  
  t1 = Sys.time()
  print(t1-t0)
  
  #### FULL MODEL #### ----------------------------------------------------
  
  
  res_full <- linear.fun(prevalent.species, 
                    exposure = atbclasses, 
                    model = full.model)
  res_full$model <-  "full.model"
  
  res <- rbind(res, res_full)
  
  fwrite(res, file='results/scapis_species_class.tsv')
  
  
  message("End Full Model")
  
  t1 = Sys.time()
  print(t1-t0)
  
  
