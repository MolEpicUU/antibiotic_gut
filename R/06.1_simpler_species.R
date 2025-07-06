# Project Antibiotic Use and the Gut Microbiota

# Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  
  t0 = Sys.time()

# This script will investigate the associations of number of previous antibiotics with species abundance
  
# LINEAR REGRESSION ####
  
  
  setwd('users/baldanzi/')

  # Import data set
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))

  clr.species <- fread('work/simpler_clr_species.tsv', na.strings = c("NA", NA, ""),  data.table = F)

  # Import basic model and species by abundance
  load('work/simpler_model_revision.Rdata')


# Create the linear regression function 

linear.fun <- function(species, exposure, model = basic.model){
  
  temp.data <- merge(simpler[, c("SIMPKEY",exposure, model) , with=F], clr.species, by = "SIMPKEY")
  
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
          cbind(df,ci,  gvif) #df.infl,
  
    }, BPPARAM = MulticoreParam(16))
  
  res <- do.call(rbind, res)
  
  setDT(res)
  
  res[,q.value := p.adjust(p.value, method = "BH")]
  res[,model := ""]
  
  
}

  # Basic model ####
  
  message("\nClass of antibiotics")
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  
  res <- linear.fun(prevalent.species ,  exposure = atbclasses )
  res$model <- "basic.model"
  
  fwrite(res, file='results/simpler_species_class.tsv')
  
  
  message("End Basic Model ")
  
  t1 = Sys.time()
  print(t1-t0)
  
  # FULL MODEL #### ----------------------------------------------------
  

   res_full <- linear.fun(prevalent.species, exposure = atbclasses, model = full.model)
   res_full$model <-  "full.model"
   
   res <- rbind(res, res_full)
  
   fwrite(res, file='results/simpler_species_class.tsv')
  
  
  t1 = Sys.time()
  print(t1-t0)
  
  message("End Full Model")
  
  
