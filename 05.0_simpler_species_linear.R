# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  

# This script will investigate the associations of number of previous antibiotics with species abundance

  # Import data set
  simpler <- fread("work/simpler_working_dataset.csv", na.strings = c("NA", NA, ""))

  clr.species <- fread('work/simpler_clr_species.tsv', na.strings = c("NA", NA, ""),  data.table = F)

  # Import model and species by abundance
  load('work/simpler_model.Rdata')


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
          cbind(df,ci,  gvif)
  
    }, BPPARAM = MulticoreParam(16))
  
  res <- do.call(rbind, res)
  
  setDT(res)
  
  res[,q.value := p.adjust(p.value, method = "BH")]
  res[,model := ""]
  
  
}

  # basic model ####
  # Number of doses ####------------------------------------------------------

  res <- linear.fun(prevalent.species , exposure = c("N1yr","N1_4yr","N4_8yr"))
  res$model <- "basic.model"
  
  
  fwrite(res, file='results/simpler_species_linear_N1N2N3.tsv')
  

  
  # By Class of antibiotics ####----------------------------------------------
  message("\nClass of antibiotics")
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  
  res <- linear.fun(prevalent.species,  exposure = atbclasses )
  res$model <- "basic.model"
  
  fwrite(res, file='results/simpler_species_linear_class_N1N2N3.tsv')
  
  
  message("End basic Model ")
  
  
  #### full MODEL #### ----------------------------------------------------
  
  res <- linear.fun(prevalent.species , exposure = c("N1yr","N1_4yr","N4_8yr") , model = full.model)
  res$model <-  "full.model"
  
  
  fwrite(res, file='results/simpler_species_linear_N1N2N3_Dis.tsv')


  # By Class of antibiotics #### ------------------------------------------------


   res <- linear.fun(prevalent.species, exposure = atbclasses, model = full.model)
   res$model <-  "full.model"
  
   fwrite(res, file='results/simpler_species_linear_class_N1N2N3_Dis.tsv')

  
  message("End full Model")
  
  
