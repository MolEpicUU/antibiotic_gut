# Project Antibiotic Use and the Gut Microbiota

 # Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  
  t0 = Sys.time()

 # This script will investigate the associations of number of previous antibiotics with species abundance
  
 # Sensitivity analysis stratified age and sex
  
# LINEAR REGRESSION ####

  setwd('users/baldanzi/')
  
  # Import data set
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))

  clr.species <- fread('work/simpler_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  

  # Import Basic model
  load('work/simpler_model_revision.Rdata')


# Create the linear regression function 

  linear.fun <- function(species, exposure = atbclasses, model = basic.model, dat = simpler){
    
    temp.data <- merge(dat[, c("SIMPKEY",exposure, model) , with=F], clr.species, by = "SIMPKEY")
    
    if(!"data.table" %in% class(temp.data)) setDF(temp.data)
    
    rownames(temp.data) <- temp.data$Subject
    
    cc <- complete.cases(temp.data)
    Nexp <- temp.data[cc ,.(exposure = names(.SD), Nexp = as.numeric(lapply(.SD, function(x) sum(x>0)))), .SDcols  = exposure]
    
    exp_equal_zero <- Nexp[Nexp<=5, exposure]
    
    exposure <- exposure[!exposure %in% exp_equal_zero]
    
    
    res <- bplapply(species, function(y){
      
      tryCatch({
        
        form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
        
        fit <- lm(form,temp.data)

        
        
        temp.res <- summary(fit)
        
        N = length(temp.res$residuals)
        res <- temp.res$coef[exposure, ]
        colnames(res) <- c("beta", "SE", "t.value","p.value")
        
        data.frame(outcome = y, exposure = rownames(res), res, N = N, message = NA)
        
      }, error = function(e){
        
        df <- data.frame(outcome = y, message = e$message)
        
        return(df)
      })
      
    }, BPPARAM = MulticoreParam(16))
    
    res <- rbindlist(res, fill =T)
    
    setDT(res)
    
    res[,q.value := p.adjust(p.value, method = "BH")]
    res[,model := ""]
    res <- merge(res, Nexp, by = c("exposure"), all=T)
    
  }

  # Full model ####

    
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  # sex stratified
  full.model_sex = full.model[-which(full.model=="sex")]
  
  ressimpler_male <- linear.fun(prevalent.species, model = full.model_sex,  dat = simpler[sex == "M", ])
  ressimpler_male$model = "male"
  ressimpler_male[, Nexp := as.numeric(Nexp)]
  
  fwrite(ressimpler_male, file='results/simpler_species_class_sa_stratified_allantibiotics__male.tsv', sep = "\t")
  rm(ressimpler_male) 
  
  ressimpler_female <- linear.fun(prevalent.species, model = full.model_sex,  dat = simpler[sex == "F", ])
  ressimpler_female$model = "female"
  ressimpler_female[, Nexp := as.numeric(Nexp)]
  
  fwrite(ressimpler_female, file='results/simpler_species_class_sa_stratified_allantibiotics__female.tsv', sep = "\t")
  rm(ressimpler_female)
  
  # No SIMPLER participants with age < 55 
  
  ressimpler_ageabove55 <- linear.fun(prevalent.species, model = full.model,  dat = simpler[age > 55, ])
  ressimpler_ageabove55$model = "ageabove55"
  ressimpler_ageabove55[, Nexp := as.numeric(Nexp)]
  
  fwrite(ressimpler_ageabove55, file='results/simpler_species_class_sa_stratified_allantibiotics__ageabove55.tsv', sep = "\t")
  rm(ressimpler_ageabove55)
  
  
  
  message("End")
  
  t1 = Sys.time()
  print(t1-t0)
  
 
  
  
  
