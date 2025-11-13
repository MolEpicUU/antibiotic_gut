# Project Antibiotic Use and the Gut Microbiota

 # Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  
  t0 = Sys.time()

 # This script will investigate the associations of number of previous antibiotics with species abundance
  
 # Sensitivity analysis stratified age and sex
  
# LINEAR REGRESSION ####

  setwd('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/')
  
  # Import data set
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  
  clr.species <- fread('work/scapis_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  
  
  # Import Basic model
  load('work/scapis_model_revision.Rdata')

# Create the linear regression function 

linear.fun <- function(species, exposure = atbclasses, model = full.model, dat = scapis){
  
  temp.data <- merge(dat[, c("Subject",exposure, model) , with=F], clr.species, by = "Subject")
  
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
          
          #ci <- data.frame(confint(fit, parm = exposure))
          #colnames(ci) <- c("LCI", "HCI")
          
          #gvif <- car::vif(fit)[exposure,"GVIF" , drop = F]
          
    
          temp.res <- summary(fit)
    
          N = length(temp.res$residuals)
          res <- temp.res$coef[exposure, ]
          colnames(res) <- c("beta", "SE", "t.value","p.value")
    
          #df <- data.frame(outcome = y, exposure = rownames(res), res, ci, N = N, gvif, message = NA)
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

  # By Class of antibiotics ####----------------------------------------------
  
  atbclasses = grep("Class_.*yr", colnames(scapis), value=T)
  full.model_sex = full.model[-which(full.model=="Sex")]
 

  # Sex stratified
  
  
  resscapis_male <- linear.fun(prevalent.species, model = full.model_sex,  dat = scapis[Sex == "male", ])
  resscapis_male$model = "male"
  resscapis_male[, Nexp := as.numeric(Nexp)]
  
  fwrite(resscapis_male, file='results/scapis_species_class_sa_stratified_allantibiotics__male.tsv', sep = "\t")
  rm(resscapis_male) 
  
  resscapis_female <- linear.fun(prevalent.species, model = full.model_sex,  dat = scapis[Sex == "female", ])
  resscapis_female$model = "female"
  resscapis_female[, Nexp := as.numeric(Nexp)]
  
  fwrite(resscapis_female, file='results/scapis_species_class_sa_stratified_allantibiotics__female.tsv', sep = "\t")
  rm(resscapis_female)
  
  resscapis_agebelow55 <- linear.fun(prevalent.species, model = full.model,  dat = scapis[age <= 55, ])
  resscapis_agebelow55$model = "agebelow55"
  resscapis_agebelow55[, Nexp := as.numeric(Nexp)]
  
  fwrite(resscapis_agebelow55, file='results/scapis_species_class_sa_stratified_allantibiotics__agebelow55.tsv', sep = "\t")
  rm(resscapis_agebelow55)
  
  resscapis_ageabove55 <- linear.fun(prevalent.species, model = full.model,  dat = scapis[age > 55, ])
  resscapis_ageabove55$model = "ageabove55"
  resscapis_ageabove55[, Nexp := as.numeric(Nexp)]
  
  fwrite(resscapis_ageabove55, file='results/scapis_species_class_sa_stratified_allantibiotics__ageabove55.tsv', sep = "\t")
  rm(resscapis_ageabove55)
  

  
  message("End - all strata")
  
  t1 = Sys.time()
  print(t1-t0)
  
  
  