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


  # Create the linear regression function 

  linear.fun <- function(species, exposure = atbclasses, model , dat){
  
    temp.data <- merge(dat[, c("lopnrMOS",exposure, model,"family") , with=F], clr.species, by = "lopnrMOS")
  
    if(!"data.table" %in% class(temp.data)) setDF(temp.data)
  
    rownames(temp.data) <- temp.data$lopnrMOS
    
    cc <- complete.cases(temp.data)
    Nexp <- temp.data[cc ,.(exposure = names(.SD), Nexp = as.numeric(lapply(.SD, function(x) sum(x>0)))), .SDcols  = exposure]
    
    exp_equal_zero <- Nexp[Nexp<=5, exposure]
    
    exposure <- exposure[!exposure %in% exp_equal_zero]
    
    
  res <- bplapply(species, function(y){
      
    tryCatch({
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"),"+(1|family)"))
          
          fit <- lmer(form,temp.data)
          
          temp.res <- summary(fit)
          exp_present <- grep("^Class", rownames(temp.res$coef), v=T)
    
          N = length(temp.res$residuals)
          res <- temp.res$coef[exp_present, ]
          colnames(res) <- c("beta", "SE","df", "t.value","p.value")
          

          data.table(outcome = y, exposure = rownames(res), res,  N = N, message = NA)
         
      }, error = function(e){
       
        df <- data.table(outcome = y, exposure =NA, beta=NA,SE=NA,df=NA,t.value=NA,p.value=NA,N=NA, message = e$message)
         return(df)
         
      })
      
    }, BPPARAM = MulticoreParam(16))
  
  res <- rbindlist(res, fill=T)
  
  setDT(res)
  
  res[!is.na(p.value), q.value := p.adjust(p.value, method = "BH")]
  res[,model := ""]
  res <- merge(res, Nexp, by = c("exposure"), all=T)
  
}

  
  
  #### Full model #### ----------------------------------------------------

   # By Class of antibiotics #### ------------------------------------------------
  

  
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
  full.model_sex = full.model[-which(full.model=="sex")]
   
  message("sex stratified")
  message("Male")

  
  
  resmos_male <- linear.fun(prevalent.species, model = full.model_sex,  dat = mos[sex == 1, ])
  resmos_male$model = "male"
  resmos_male[, Nexp := as.numeric(Nexp)]
  
  fwrite(resmos_male, file='results/mos_species_class_sa_stratified_allantibiotics__male.tsv', sep = "\t")
  rm(resmos_male) 
  
  message("Female")
  resmos_female <- linear.fun(prevalent.species, model = full.model_sex,  dat = mos[sex == 2, ])
  resmos_female$model = "female"
  resmos_female[, Nexp := as.numeric(Nexp)]
  
  fwrite(resmos_female, file='results/mos_species_class_sa_stratified_allantibiotics__female.tsv', sep = "\t")
  rm(resmos_female)
  

  message("\nAge stratified")
   
  message("Below 55")

  resmos_agebelow55 <- linear.fun(prevalent.species, model = full.model,  dat = mos[age <= 55, ])
  resmos_agebelow55$model = "agebelow55"
  resmos_agebelow55[, Nexp := as.numeric(Nexp)]
  
  fwrite(resmos_agebelow55, file='results/mos_species_class_sa_stratified_allantibiotics__agebelow55.tsv', sep = "\t")
  rm(resmos_agebelow55)
  
  message("Above 55")
  
  resmos_ageabove55 <- linear.fun(prevalent.species, model = full.model,  dat = mos[age > 55, ])
  resmos_ageabove55$model = "ageabove55"
  resmos_ageabove55[, Nexp := as.numeric(Nexp)]
  
  fwrite(resmos_ageabove55, file='results/mos_species_class_sa_stratified_allantibiotics__ageabove55.tsv', sep = "\t")
  rm(resmos_ageabove55)
  

  t1 = Sys.time()
  print(t1-t0)
  message("End")
  
  