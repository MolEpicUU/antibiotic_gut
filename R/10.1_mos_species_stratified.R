# Project Antibiotic Use and the Gut Microbiota

# Script created to use the MOS data 

# Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  library(lmerTest)
  
  t0 = Sys.time()

# This script will investigate the associations of number of previous antibiotics with species abundance stratified 
# by sex and age
  
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
    
    exp_equal_zero <- Nexp[Nexp<4, exposure]
    
    exposure <- exposure[!exposure %in% exp_equal_zero]
    
    
  res <- bplapply(species, function(y){
      
    tryCatch({
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"),"+(1|family)"))
          
          fit <- lmer(form,temp.data)
          
          temp.res <- summary(fit)
          exp_present <- grep("^Class", rownames(temp.res$coef), v=T)
          
          gvif <- car::vif(fit)[exp_present, "GVIF", drop = F]
    
          N = length(temp.res$residuals)
          res <- temp.res$coef[exp_present, ]
          colnames(res) <- c("beta", "SE","df", "t.value","p.value")
          
          listCI <- tryCatch({
              suppressMessages(ci <- confint(fit,parm=exp_present))
              LCI = ci[,1]
              HCI = ci[,2]
              
              data.frame(LCI = LCI, HCI = HCI, message = NA)
              
            }, error = function(e){
              
              data.frame(LCI = NA, HCI = NA, message = e$message)
              
            })
          
         data.table(outcome = y, exposure = rownames(res), res, LCI=listCI$LCI, HCI=listCI$HCI, N = N, gvif, message = listCI$message)
          
         
      }, error = function(e){
       
         df <- data.table(outcome = y, exposure =NA, beta=NA,SE=NA,df=NA,t.value=NA,p.value=NA,LCI=NA,HCI=NA,N=NA,GVIF=NA, message = e$message)
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

  message("Sex stratified")
  message("Male")
  
  atbclasses <-  grep("Class_.*yr", colnames(mos), value=T)
  atbclasses <-  atbclasses[!grepl("NIT", atbclasses)] 
  atbclasses <- grep("Peni_Comb", atbclasses , value=T, invert=T) 
  atbclasses <- grep("SMZTMP", atbclasses , value=T, invert=T) 
   

  full.model_sex = full.model[-which(full.model=="sex")]
  resmos_male <- linear.fun(prevalent.species, model = full.model_sex,  dat = mos[sex == 1, ])
  resmos_male$model = "male"

  resmos_male[, Nexp := as.numeric(Nexp)]

  message("Female")
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
  resmos_female <- linear.fun(prevalent.species, model = full.model_sex,  dat = mos[sex == 2, ])
  resmos_female$model = "female"

  fwrite(rbindlist(list(resmos_male, resmos_female), fill=T) , file='results/mos_species_class_sa_stratified.tsv', sep = "\t")

  message("\nAge stratified")
  
  message("Below 55")
  
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
  atbclasses <- grep("Peni_Comb", atbclasses , value=T, invert=T) 
  
  resmos_agebelow55 <- linear.fun(prevalent.species, model = full.model,  dat = mos[age <= 55, ])
  resmos_agebelow55$model = "age <= 55 years"

  fwrite(rbindlist(list(resmos_male, resmos_female, resmos_agebelow55), fill=T) , file='results/mos_species_class_sa_stratified.tsv', sep = "\t")

  message("Above 55")
  
  atbclasses <-  grep("Class_.*yr", colnames(mos), value=T)
  
  resmos_ageabove55 <- linear.fun(prevalent.species, model = full.model,  dat = mos[age > 55 & !lopnrMOS %in% ind_to_exclude, ])
  resmos_ageabove55$model = "age > 55 years"
  
   
   resmos_age_sex <- rbindlist(list(resmos_male, resmos_female, resmos_agebelow55, resmos_ageabove55), fill =T)
   resmos_age_sex$cohort <- "MOS"

  fwrite(resmos_age_sex, file='results/mos_species_class_sa_stratified.tsv', sep = "\t")

  t1 = Sys.time()
  print(t1-t0)
  message("End")
  
