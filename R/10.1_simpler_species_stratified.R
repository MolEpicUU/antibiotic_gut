# Project Antibiotic Use and the Gut Microbiota

 # Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  
  t0 = Sys.time()

  # SIMPLER
  
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
  
  rownames(temp.data) <- temp.data$SIMPKEY
  
  
  res <- bplapply(species, function(y){
    
    tryCatch({
    
          cc <- complete.cases(temp.data)
          Nexp <- temp.data[cc ,.(exposure = names(.SD), outcome = y, Nexp = lapply(.SD, function(x) sum(x>0))), .SDcols  = exposure]
          
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
          
          fit <- lm(form,temp.data)
          
          ci <- data.frame(confint(fit, parm = exposure))
          colnames(ci) <- c("LCI", "HCI")
          
          gvif <- car::vif(fit)[exposure,"GVIF" , drop = F]
          
    
          temp.res <- summary(fit)
    
          N = length(temp.res$residuals)
          res <- temp.res$coef[exposure, ]
          colnames(res) <- c("beta", "SE", "t.value","p.value")
    
          df <- data.frame(outcome = y, exposure = rownames(res), res, ci, N = N, gvif, message = NA)
          merge(df, Nexp, by = c("exposure", "outcome"), all.y=T)
          
      }, error = function(e){
        
        df <- data.frame(outcome = y, exposure = NA, beta=NA, SE=NA,t.value=NA,p.value=NA,LCI=NA,HCI=NA,N=NA,GVIF=NA,Nexp=NA, message = e$message)
        
        return(df)
      })
  
    }, BPPARAM = MulticoreParam(16))
  
  res <- rbindlist(res, fill =T)
  
  setDT(res)
  
  res[,q.value := p.adjust(p.value, method = "BH")]
  res[,model := ""]
  
  
}

  # Full model ####

  # By Class of antibiotics ####----------------------------------------------
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  
  # sex stratified
  full.model_sex = full.model[-which(full.model=="sex")]
  
  ressimpler_male <- linear.fun(prevalent.species, model = full.model_sex,  dat = simpler[sex == "M", ])
  ressimpler_male$model = "male"
  
  ressimpler_female <- linear.fun(prevalent.species, model = full.model_sex,  dat = simpler[sex == "F", ])
  ressimpler_female$model = "female"
  
  
  ressimpler_ageabove55 <- linear.fun(prevalent.species, model = full.model,  dat = simpler[age > 55, ])
  ressimpler_ageabove55$model = "age > 55 years"
  
  
  ressimpler_age_sex <- rbind(ressimpler_male, ressimpler_female, ressimpler_ageabove55)
  ressimpler_age_sex$cohort = "SIMPLER"
  
  fwrite(ressimpler_age_sex, file='results/simpler_species_sa_stratified.tsv')
  
  
  message("End")
  
  t1 = Sys.time()
  print(t1-t0)
  
  
