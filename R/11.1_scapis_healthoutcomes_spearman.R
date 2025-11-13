# Project Antibiotic Use and the Gut Microbiota

# Script created to use the SCAPIS data 

# Correlation between species associated with antibiotic use and cardiometabolic biomarkers. 

# Load packages
library(data.table)
library(dplyr)
library(tidyr)
library(BiocParallel)
library(ppcor)

setwd('nobackup/users/baldanzi')
  
  meta_res <- fread("results/meta_species_class_clean.tsv")
  
  meta_res <- meta_res[grepl("lincos|FQs|BetaR", exposure) & model == "full.model" & q.value < 0.05, ]
  meta_res <- meta_res %>% separate(exposure, into = c("antibiotic", "period"), sep = "_(?=\\d)", extra = "merge", remove = F)
  setDT(meta_res)
  
  species_overlap <- lapply(unique(meta_res$antibiotic), function(abx){
   
    sp <- meta_res[antibiotic == abx, ][["outcome"]] 
    return(sp)
    
  })
  
  
  species_overlap <- Reduce(x = species_overlap, f = intersect)
  
  length(species_overlap)


  
  # Import data set
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  scapis <- scapis[diabd == "no"]  # Remove 743 individuals with diabetes. 
  
  # non-HDL-c
  scapis_species <- fread('../SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_mgs_relative_abundances_v1.0.tsv', data.table=F)
  setnames(scapis_species, "scapis_id", "Subject")
  scapis_species <- scapis_species[, c("Subject", species_overlap)]
  scapis_species <- scapis_species[scapis_species$Subject %in% scapis$Subject, ]
  
  
  # Import Basic model
  load('work/scapis_model_revision.Rdata')
  
  
  setnames(scapis, c("WaistHip", "HdlFormattedResult", "LdlFormattedResult", "TgFormattedResult", 
                     "Hba1cFormattedResult",  
                     "CholesterolFormattedResult", "CrpFormattedResult", "SBP_Mean"),
                    c("WHR", "HDL", "LDL", "TG","HbA1c", "TC", "CRP", "SBP"))
  
  scapis[, "non-HDL" := TC - HDL]
  
  healthoutcomes <- c("WHR", "TG", "non-HDL", "HbA1c", "SBP", "CRP")
  
  
  scapis <- merge(scapis[, c("Subject", full.model,  healthoutcomes) , with=F], scapis_species, by = "Subject")
  
  if(!"data.table" %in% class(scapis)) setDF(scapis)
  
  rownames(scapis) <- scapis$Subject
  
  
  # Spearman correlation function 
  spearman.fun <- function(outcome, species =  species_overlap , cov){
    
    res <- bplapply(species, function(y, model = cov){
      
      if(outcome!="BMI"){
        model <- c(model, "BMI")
      } 
      
      cc <- complete.cases(scapis[, c(outcome, y, model)] )
      
      temp.data <- scapis[cc, ]
      
      x1 <- temp.data[[outcome]]
      x2 <- temp.data[[y]]
      cov <- model.matrix(~.  , temp.data[, model])[, -1]
      
      res=pcor.test(x=x1,
                    y=x2,
                    z=cov,
                    method = "spearman")
      
      temp <- data.frame(species = y, outcome = outcome ,
                         rho=res$estimate,
                         p.value=res$p.value,
                         N=res$n,
                         method = res$Method,
                         covariates=paste(model,collapse = "+"))
      
      
      return(temp)
      
    }, BPPARAM = MulticoreParam(16))
    
    res <- do.call(rbind, res)
    
    setDT(res)
    
    return(res)
    
  }
  
  
  # Run model ####
  setDF(scapis)
  covariates <- c("age", "Sex", "site_plate", "placebirth", "smokestatus", "education")
  res <- lapply(c("BMI", healthoutcomes), spearman.fun, cov = covariates)
  res <- rbindlist(res)
  res[, q.value := p.adjust(p.value, method="BH")]
  
  res[q.value<0.05,]
  
  fwrite(res, file = "results/scapis_healthoutcomes_spearman.tsv", sep = "\t")
  
  message("End - all")
  
   # Run model - males ####
  scapis_backup <- scapis 
  scapis <- scapis_backup[scapis_backup$Sex=="male", ]  

  covariates <- c("age", "site_plate", "placebirth", "smokestatus", "education")
  res <- lapply(c("BMI", healthoutcomes), spearman.fun, cov = covariates)
  res <- rbindlist(res)
  res[, q.value := p.adjust(p.value, method="BH")]

  res[q.value<0.05,]

  fwrite(res, file = "results/scapis_healthoutcomes_spearman_males.tsv", sep = "\t")

  message("End - males")
  
  
  
  # Run model - females ####
  scapis <- scapis_backup[scapis_backup$Sex=="female", ]


  res <- lapply(c("BMI", healthoutcomes), spearman.fun, cov = covariates)
  res <- rbindlist(res)
  res[, q.value := p.adjust(p.value, method="BH")]

  res[q.value<0.05,]

  fwrite(res, file = "results/scapis_healthoutcomes_spearman_females.tsv", sep = "\t")

  message("End - females")
