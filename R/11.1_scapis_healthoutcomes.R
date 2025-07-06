# Project Antibiotic Use and the Gut Microbiota

# Script created to use the SCAPIS data 

# This script will investigate the association between species and 
# cardiometabolic markers. 

# Load packages
library(data.table)
library(dplyr)
library(tidyr)
library(BiocParallel)
library(car)

  setwd("nobackup/users/baldanzi/atb_gut/")
  
  meta_res <- fread("results/meta_species_class_clean.tsv")
  
  meta_res <- meta_res[grepl("lincos|FQs|BetaR", exposure) & model == "full.model" & q.value<0.05, ]
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
  
  clr.species <- fread('work/scapis_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  clr.species <- clr.species[, c("Subject", species_overlap), with=F]
  
  
  # Import Basic model
  load('work/scapis_model_revision.Rdata')
  
  
  setnames(scapis, c("WaistHip", "HdlFormattedResult", "LdlFormattedResult", "TgFormattedResult", "GlucoseFormattedResult", "Hba1cFormattedResult", 
                     "CreatinineFormattedResult", "CholesterolFormattedResult", "CrpFormattedResult"),
                    c("WHR", "HDL", "LDL", "TG", "Glucose", "HbA1c", "TC", "CRP"))
  
  healthoutcomes <- c("WHR", "HDL", "LDL", "TG", "Glucose", "HbA1c", "TC", "CRP")
  
  
  
  # Spearman correlation function 
  spearman.fun <- function(outcome, species =  species_overlap ){
    
    res <- bplapply(species, function(y){
      
      model <- c("age", "Sex", "site_plate", "placebirth", "smokestatus", "education")
      
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
      
    }, BPPARAM = MulticoreParam(4))
    
    res <- do.call(rbind, res)
    
    setDT(res)
    
    return(res)
    
  }
  
  
  # Run model ####
  setDF(scapis)
  res <- lapply(c("BMI", healthoutcomes), spearman.fun)
  res <- rbindlist(res)
  res[, q.value := p.adjust(p.value, method="BH")]
  
  res[q.value<0.05,]
  
  fwrite(res, file = "results/scapis_healthoutcomes_spearman.tsv", sep = "\t")
  
  message("End")
  
  
  
  
  
  
  