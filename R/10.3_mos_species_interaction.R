# Project Antibiotic Use and the Gut Microbiota

# Script created to use the MOS data 


# Load packages
  library(data.table)
  library(BiocParallel)
  library(lmerTest)
  
  t0 = Sys.time()
  
  args <- commandArgs(trailingOnly = TRUE)
  # args <- c("Class_FQs_", "sex") 
  
  abx_it <- args[1]
  interaction <- args[2]
  
  cat("\nantibiotic =", abx_it)
  cat("\ninteraction =", interaction)

  
# MIXED MODEL REGRESSION ####

  setwd('users/baldanzi/')

  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  clr.species <- fread('work/mos_clr_species.tsv', na.strings = c("NA", NA, ""),  data.table = F)

  # Import Basic model and species by abundance
  load('work/mos_model_revision.Rdata')
  
  mos[, sex := factor(sex, c(2, 1), c("female", "male"))]


  # Create the linear regression function 

  linear.fun <- function(species, exposure = atbclasses, model , abx = abx_it, it = interaction, dat){
  
    temp.data <- merge(dat[, c("lopnrMOS",exposure, model,"family") , with=F], clr.species, by = "lopnrMOS")
  
    if(!"data.table" %in% class(temp.data)) setDF(temp.data)
  
    rownames(temp.data) <- temp.data$lopnrMOS
    
    
  res <- bplapply(species, function(y){
      
    tryCatch({
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"),
                                    "+", paste(paste0(abx, "*", it), collapse = "+"), "+(1|family)"))
          
          fit <- lmer(form,temp.data)
          
          temp.res <- summary(fit)

    
          N = length(temp.res$residuals)
          res <- temp.res$coef[grepl(":", rownames(temp.res$coef)), , drop=F ]
          colnames(res) <- c("beta", "SE","df", "t.value","p.value")
          
         
          
          # Run linearHypothesis test
          interaction_terms <- grep(":", rownames(temp.res$coef), value = TRUE)
          hypotheses <- paste0(interaction_terms, " = 0")
          lh <- car::linearHypothesis(fit, hypotheses)
          p_value <- lh[[grep("Pr\\(>", colnames(lh), value = TRUE)]][2]
          
          # LRT 
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"), "+(1|family)"))
          fit_reduced <- lmer(form,temp.data)
          plrt <- anova(fit, fit_reduced)
          p_lrtest <- plrt$`Pr(>Chisq)`[2]
          
          
         data.table(outcome = y, exposure = rownames(res), res,  N = N, p_omnibus = p_value, p_lrt = p_lrtest, message = NA)
          
         
      }, error = function(e){
       
         df <- data.table(outcome = y, exposure =NA, beta=NA,SE=NA,df=NA,t.value=NA,p.value=NA,N=NA, message = e$message)
         return(df)
         
      })
      
    }, BPPARAM = MulticoreParam(16))
  
  res <- rbindlist(res, fill=T)
  
  setDT(res)
  

  res[,model := it]
  return(res)

  
}

  
  
  #### Full model #### ----------------------------------------------------

   # By Class of antibiotics #### ------------------------------------------------
  
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
  
  
  if(interaction == "sex"){
    
    cc <- complete.cases(mos[, c("shannon", full.model), with=F])
    to_remove <- mos[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses, by = sex]
    to_remove <- to_remove[Nexp<=5, unique(names)]
    if(any(grepl(abx_it, to_remove))) {message(paste(abx_it , "is too sparse. Exiting.")); quit(save="no", status = 0)} 
    
  }
  
  if(interaction == "age"){
    
    cc <- complete.cases(mos[, c("shannon", full.model), with=F])
    to_remove <- mos[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses]
    to_remove <- to_remove[Nexp<=10, unique(names)]
    if(any(grepl(abx_it, to_remove))) {message(paste(abx_it , "is too sparse. Exiting.")); quit(save="no", status = 0)} 
    
  }

  abx_it <- grep(abx_it, atbclasses, v=T)
  
  # Fitting the model 
  #   y = "shannon"; model = full.model; abx = abx_it; it = interaction; temp.data = mos; exposure = atbclasses
  resmos <- linear.fun(prevalent.species, model = full.model, abx = abx_it, it = interaction, dat = mos)

  fwrite(resmos, file=paste0('results/interaction/mos_species__interaction_', interaction, '_', args[1], '.tsv'), sep = "\t")

  t1 = Sys.time()
  print(t1-t0)
  message("End")
  
  