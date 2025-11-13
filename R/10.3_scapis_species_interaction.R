# Project Antibiotic Use and the Gut Microbiota


 # Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  
  t0 = Sys.time()
  
  args <- commandArgs(trailingOnly = TRUE)
# args <- c("Class_FQs_1yr", "sex") 

  abx_it <- args[1]
  interaction <- args[2]

 # This script will investigate the associations of number of previous antibiotics with species abundance

# LINEAR REGRESSION ####

  setwd('users/baldanzi/')
  
  # Import data set
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  setnames(scapis, "Sex", "sex")

  clr.species <- fread('work/scapis_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  

  # Import Basic model
  load('work/scapis_model_revision.Rdata')
  full.model[full.model == "Sex"] <- "sex"


# Create the linear regression function 

linear.fun <- function(species, exposure = atbclasses, model = full.model, abx = abx_it, it = interaction, dat = scapis){
  
  temp.data <- merge(dat[, c("Subject",exposure, model) , with=F], clr.species, by = "Subject")
  
  if(!"data.table" %in% class(temp.data)) setDF(temp.data)
  
  rownames(temp.data) <- temp.data$Subject
  
  
  res <- bplapply(species, function(y){
    
    tryCatch({
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"), "+",
                                    paste(paste0(abx, "*", it), collapse = "+")))
          
          fit <- lm(form,temp.data)
          
    
          temp.res <- summary(fit)
    
          N = length(temp.res$residuals)
          res <- temp.res$coef[grepl(":", rownames(temp.res$coef)), , drop=F ]
          colnames(res) <- c("beta", "SE", "t.value","p.value")
          
          

          # Run linearHypothesis test
          interaction_terms <- grep(":", rownames(temp.res$coef), value = TRUE)
          hypotheses <- paste0(interaction_terms, " = 0")
          lh <- linearHypothesis(fit, hypotheses)
          p_value <- lh[[grep("Pr\\(>", colnames(lh), value = TRUE)]][2]
          
          # LRT 
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
          fit_reduced <- lm(form,temp.data)
          plrt <- lmtest::lrtest(fit, fit_reduced)
          p_lrtest <- plrt$`Pr(>Chisq)`[2]
    
          df <- data.frame(outcome = y, exposure = rownames(res), res, N = N, p_omnibus = p_value, p_lrt = p_lrtest, message = NA)
          return(df)
          
      }, error = function(e){
        
        df <- data.frame(outcome = y, message = e$message)
        
        return(df)
      })
  
    }, BPPARAM = MulticoreParam(16))
  
  res <- rbindlist(res, fill =T)
  
  setDT(res)
  
  res[, model := it]
  
  return(res)
  
}

  # Full model ####

  
  atbclasses = grep("Class_.*yr", colnames(scapis), value=T)

  if(interaction == "sex"){
    
    cc <- complete.cases(scapis[, c("shannon", full.model), with=F])
    to_remove <- scapis[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses, by = sex]
    to_remove <- to_remove[Nexp<=5, unique(names)]
    if(any(grepl(abx_it, to_remove))) {message(paste(abx_it , "is too sparse. Exiting.")); quit(save="no", status = 0)} 
    
  }
  
  if(interaction == "age"){
    
    cc <- complete.cases(scapis[, c("shannon", full.model), with=F])
    to_remove <- scapis[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses]
    to_remove <- to_remove[Nexp<=10, unique(names)]
    if(any(grepl(abx_it, to_remove))) {message(paste(abx_it , "is too sparse. Exiting.")); quit(save="no", status = 0)} 
    
  }
  
  abx_it <- grep(abx_it, atbclasses, v=T)

  # Fitting the model 
  
  resscapis <- linear.fun(prevalent.species, model = full.model, abx = abx_it, it = interaction,  dat = scapis)

  
  fwrite(resscapis, file=paste0('results/interaction/scapis_species__interaction_', interaction, '_', args[1], '.tsv'), sep = "\t")
  
  
  message("End - main analysis")
  
  t1 = Sys.time()
  print(t1-t0)
  
  
  