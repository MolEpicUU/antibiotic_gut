# Project Antibiotic Use and the Gut Microbiota

# Script created to use the SIMPLER data 


# Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  
  t0 = Sys.time()
  

# Interaction term analysis - alpha diversity 
  

  setwd('users/baldanzi/')

  # Import data set
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  simpler[, sex := factor(sex, c("F", "M"), c("female", "male"))]

  # Import models and species by abundance
  load('work/simpler_model_revision.Rdata')


  # Create the linear regression function 

  linear.fun <- function(abx , exposure = atbclasses, model = full.model,  it , dat){
  
    res <- lapply(c("shannon", "richness", "invsimpson"), function(y){
      
    tryCatch({
      
          abx_it <- grep(abx, atbclasses, v=T)
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"), 
                                    "+", paste(paste0(abx_it, "*", it), collapse = "+")))
          
          fit <- lm(form, data = dat)
          
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
          fit_reduced <- lm(form, data = dat)
          plrt <- lmtest::lrtest(fit, fit_reduced)
          p_lrtest <- plrt$`Pr(>Chisq)`[2]
          
          
          data.table(outcome = y, exposure = rownames(res), res,  N = N, p_omnibus = p_value, p_lrt = p_lrtest, message = NA)
          
         
      }, error = function(e){
       
         df <- data.table(outcome = y, exposure =NA, beta=NA,SE=NA,t.value=NA,p.value=NA,N=NA, message = e$message)
         return(df)
         
      })
      
    })
  
  res <- rbindlist(res, fill=T)
  
  setDT(res)
  

  res[,model := it]
  return(res)

  
}

  
  
  #### Full model #### ----------------------------------------------------
  
  # AGE
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  classes <- c("Class_Peni_Ext",  "Class_Peni_BetaS", "Class_lincosamides", "Class_SMZTMP",  "Class_FQs" ,  "Class_Peni_Comb", 
               "Class_NIT", "Class_Peni_BetaR",  "Class_TCLs", "Class_macrolides",  "Class_cephalosporins")
  
    cc <- complete.cases(simpler[, c("shannon", full.model), with=F])
    to_remove <- simpler[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses]
    to_remove <- to_remove[Nexp<=10, unique(names)]
    
  to_remove <- unique(gsub("_\\d.*", "", to_remove))
  
  if(length(to_remove)>0) classes <- classes[!classes %in% to_remove]
  
  ressimpler_age <- lapply(classes, linear.fun, model = full.model, it = "age", dat = simpler) # BPPARAM = MulticoreParam(16))
  ressimpler_age <- rbindlist(ressimpler_age)
  
  # SEX 
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  classes <- c("Class_Peni_Ext",  "Class_Peni_BetaS", "Class_lincosamides", "Class_SMZTMP",  "Class_FQs" ,  "Class_Peni_Comb", 
               "Class_NIT", "Class_Peni_BetaR",  "Class_TCLs", "Class_macrolides",  "Class_cephalosporins")
  
    cc <- complete.cases(simpler[, c("shannon", full.model), with=F])
    to_remove <- simpler[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses, by = sex]
    to_remove <- to_remove[Nexp<=5, unique(names)]
    to_remove <- unique(gsub("_\\d.*", "", to_remove))
  
  print(to_remove)
  
  if(length(to_remove)>0) classes <- classes[!classes %in% to_remove]
  
  ressimpler_sex <- lapply(classes, linear.fun, model = full.model, it = "sex", dat = simpler) # BPPARAM = MulticoreParam(16))
  ressimpler_sex <- rbindlist(ressimpler_sex)
  
  ressimpler <- rbind(ressimpler_age, ressimpler_sex)
  ressimpler[, cohort:="SIMPLER"]

  fwrite(ressimpler, file='results/simpler_alpha__interaction.tsv', sep = "\t")
  fwrite(ressimpler, file='/proj/nobackup/simp2023007/wharf/baldanzi/baldanzi-simp2023007/atbgut/results/simpler_alpha__interaction.tsv', sep = "\t")

  t1 = Sys.time()
  print(t1-t0)
  message("End")
  
  