# Project Antibiotic Use and the Gut Microbiota

# Script created to use the MOS data 

# Load packages
  library(data.table)
  library(BiocParallel)
  library(lmerTest)
  library(car)
  
  t0 = Sys.time()
  

# Interaction term analysis - alpha diversity 
  
# MIXED MODEL REGRESSION ####

  setwd('users/baldanzi/')

  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))

  # Import Basic model and species by abundance
  load('work/mos_model_revision.Rdata')
  
  mos[, sex := factor(sex, c(2, 1), c("female", "male"))]


  # Create the linear regression function 

  linear.fun <- function(abx, exposure = atbclasses, model = full.model,  it, dat){
  
    res <- lapply(c("shannon", "richness", "invsimpson"), function(y){
      
    tryCatch({
    
          abx_it <- grep(abx, atbclasses, v=T)
      
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"),
                                    "+", paste(paste0(abx_it, "*", it), collapse = "+"), "+(1|family)"))
          
          fit <- lmer(form, data = dat)
          
          temp.res <- summary(fit)

    
          N = length(temp.res$residuals)
          res <- temp.res$coef[grepl(":", rownames(temp.res$coef)), , drop=F ]
          colnames(res) <- c("beta", "SE","df", "t.value","p.value")
          
          # Run linearHypothesis test
          interaction_terms <- grep(":", rownames(temp.res$coef), value = TRUE)
          hypotheses <- paste0(interaction_terms, " = 0")
          lh <- linearHypothesis(fit, hypotheses)
          p_value <- lh[[grep("Pr\\(>", colnames(lh), value = TRUE)]][2]
          
          # LRT 
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"), "+(1|family)"))
          fit_reduced <- lmer(form, data = dat)
          plrt <- anova(fit, fit_reduced)
          p_lrtest <- plrt$`Pr(>Chisq)`[2]
          
          
         data.table(outcome = y, exposure = rownames(res), res,  N = N, p_omnibus = p_value, p_lrt = p_lrtest, message = NA)
          
         
      }, error = function(e){
       
         df <- data.table(outcome = y, exposure =NA, beta=NA,SE=NA,df=NA,t.value=NA,p.value=NA,N=NA, message = e$message)
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
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
  classes <- c("Class_Peni_Ext",  "Class_Peni_BetaS", "Class_lincosamides", "Class_SMZTMP",  "Class_FQs" ,  "Class_Peni_Comb", 
               "Class_NIT", "Class_Peni_BetaR",  "Class_TCLs", "Class_macrolides",  "Class_cephalosporins")
  
    cc <- complete.cases(mos[, c("shannon", full.model), with=F])
    to_remove <- mos[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses]
    to_remove <- to_remove[Nexp<=10, unique(names)]
  
  to_remove <- unique(gsub("_\\d.*", "", to_remove))
  print(to_remove)
  
  if(length(to_remove)>0) classes <- classes[!classes %in% to_remove]
  
  resmos_age <- bplapply(classes, linear.fun, model = full.model, it = "age", dat = mos, BPPARAM = MulticoreParam(16))
  resmos_age <- rbindlist(resmos_age)
  
  # SEX 
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
  classes <- c("Class_Peni_Ext",  "Class_Peni_BetaS", "Class_lincosamides", "Class_SMZTMP",  "Class_FQs" ,  "Class_Peni_Comb", 
               "Class_NIT", "Class_Peni_BetaR",  "Class_TCLs", "Class_macrolides",  "Class_cephalosporins")
  
    cc <- complete.cases(mos[, c("shannon", full.model), with=F])
    to_remove <- mos[cc , .(names = names(.SD), Nexp= colSums(.SD>0)) , .SDcols = atbclasses, by = sex]
    to_remove <- to_remove[Nexp<=5, unique(names)]
    atbclasses <- atbclasses[! atbclasses %in% to_remove]
    
  to_remove <- unique(gsub("_\\d.*", "", to_remove))
  print(to_remove)
    
  if(length(to_remove)>0) classes <- classes[!classes %in% to_remove]
  
  resmos_sex <- bplapply(classes, linear.fun, model = full.model, it = "sex", dat = mos, BPPARAM = MulticoreParam(16))
  resmos_sex <- rbindlist(resmos_sex, fill = T)
  
  resmos <- rbind(resmos_age, resmos_sex)
  resmos[, cohort := "MOS"]

  fwrite(resmos, file='results/mos_alpha__interaction.tsv', sep = "\t")

  t1 = Sys.time()
  print(t1-t0)
  message("End")
  
  