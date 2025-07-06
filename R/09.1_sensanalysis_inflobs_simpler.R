# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  library(BiocParallel)

  # This script will import the meta-analysis results. 
  # For the results with significant heterogeneity, we will investigate heterogeneity in the individuals studies
  
  setwd('users/baldanzi/')
  
  
  # Import meta-analysis results 
  hetero_species <- readRDS("work/hetero_species.rds")
  
  
  # SIMPLER 
  # Import data
  simpler <- fread("work/simpler_working_dataset_revision.csv", na.strings = c("NA", NA, ""))
  clr.species <- fread('work/simpler_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  load('work/simpler_model_revision.Rdata') # load models 
  

# Create the linear regression function 

  linear.fun <- function(species, exposure, model = full.model){
    
    temp.data <- merge(simpler[, c("SIMPKEY",exposure, model) , with=F], clr.species, by = "SIMPKEY")
    
    if("data.table" %in% class(temp.data)) setDF(temp.data)
    
    rownames(temp.data) <- temp.data$SIMPKEY
    
    res <- bplapply(species, function(y){
      
      form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
      
      fit <- lm(form,temp.data)
      
      dfbetas_fit <-  as.data.frame(dfbetas(fit)[,exposure])
      dfbetas_fit$SIMPKEY <- rownames(dfbetas_fit)
      
      influential.subjects <- apply(dfbetas_fit[,exposure], 2, function(x) dfbetas_fit$SIMPKEY[which.max(abs(x))])
      
      df.infl <- lapply(1:length(influential.subjects), function(x) {
        
        res.infl <- summary(lm(form, temp.data[-which(temp.data$SIMPKEY == influential.subjects[x]),]))
        res.infl <- data.frame(beta_infl = res.infl$coef[names(influential.subjects)[x], c(1)],
                               se_infl = res.infl$coef[names(influential.subjects)[x], c(2)],
                               p.value_inf = res.infl$coef[names(influential.subjects)[x], c(4)],
                               N_inf = length(res.infl$residuals))
        
        res.infl$dfbeta <- dfbetas_fit[dfbetas_fit$SIMPKEY==influential.subjects[x] , names(influential.subjects)[x]]
        return(res.infl)
        
      })
      
      df.infl <- do.call(rbind, df.infl)
      
      
      temp.res <- summary(fit)
      
      N = length(temp.res$residuals)
      res <- temp.res$coef[exposure, ]
      colnames(res) <- c("beta", "SE", "t.value","p.value")
      
      data.frame(outcome = y, exposure = rownames(res), res, N = N, df.infl)
      
    }, BPPARAM = MulticoreParam(16))
    
    res <- do.call(rbind, res)
    
    setDT(res)
    
    res[,q.value := p.adjust(p.value, method = "BH")]
    res[,model := "full.model"]
    
    
  }
  
  message("SIMPLER")
  
  atbclasses = grep("Class_.*yr", colnames(simpler), value=T)
  res <- linear.fun(hetero_species , exposure = atbclasses )
  res$model <- "full.model"
  
  fwrite(res, file='results/simpler_infobs_byclass.tsv', sep = "\t")
  
  message("\nSIMPLER - END")
  
 