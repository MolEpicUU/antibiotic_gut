# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  library(BiocParallel)
  library(lmerTest)

  # This script will import the meta-analysis results. 
  # For the results with significant heterogeneity (Cochrans p-value <0.05), we will investigate heterogeneity in the individuals studies
  
  setwd('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/')
  
  
  # Import meta-analysis results 
  resmeta <- fread('results/meta_species_class.tsv')
  
  hetero_species <- resmeta[Qpval<.05 & model == "full.model" & q.value<.05, unique(outcome)] # 179 species 
  saveRDS(hetero_species, file = "work/hetero_species.rds")
  saveRDS(hetero_species, file = "/proj/sens2019512/nobackup/wharf/baldanzi/baldanzi-sens2019512/atbgut/hetero_species.rds")
  
  
  # SCAPIS 
  # Import data
  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  clr.species <- fread('work/scapis_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  load('work/scapis_model_revision.Rdata') # load models 
  

# Create the linear regression function 

  linear.fun <- function(species, exposure, model = full.model){
    
    temp.data <- merge(scapis[, c("Subject",exposure, model) , with=F], clr.species, by = "Subject")
    
    if("data.table" %in% class(temp.data)) setDF(temp.data)
    
    rownames(temp.data) <- temp.data$Subject
    
    res <- bplapply(species, function(y){
      
      form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
      
      fit <- lm(form,temp.data)
      
      dfbetas_fit <-  as.data.frame(dfbetas(fit)[,exposure])
      dfbetas_fit$Subject <- rownames(dfbetas_fit)
      
      influential.subjects <- apply(dfbetas_fit[,exposure], 2, function(x) dfbetas_fit$Subject[which.max(abs(x))])
      
      df.infl <- lapply(1:length(influential.subjects), function(x) {
        
        res.infl <- summary(lm(form, temp.data[-which(temp.data$Subject == influential.subjects[x]),]))
        res.infl <- data.frame(beta_infl = res.infl$coef[names(influential.subjects)[x], c(1)],
                               se_infl = res.infl$coef[names(influential.subjects)[x], c(2)],
                               p.value_inf = res.infl$coef[names(influential.subjects)[x], c(4)],
                               N_inf = length(res.infl$residuals))
        
        res.infl$dfbeta <- dfbetas_fit[dfbetas_fit$Subject==influential.subjects[x] , names(influential.subjects)[x]]
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
  
  message("SCAPIS")
  
  atbclasses = grep("Class_.*yr", colnames(scapis), value=T)
  res <- linear.fun(hetero_species , exposure = atbclasses )
  res$model <- "full.model"
  
  fwrite(res, file='results/scapis_infobs_byclass.tsv', sep = "\t")
  
  message("SCAPIS - end")
  
  ###############################
  ## MOS ########################
  ###############################
  
  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  clr.species <- fread('work/mos_clr_species.tsv', na.strings = c("NA", NA, ""),  data.table = F)
  load('work/mos_model_revision.Rdata')
  
  
  # Create the mix-model function 
  
  mix.fun <- function(species, exposure, model = full.model){
    
    temp.data <- merge(mos[, c("lopnrMOS",exposure, model,"family") , with=F], clr.species, by = "lopnrMOS")
    
    if(!"data.table" %in% class(temp.data)) setDF(temp.data)
    
    rownames(temp.data) <- temp.data$lopnrMOS
    
    
    res <- bplapply(species, function(y){
      
      
      form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
      formmix <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+"),"+(1|family)"))
      
      fit <- lm(form,temp.data)
      
      dfbetas_fit <-  as.data.frame(dfbetas(fit)[,exposure])
      dfbetas_fit$lopnrMOS <- rownames(dfbetas_fit)
      
      influential.subjects <- apply(dfbetas_fit[,exposure], 2, function(x) dfbetas_fit$lopnrMOS[which.max(abs(x))])
      
      df.infl <- lapply(1:length(influential.subjects), function(x) {
        
        res.infl <- summary(lmer(formmix, temp.data[-which(temp.data$lopnrMOS == influential.subjects[x]),]))
        res.infl <- data.frame(beta_infl = res.infl$coef[names(influential.subjects)[x], c(1)],
                               se_infl = res.infl$coef[names(influential.subjects)[x], c(2)],
                               p.value_inf = res.infl$coef[names(influential.subjects)[x], c(5)],
                               N_inf = length(res.infl$residuals))
        
        res.infl$dfbeta <- dfbetas_fit[dfbetas_fit$lopnrMOS==influential.subjects[x] , names(influential.subjects)[x]]
        rownames(res.infl) <- names(influential.subjects)[x]
        return(res.infl)
        
      })
      
      df.infl <- do.call(rbind, df.infl)
      
      
      temp.res <- summary(fit)
      
      N = length(temp.res$residuals)
      res <- temp.res$coef[exposure, ]
      colnames(res) <- c("beta", "SE","t.value","p.value")

      
      data.frame(outcome = y, exposure = rownames(res), res, N = N, df.infl)
      
      
    }, BPPARAM = MulticoreParam(16))
    
    res <- do.call(rbind, res)
    
    setDT(res)
    
    res[,q.value := p.adjust(p.value, method = "BH")]
    res[,model := ""]
    
    
  }
  
  message("MOS")

  
  atbclasses = grep("Class_.*yr", colnames(mos), value=T)
  res <- mix.fun(hetero_species , exposure = atbclasses )
  res$model <- "full.model"

  fwrite(res, file='results/mos_infobs_byclass.tsv', sep = "\t")
  
  message("MOS - END")
  
  
  
  
  