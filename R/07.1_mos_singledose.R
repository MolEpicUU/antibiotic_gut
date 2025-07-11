# Project Antibiotic Use and the Gut Microbiota
  
  # MOS 
  
  # This script will investigate the associations a single previous antibiotics course with species abundance
  
  # Load packages
library(data.table)
  library(BiocParallel)
  library(car)
  library(lmerTest)
  
  t0 = Sys.time()
  
  
  # MIXED MODEL REGRESSION ####
  
  setwd('nobackup/users/baldanzi/atb_gut/')
  
  # Import data set
  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))
  registerdata   <-  fread('nobackup/users/baldanzi/atb_gut/Data/MOS/MOS_Halftime_n2644_Prescribed_Drug_register_until_2019.tsv')
  clr.species <- fread('work/mos_clr_species.tsv', na.strings = c("NA", NA, ""),  data.table = F)
  
  # Import model and species by abundance
  load('work/mos_model_revision.Rdata')
  
  # Filter register data
  registerdata <- registerdata[grep("^J01", ATC),]
  registerdata <- merge(mos[, .(lopnrMOS, Visit1, Visit2)], registerdata, by="lopnrMOS")
  
  # Atb after 
  registerdata_N0 <- registerdata[lopnrMOS %in% mos$lopnrMOS[mos$N_all==0], ]
  registerdata_N0 <- registerdata_N0[edatum>Visit1 & edatum<(Visit2+365.25),]
  atb_after_participants <- registerdata_N0[,.N, by = lopnrMOS][N==1, lopnrMOS]
  
  # Atb before 
  registerdata_N1 <- registerdata[lopnrMOS %in% mos$lopnrMOS[mos$N_all==1], ]
  registerdata_N1 <- registerdata_N1[edatum<Visit1,]
  registerdata_N1 <- registerdata_N1[which(Visit1-edatum<8*365.2),]
  all(registerdata_N1[,.N, by = lopnrMOS][,N]==1)
  
  class_atb_single <- registerdata_N1[ , .(lopnrMOS, ATC)]
  class_atb_single[grep("^J0",ATC), class:="other"]
  class_atb_single[grep("^J01D[B,C,D,E]",ATC), class:="Class_cephalosporins"]
  class_atb_single[grep("^J01FA",ATC), class:="Class_macrolides"]
  class_atb_single[grep("^J01FF",ATC), class:="Class_lincosamides"]
  class_atb_single[grep("^J01MA",ATC), class:="Class_FQs"]
  class_atb_single[grep("^J01A",ATC), class:="Class_TCLs"]
  class_atb_single[grep("^J01E",ATC), class:="Class_SMZTMP"]
  class_atb_single[grep("^J01XE",ATC), class:="Class_NIT"]
  class_atb_single[grep("^J01CA",ATC), class:="Class_Peni_Ext"]
  class_atb_single[grep("^J01CE",ATC), class:="Class_Peni_BetaS"]
  class_atb_single[grep("^J01CF",ATC), class:="Class_Peni_BetaR"]
  class_atb_single[grep("^J01CR",ATC), class:="Class_Peni_Comb"]
  
  no_atb <- mos[N_all==0, .(lopnrMOS, class = last_atb_period)]
  class_atb_single <- rbind(class_atb_single[, -"ATC"], no_atb)
  
  alphas <- c("richness", "shannon", "invsimpson")
  
  class_atb_single <- merge(class_atb_single, mos[, c("lopnrMOS",full.model,"family", alphas, "last_atb_period"), with=F], by = "lopnrMOS")
  class_atb_single[, last_atb_period := factor(last_atb_period, c("0","1yr","1_4yr","4_8yr"), c("0","1_4yr","1_4yr","4_8yr"))]
  class_atb_single[last_atb_period!=0, class:=paste0(class,last_atb_period )]
  class_atb_single <- merge(class_atb_single, clr.species, by = "lopnrMOS")
  
  
    # Create the mixed-model regression function 
  
  linear.fun <- function(species, exposure, model = basic.model, df= class_atb_single){
    
    cc <- complete.cases(df[, c(model, exposure), with=F])
    tab <- df[cc, .N, by = "class"]
    tab$class <- paste0("class",tab$class)
    
    res <- bplapply(species, function(y){
      
      form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", exposure,"+(1|family)"))
      
      fit <- lmer(form, data = class_atb_single)
      exp_names = grep(exposure,rownames(summary(fit)$coef),value=T)
      suppressMessages(ci <- confint(fit, parm = exp_names))
      
      temp.res <- summary(fit)
      
      N = length(temp.res$residuals)
      res <- as.data.frame(temp.res$coef[exp_names, ])
      colnames(res) <- c("beta", "SE","df", "t.value","p.value")
      res$df <- NULL
      
      data.frame(outcome = y, exposure = exp_names, res, lci=ci[,1],hci=ci[,2], N = N)
      
    }, BPPARAM = MulticoreParam(16))
    
    res <- do.call(rbind, res)
    
    setDT(res)
    
    res <- merge(res, tab, by.x = "exposure", by.y = "class", all.x=T)
    res[,q.value := p.adjust(p.value, method = "BH")]
    res[,model := ""]
    
  }
  
  
  
  # Species -----------------------------------------------------------
  
  res <- linear.fun(prevalent.species , exposure = "class",  model = basic.model, df = class_atb_single)
  res$model <- "basic.model"
  res$cohort <- "MOS"
  
  res_full <- linear.fun(prevalent.species , exposure = "class",  model = full.model, df = class_atb_single)
  res_full$model <- "full.model"
  res_full$cohort <- "MOS"
  
  res <- rbind(res, res_full)
  
  fwrite(res, file='results/mos_species_singledose.tsv', sep = '\t')
  
  # Alphas --------------------------------------------------------------
  res <- linear.fun(alphas , exposure = "class",  model = basic.model, df = class_atb_single)
  res$model <- "basic.model"
  res$cohort <- "MOS"
  
  res_full <- linear.fun(alphas , exposure = "class",  model = full.model, df = class_atb_single)
  res_full$model <- "full.model"
  res_full$cohort <- "MOS"
  
  res <- rbind(res, res_full)
  
  fwrite(res, file='results/mos_alpha_singledose.tsv', sep = "\t")
  
  
  message("End")
  
  t1 = Sys.time()
  print(t1-t0)
  
  