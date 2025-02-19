# Project Antibiotic Use and the Gut Microbiota


# Load packages
  library(data.table)
  library(BiocParallel)
  library(car)
  

# This script will investigate the associations of single antibiotic course
  
  # Import data set
  simpler <- fread("work/simpler_working_dataset.csv", na.strings = c("NA", NA, ""))
  registerdata   <- fread('/proj/simp2023007/data/phenotypes/prescibed_drugs.csv', na.strings = c("NA",NA,"")) 
  clr.species <- fread('work/simpler_clr_species.tsv', na.strings = c("NA", NA, ""), data.table = F)
  
  
  # Import model and species by abundance
  load('work/simpler_model.Rdata')
  
  # Filter register data
  registerdata <- registerdata[grep("^J01", ATC),]
  registerdata <- merge(simpler[, .(SIMPKEY, Visit1)], registerdata, by="SIMPKEY")

  # Atb before 
  registerdata_N1 <- registerdata[SIMPKEY %in% simpler$SIMPKEY[simpler$N_all==1], ]
  registerdata_N1 <- registerdata_N1[EDATUM<Visit1,]
  registerdata_N1 <- registerdata_N1[which(Visit1-EDATUM<8*365.2),]
  all(registerdata_N1[,.N, by = SIMPKEY][,N]==1)
  
  class_atb_single <- registerdata_N1[ , .(SIMPKEY, ATC)]
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
  
  no_atb <- simpler[N_all==0, .(SIMPKEY, class = last_atb_period)]
  class_atb_single <- rbind(class_atb_single[, -"ATC"], no_atb)

  
  class_atb_single <- merge(class_atb_single, simpler[, c("SIMPKEY",full.model, alphas, "last_atb_period"), with=F], by = "SIMPKEY")
  class_atb_single[, last_atb_period := factor(last_atb_period, c("0","1yr","1_4yr","4_8yr"), c("0","1_4yr","1_4yr","4_8yr"))]
  class_atb_single[last_atb_period!=0, class:=paste0(class,last_atb_period )]
  class_atb_single <- merge(class_atb_single, clr.species, by = "SIMPKEY")

# Create the linear regression function 

linear.fun <- function(species, exposure, model = basic.model, df = class_atb_single){
  
  cc <- complete.cases(df[, c(model, exposure), with=F])
  tab <- df[cc, .N, by = "class"]
  tab$class <- paste0("class",tab$class)
  
  res <- bplapply(species, function(y){
    
          form <- as.formula(paste0(y, "~ ", paste0(model, collapse="+"),"+", paste0(exposure, collapse="+")))
          
          fit <- lm(form, data = df)
          exp_names <- grep(exposure, rownames(summary(fit)$coef), value = T)
          ci <- confint(fit, parm = exp_names)
    
          temp.res <- summary(fit)
    
          N = length(temp.res$residuals)
          res <- temp.res$coef
          res <- res[grep(exposure,rownames(res),value=T), ]
          colnames(res) <- c("beta", "SE", "t.value","p.value")
    
           data.frame(outcome = y, exposure = rownames(res), res, N = N, lci = ci[,1], hci = ci[,2])
  
    }, BPPARAM = MulticoreParam(16))
  
  res <- do.call(rbind, res)
  
  setDT(res)
  
  res <- merge(res, tab, by.x = "exposure",by.y = "class", all.x=T)
  res[,q.value := p.adjust(p.value, method = "BH")]
  res[,model := ""]
  
  
}


  # alphas ####------------------------------------------------------
  alphas <- c("richness", "shannon", "invsimpson")
  
  
  res <- linear.fun(alphas, exposure = "class" , model = basic.model, df = class_atb_single)
  res$model <- "basic.model"
  res$cohort <- "SIMPLER"
  
  
  res_dis <- linear.fun(alphas, exposure = "class" , model = full.model, df = class_atb_single)
  res_dis$model <-  "full.model"
  res_dis$cohort <- "SIMPLER"
  
  res <- rbind(res, res_dis)
  
  fwrite(res, file='results/simpler_alpha_singledose_class.tsv')
  
  
  message("End ")

  # Species ####------------------------------------------------------

  res <- linear.fun(prevalent.species, exposure = "class" , model = basic.model, df = class_atb_single)
  res$model <- "basic.model"
  res$cohort <- "SIMPLER"
  
  
  res_dis <- linear.fun(prevalent.species , exposure = "class" , model = full.model, df = class_atb_single)
  res_dis$model <-  "full.model"
  res_dis$cohort <- "SIMPLER"
  
  res <- rbind(res, res_dis)
  
  fwrite(res, file='results/simpler_singledose_class.tsv')
  

   
