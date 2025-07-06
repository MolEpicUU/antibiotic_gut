# Project Antibiotic Use and the Gut Microbiota

rm(list=ls())

  # This script will process the MOS phenotype data
  
  # Load packages
  library(data.table)
  library(dplyr)
  library(tidyr)


  simpler <- fread("nobackup/simp2023007/users/baldanzi/work/simpler_noexclusions.tsv", , na.strings = c("NA",NA,""))
  drugs   <- fread('data/phenotypes/prescibed_drugs.csv', na.strings = c("NA",NA,""))
  microb  <- fread('data/microbiome/processed/simpler_metagenomics_mgs_relative_abundances_v2.0.tsv') 
  downmgs <- fread('data/microbiome/processed/simpler_metagenomics_mgs_ds_relative_abundances_v2.0.tsv', data.table=F, na.strings = c("NA",NA,""))
  alpha   <- fread('data/microbiome/processed/simpler_metagenomics_alpha_diversity_v2.0.tsv', na.strings=c("NA","",NA)) 
  
  
  # Restrict to participants with gut microbiota data 
  simpler <- simpler[SIMPKEY %in% microb$SIMPKEY, ] # 5888 # It will increase when using microbiota_v1.1
 
  # Ensure that all participants have a minimum follow up of 8 years. 
  simpler[ as.Date(Visit1)  <= as.Date("2013-07-01"), .N]
  simpler <- simpler[ as.Date(Visit1)  > as.Date("2013-07-01"), ] # 5442
  
  # Exclusion based on the phenotype data 
  simpler <- simpler[ COPD == "no", ]
  simpler <- simpler[ IBD  == "no", ]  
  
  # Antibiotic drug management --------------------------------------------------
  atb <-  copy(drugs[grep("^J01", ATC), ])
  length(unique(atb$SIMPKEY)) 
  
  # Restricting to participants included in the phenotype data. 
  atb <- atb[SIMPKEY %in% simpler$SIMPKEY,] 
  length(unique(atb$SIMPKEY)) 
  
  atb[grep("^J0"   ,ATC), class:="other"]
  atb[grep("^J01CA",ATC), class:="Class_Peni_Ext"]
  atb[grep("^J01CE",ATC), class:="Class_Peni_BetaS"]
  atb[grep("^J01CF",ATC), class:="Class_Peni_BetaR"]
  atb[grep("^J01CR",ATC), class:="Class_Peni_Comb"]
  atb[grep("^J01D[B,C,D,E]",ATC), class:="Class_cephalosporins"]
  atb[grep("^J01FA",ATC), class:="Class_macrolides"] 
  atb[grep("^J01FF",ATC), class:="Class_lincosamides"]
  atb[grep("^J01MA",ATC), class:="Class_FQs"]
  atb[grep("^J01A" ,ATC),  class:="Class_TCLs"]
  atb[grep("^J01E" ,ATC),  class:="Class_SMZTMP"]
  atb[grep("^J01XE" ,ATC), class:="Class_NIT"]
  
  
  atb <- merge(atb, simpler[,c("SIMPKEY","Visit1")], by = "SIMPKEY",  all=T)
  
  # Create object with the prescription after Visit1
  atb_after_microb <- copy(atb)
  atb_after <- atb[which( as.Date(EDATUM) > (Visit1) ), ]

  
  # Ignore prescriptions that occurred after fecal collection 
  atb[,AtbDay_Day1 := as.Date(Visit1) - as.Date(EDATUM)] # Calculate intervals between prescription and visit 1
  atb[AtbDay_Day1<=0, AtbDay_Day1 := NA]
  atb[EDATUM>=(Visit1), ATC := NA ]
  atb[EDATUM>=(Visit1), EDATUM := NA ]
  atb[EDATUM>=(Visit1), class := NA ]
  atb[EDATUM>=(Visit1), AtbDay_Day1 := NA] # 33278

  
  # Remove long-term antibiotic users during fecal sampling ####
  message("Long term users")

  doxy <- atb[ ATC %in% c("J01AA02") & AtbDay_Day1<=90 , .(doxy=sum(forpddd),ATC,EDATUM,forpddd,AtbDay_Day1), by=SIMPKEY]
  doxy.users <- doxy[ (AtbDay_Day1<=42 & doxy>=42) | (AtbDay_Day1<=56 & doxy>=56) | (AtbDay_Day1<=84 & doxy>=84)  , unique(SIMPKEY)]
  doxy40.users <- atb[ grep("40", lnmn), ] %>% filter(ATC %in% c("J01AA02") & AtbDay_Day1<=90)  %>% pull(SIMPKEY)
  
  tetra <- atb[ ATC %in% c("J01AA07") & AtbDay_Day1<=90 , .(tetra=sum(forpddd),ATC,EDATUM,forpddd,AtbDay_Day1), by=SIMPKEY]
  tetracycline.users <- tetra[ (AtbDay_Day1<=42 & tetra>=21) | (AtbDay_Day1<=56 & tetra>=28) | (AtbDay_Day1<=84 & tetra>=42), unique(SIMPKEY)]
  
  lyme <- atb[ ATC %in% c("J01AA04") & AtbDay_Day1<=90 , .(lyme=sum(forpddd),ATC,EDATUM,forpddd,AtbDay_Day1), by=SIMPKEY]
  lymecycline.users <- lyme[ (AtbDay_Day1<=42 & lyme>=21) | (AtbDay_Day1<=56 & lyme>=28) | (AtbDay_Day1<=84 & lyme>=42), unique(SIMPKEY)]
  
  methenamine.users <- atb[AtbDay_Day1<=90 & ATC == "J01XX05", unique(SIMPKEY)] 
  
  nitro_trime.users <- atb[AtbDay_Day1<=90 & ATC %in% c("J01XE01","J01EA01") & forpddd>=22.5, unique(SIMPKEY)]
  
  longterm.users <- c(doxy.users, tetracycline.users, lymecycline.users, methenamine.users, nitro_trime.users) # 27
  
  atb <- atb[-which(SIMPKEY %in% longterm.users),] 
  simpler <- simpler[-which(SIMPKEY %in% longterm.users),] 
  
  message(paste(length(longterm.users), "long term antibiotic users during fecal sampling")) # 27  nrow(atb) 32052
  
  
  # Restrict the atb data to max 8 years before visit
  atb[AtbDay_Day1>8*365.2,ATC:=NA] 
  atb[AtbDay_Day1>8*365.2,EDATUM:=NA]
  atb[is.na(ATC), AtbDay_Day1 := NA]
  
  # Final sample size calculations
  atb <- merge(simpler[,"SIMPKEY"], atb, by = "SIMPKEY", all=T)
  n_final <- length(unique(atb$SIMPKEY)) 
  atb[,atb:=1]
  atb[is.na(ATC),atb:=0]
  temp_atb <- atb[,.(N_atb=sum(atb)),by=SIMPKEY]
  n_noatb <- nrow(temp_atb[N_atb==0,]) 
  
  #### Number of participants left ####
  message(paste0("Number of participants left = ",length(unique(atb$SIMPKEY))))
  message(paste0("Number of participants that did not have any antibiotic = ",n_noatb))
  
  
  # Categorize the prescriptions according to the period when it was dispensed ####
  atb[ ,t1yr :=  0 ]
  atb[ ,t1_4yr :=  0 ]
  atb[ ,t4_8yr :=  0 ]
  
  atb[AtbDay_Day1<365, t1yr := 1 ]
  atb[AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2), t1_4yr := 1 ]
  atb[AtbDay_Day1>=(4*365.2), t4_8yr := 1 ]
  
  
  # Add the number for atb courses per participant per period
  N_atb <- atb[,.(N1yr = sum(t1yr) , N1_4yr=sum(t1_4yr) , N4_8yr = sum(t4_8yr), N_all = sum(t1yr,t1_4yr,t4_8yr)) , by = SIMPKEY]
  
  
  classes <- unique(atb$class[!is.na(atb$class)])
  
  for(c in classes){
    
    atb[!is.na(ATC) ,t1yr   :=  0 ]
    atb[!is.na(ATC) ,t1_4yr :=  0 ]
    atb[!is.na(ATC) ,t4_8yr :=  0 ]
    
    atb[class==c & AtbDay_Day1<365, t1yr := 1 ]
    atb[class==c & AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2), t1_4yr := 1 ]
    atb[class==c & AtbDay_Day1>=(4*365.2), t4_8yr := 1 ]
    
    
    # Add the number for atb courses per participant per period
    N_atb_class <- atb[,.(N1yr=sum(t1yr) , N1_4yr=sum(t1_4yr) , N4_8yr = sum(t4_8yr), N_all = sum(t1yr,t1_4yr,t4_8yr)) , by = SIMPKEY]
    
    names(N_atb_class) <- c("SIMPKEY", paste0(c, c("_1yr","_1_4yr","_4_8yr", "_all"))) 
    
    if(all(N_atb$SIMPKEY == N_atb_class$SIMPKEY)){N_atb <- cbind(N_atb,N_atb_class[,-"SIMPKEY"])} else {
      N_atb <- merge(N_atb,N_atb_class, by = "SIMPKEY")
    }
    
  }
  
  
  # Most recent antibiotic per participant
  
  atb_lastdose <- atb[,.SD[which.min(AtbDay_Day1)], by=.(SIMPKEY)]
  atb_lastdose <- atb_lastdose[,.(SIMPKEY,EDATUM,ATC,class,AtbDay_Day1)]
  
  # Categorize period of the most recent antibiotic 
  atb_lastdose[AtbDay_Day1<365 , last_atb_period:="1yr"]
  atb_lastdose[AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2) ,last_atb_period:="1_4yr"]
  atb_lastdose[AtbDay_Day1>=(4*365.2) ,last_atb_period:="4_8yr"]
  
  # Rename variables for most recent antibiotic
  
  a <- c("EDATUM","ATC","class","AtbDay_Day1")
  setnames(atb_lastdose, a, paste0("last_",a))
  
  # Merge all antibiotic data (now in wide format)
  N_atb <- merge(N_atb,atb_lastdose, by="SIMPKEY", all.x = T )
  N_atb$last_atb_period[is.na(N_atb$last_ATC)] <- "0"
  
  
  # Single dose ----------------------------------------------------------------
  singledoseSIMPKEY <-  N_atb[N_all==1,"SIMPKEY"] %>% pull(SIMPKEY)
  singledose <- atb[SIMPKEY %in% singledoseSIMPKEY & !is.na(ATC),.(SIMPKEY,AtbDay_Day1)]
  setnames(singledose, "AtbDay_Day1","singledoseAtbDay_Day1")
  
  simpler[,singledose := 0]
  simpler[SIMPKEY %in% singledoseSIMPKEY ,singledose := 1]
  simpler <- merge(simpler, singledose, by = "SIMPKEY", all=T)
  
  # Merge all antibiotic data with the simpler 
  simpler <- merge(simpler, N_atb, by="SIMPKEY") # 4609
  
  
  # Alpha diversity -----------------------------------------------------------
  richness <- data.frame(SIMPKEY = downmgs$SIMPKEY, richness = rowSums(downmgs[,-1] > 0))
  
  simpler <- merge(simpler, richness, by="SIMPKEY")
  simpler <- merge(simpler, alpha , by="SIMPKEY")
  
  # Save final simpler data ####
  fwrite(simpler, file = "users/baldanzi/work/simpler_working_dataset_revision_noatbexclusion.tsv", sep="\t")
  
  
  simpler[ last_EDATUM >= (Visit1-30), .N]
  simpler[is.na(Visit1), .N]
  simpler <- simpler[is.na(last_EDATUM) | last_EDATUM < (Visit1-30), ]
  fwrite(simpler, file = "users/baldanzi/work/simpler_working_dataset_revision.tsv", sep="\t")
  
  # Calculate species prevalence -------------------------------------------------------------------------------------------------
  # Main model covariates 
  basic.model <- c("age", "sex", "placebirth", "smokestatus", "education", "center","dna_extraction_plate")
  #full.model <- c(main.model, "BMI","CCIw", "diabd","rheumatic","cancer", "ppi")
  
  full.model <- c(basic.model, "BMI", "CCIw", "ppi", "statins", "metformin", "betablock", "ssri", "polypharmacy12m_cat")
  
  cc <- complete.cases(simpler[, basic.model, with=F])
  
  species <- ifelse(as.matrix(downmgs[downmgs$SIMPKEY %in% simpler[cc,SIMPKEY], -1]) == 0, 0, 1)
  prevalence <- colSums(species)/nrow(species)
  
  prevalent.species <- names(prevalence[prevalence >.02]) 
  non.prevalent.species <- names(prevalence[prevalence <=.02])
  
  
  save(list=c("basic.model","prevalent.species","non.prevalent.species","full.model") , 
       file='/proj/nobackup/simp2023007/users/baldanzi/work/simpler_model_revision.Rdata')
  
  message(paste0("Basic model saved = ",paste0(basic.model,collapse = ", ")))
  message(paste0("Full model saved = ",paste0(full.model,collapse = ", ")))
  message(paste0("Length of prevalent species = ",length(prevalent.species)))
  
  
  message("End")
  
  
  