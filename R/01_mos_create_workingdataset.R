# Project antibiotic use and gut microbiome

# This script will process the MOS phenotype data

  rm(list=ls())

  library(data.table)
  library(dplyr)
  library(tidyr)
  
  
  # Import data ####
  atb   <-  fread('nobackup/users/baldanzi/atb_gut/Data/MOS/MOS_Halftime_n2644_Prescribed_Drug_register_until_2019.tsv')
  overlapscapis <-  fread('MOS/Phenotypes/overlap_scapis_mods_ids.csv')
  microb  <-  fread('MOS/Metagenomics/mos_metagenomics_mgs_relative_abundances_v1.0.tsv')
  mos <- fread("nobackup/users/baldanzi/atb_gut/work/mos_noexclusions.tsv", na.strings=c("NA", NA, ""))
  
  
  ### Exclusion of participants based on the phenotype data ###
  
  # Remove participants that were also in SCAPIS 
  mos <- mos[-which(lopnrMOS_str9 %in% overlapscapis$mosnr) , ]
  
  # Ensure that all participants have a minimum follow up of 8 years. 
  mos <- mos[as.Date(Visit1)  > as.Date("2013-07-01"), ] # 97 are excluded 
  
  # Remove participants with interval between visit1 and 2 greater than 60 days. 
  mos <- mos[-(which(Visit1_Visit2>60)),] # 31 are excluded 
  
  # Exclusion of participants with fecal samples collected more than 7 days after visit 2
  mos <- mos[ -which(as.Date(faeces_datum) > (as.Date(Visit2)+7)),] # 9 are excluded 
  
  # Exclusion of chronic pulmonary lung disease or IBD 
  mos <- mos[-which(LungDisease == "yes"),] # 39 are excluded 
  mos <- mos[-which(IBD == "yes"),] # 23 are excluded 
  
  # -----------------------------------------------------------------------------------------------------------------------
  
  # Drug Register Data ####
  # This dataset is in long format, meaning that every row is a prescription. Participants are represented in multiple lines

  # Restricting to participants included in the phenotype data. 
  atb <- atb[lopnrMOS %in% mos$lopnrMOS,] 
  
  # Antibiotic classes ####
  atb[grep("^J0",ATC), class:="other"]
  atb[grep("^J01CA",ATC), class:="Class_Peni_Ext"]
  atb[grep("^J01CE",ATC), class:="Class_Peni_BetaS"]
  atb[grep("^J01CF",ATC), class:="Class_Peni_BetaR"]
  atb[grep("^J01CR",ATC), class:="Class_Peni_Comb"]
  atb[grep("^J01D[B,C,D,E]",ATC), class:="Class_cephalosporins"]
  atb[grep("^J01FA",ATC), class:="Class_macrolides"] 
  atb[grep("^J01FF",ATC), class:="Class_lincosamides"]
  atb[grep("^J01MA",ATC), class:="Class_FQs"]
  atb[grep("^J01A",ATC), class:="Class_TCLs"]
  atb[grep("^J01E",ATC), class:="Class_SMZTMP"]
  atb[grep("^J01XE",ATC), class:="Class_NIT"]
  
  
  # Merge antibiotic data with dates and intervals from phenotype
  atb <- merge(atb, mos[,.(lopnrMOS, Visit1, Visit2, faeces_datum)], by='lopnrMOS')
  setnames(atb, "edatum", "EDATUM")
  atb[, c("Visit1", "Visit2","faeces_datum","EDATUM") := lapply(.SD,as.Date), .SDcols = c("Visit1", "Visit2","faeces_datum","EDATUM") ] 
  
  # Create object with the prescription after 7 days of visit 2 
  atb_after <- atb[which( EDATUM>=(Visit2+7) ), ]
  
  # Ignore prescriptions that occurred after 7 days of visit 2 
  atb[EDATUM>=(Visit2+7), ATC := NA ]
  atb[EDATUM>=(Visit2+7), EDATUM := NA ]
  atb[EDATUM>=(Visit2+7), class := NA ]
  
  # Remove participants with antibiotic prescription between visits 
  ids.to.remove <- atb[EDATUM<=(Visit2+7) & EDATUM >= Visit1, unique(lopnrMOS)]
  
  atb <- atb[-which(lopnrMOS %in% ids.to.remove),]
  mos <- mos[-which(lopnrMOS %in% ids.to.remove),] 
  
  # Remove long-term antibiotic users during fecal sampling ####
  atb[,AtbDay_Day1 := Visit1 - as.Date(EDATUM)] # Calculate intervals between prescription and visit 1
  atb[AtbDay_Day1<0, AtbDay_Day1 := NA]
  
  atb[, forpdddtotal :=  antal * forpddd]
  
  doxy <- atb[ ATC %in% c("J01AA02") & AtbDay_Day1<=90 , .(doxy=sum(forpdddtotal),ATC,EDATUM,forpddd,AtbDay_Day1), by=lopnrMOS]
  doxy.users <- doxy[ (AtbDay_Day1<=42 & doxy>=42) | (AtbDay_Day1<=56 & doxy>=56) | (AtbDay_Day1<=84 & doxy>=84)  , unique(lopnrMOS)]
  
  tetra <- atb[ ATC %in% c("J01AA07") & AtbDay_Day1<=90 , .(tetra=sum(forpdddtotal),ATC,EDATUM,forpddd,AtbDay_Day1), by=lopnrMOS]
  tetracycline.users <- tetra[ (AtbDay_Day1<=42 & tetra>=21) | (AtbDay_Day1<=56 & tetra>=28) | (AtbDay_Day1<=84 & tetra>=42), unique(lopnrMOS)]
  
  lyme <- atb[ ATC %in% c("J01AA04") & AtbDay_Day1<=90 , .(lyme=sum(forpdddtotal),ATC,EDATUM,forpddd,AtbDay_Day1), by=lopnrMOS]
  lymecycline.users <- lyme[ (AtbDay_Day1<=42 & lyme>=21) | (AtbDay_Day1<=56 & lyme>=28) | (AtbDay_Day1<=84 & lyme>=42), unique(lopnrMOS)]
  
  methenamine.users <- atb[AtbDay_Day1<=90 & ATC == "J01XX05", unique(lopnrMOS)] 
  
  nitro_trime.users <- atb[AtbDay_Day1<=90 & ATC %in% c("J01XE01","J01EA01") & forpdddtotal>=22.5, unique(lopnrMOS)]
  
  longterm.users <- c(doxy.users, tetracycline.users, lymecycline.users, methenamine.users, nitro_trime.users) # 28
  
  atb <- atb[-which(lopnrMOS %in% longterm.users),]
  mos <- mos[-which(lopnrMOS %in% longterm.users),]
  
  

  # Restrict the atb data to max 8 years before visit
  atb[AtbDay_Day1>8*365.2,ATC:=NA] 
  atb[AtbDay_Day1>8*365.2,EDATUM:=NA]
  atb[AtbDay_Day1>8*365.2,class:=NA]
  atb[is.na(ATC), AtbDay_Day1 := NA]
  
  
  # Final sample size calculations 
  atb <- merge(mos[,"lopnrMOS"], atb, by = "lopnrMOS", all.x=T, all.y=T)
  n_final <- length(unique(atb$lopnrMOS)) 
  atb[,atb:=1]
  atb[is.na(ATC),atb:=0]
  temp_atb <- atb[,.(N_atb=sum(atb)),by=lopnrMOS]
  n_noatb <- nrow(temp_atb[N_atb==0,]) 
  
  
  #### Number of participants left ####
  
  message(paste0("Number of participants left = ",length(unique(atb$lopnrMOS)))) 
  message(paste0("Number of participants that did not have any antibiotic = ",n_noatb)) 
  
  
  # Categorize the prescriptions according to the period when it was dispensed ####
  atb[ , t1yr :=  0 ]
  atb[ , t1_4yr :=  0 ]
  atb[ , t4_8yr :=  0 ]
  
  atb[AtbDay_Day1<365, t1yr := 1 ]
  atb[AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2), t1_4yr := 1 ]
  atb[AtbDay_Day1>=(4*365.2), t4_8yr := 1 ]
  
  
  # Add the number for atb courses per participant per period
  N_atb <- atb[,.(N1yr = sum(t1yr) , N1_4yr=sum(t1_4yr) , N4_8yr = sum(t4_8yr), N_all = sum(t1yr,t1_4yr,t4_8yr)) , by = lopnrMOS]
  
  
  # Create loop to  count the number of antibiotics for every class ####
  
  classes <- unique(atb$class[!is.na(atb$class)])
  
  for(c in classes){
    
    atb[ ,t1yr   :=  0 ]
    atb[ ,t1_4yr :=  0 ]
    atb[ ,t4_8yr :=  0 ]
    
    atb[class==c & AtbDay_Day1<365, t1yr := 1 ]
    atb[class==c & AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2), t1_4yr := 1 ]
    atb[class==c & AtbDay_Day1>=(4*365.2), t4_8yr := 1 ]
    
    
    # Add the number for atb courses per participant per period
    N_atb_class <- atb[,.(N1yr = sum(t1yr) , N1_4yr=sum(t1_4yr) , N4_8yr = sum(t4_8yr), N_all = sum(t1yr,t1_4yr,t4_8yr)) , by = lopnrMOS]
    
    names(N_atb_class) <- c("lopnrMOS", paste0(c, c("_1yr","_1_4yr","_4_8yr","_all"))) 
    
    if(all(N_atb$lopnrMOS == N_atb_class$lopnrMOS)){N_atb <- cbind(N_atb,N_atb_class[,-"lopnrMOS"])} else {
      N_atb <- merge(N_atb,N_atb_class, by = "lopnrMOS")
    }
    
  }
  
  # More recent antibiotic per participant
  atb_lastdose <- atb[,.SD[which.min(AtbDay_Day1)], by=.(lopnrMOS)]
  atb_lastdose <- atb_lastdose[,.(lopnrMOS,EDATUM,ATC,class,AtbDay_Day1)]
  
  # Categorize period of the most recent antibiotic 
  atb_lastdose[AtbDay_Day1<365 , last_atb_period:="1yr"]
  atb_lastdose[AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2) ,last_atb_period:="1_4yr"]
  atb_lastdose[AtbDay_Day1>=(4*365.2) ,last_atb_period:="4_8yr"]
  
  # Rename variables for most recent antibiotic
  
  a <- c("EDATUM","ATC","class","AtbDay_Day1")
  setnames(atb_lastdose, a, paste0("last_",a))
  
  # Merge all antibiotic data (now in wide format)
  N_atb <- merge(N_atb,atb_lastdose, by="lopnrMOS", all=T )
  N_atb$last_atb_period[is.na(N_atb$last_ATC)] <- "0"
  
  # Merge all antibiotic data with the phenotype data 
  mos <- merge(mos, N_atb, by="lopnrMOS")

  
  # Antibiotics after #### ----------------------------------------------------
  atb_after[ ,after_2yr :=  0 ]
  atb_after[, Day2_EDATUM := as.Date(EDATUM) - (Visit2+7)]
  atb_after[Day2_EDATUM<2*365.2, after_2yr := 1 ]
  
  
  # Add the number for atb courses per participant per period
  N_atb_after <- atb_after[,.(N_after_2yr = sum(after_2yr)) , by = lopnrMOS]
  
  mos <- merge(mos, N_atb_after, by="lopnrMOS", all.x=T)
  
  
  
  ## Alpha diversity metrics  ----------------------------------------------
  downmgs <- fread("MOS/Metagenomics/mos_metagenomics_mgs_ds_relative_abundances_v1.0.tsv", data.table=F)
  alpha <- fread('MOS/Metagenomics/mos_metagenomics_alpha_diversity_v1.0.tsv')
  
  richness <- data.frame(lopnrMOS = downmgs$lopnrMOS, richness =  rowSums(downmgs[,-1] > 0))
  
  mos <- merge(mos, alpha, by="lopnrMOS")
  mos <- merge(mos, richness, by="lopnrMOS")
  
  # Save final phenotype MOS data ####
  fwrite(mos, file = "nobackup/users/baldanzi/atb_gut/work/mos_working_dataset_revision_noatbexclusion.tsv")
  
  mos[last_EDATUM >= (Visit1-30), .N]
  mos[is.na(Visit1), .N]
  mos <- mos[is.na(last_EDATUM) | last_EDATUM < (Visit1-30), ]
  fwrite(mos, file = "nobackup/users/baldanzi/atb_gut/work/mos_working_dataset_revision.tsv")
  
  
  ## Calculate species prevalence ----------------------------------------------
  basic.model <- c("age","sex","country_birth","smoking","education","extraction_plate")
  full.model <- c(basic.model,"BMI","CCIw", "ppi", "statins", "metformin", "betablock", "ssri", "polypharmacy_cat", "antipsycho")
  
  cc <- complete.cases(mos[, basic.model, with=F])
  
  
  species <- ifelse(as.matrix(downmgs[downmgs$lopnrMOS %in% mos[cc, lopnrMOS],-1]) == 0, 0, 1)
  prevalence <- colSums(species)/nrow(species)
  
  prevalent.species <- names(prevalence[prevalence >.02]) # 1368
  non.prevalent.species <- names(prevalence[prevalence <=.02])
  
  
  save(list=c("basic.model","prevalent.species","non.prevalent.species", "full.model"),
       file='nobackup/users/baldanzi/atb_gut/work/mos_model_revision.Rdata')
  
  message("End")
  