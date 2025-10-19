# Project antibiotic use and gut microbiota 

# This script will process the MOS phenotype data, merge with the MGS 2.0 data
# and process  the drug regsister data. 

# last update: 18 April 2025

  rm(list=ls())
  library(data.table)
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  
  
  # Import data ####
  
  mos1  <-  rio::import('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/data_MOS.xlsx') # Phenotypes for MOS participants
  mos2  <-  rio::import('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/data39only_MOS.xlsx') # Phtenoypes for MOS participants that did not have 16S data
  mos3  <-  rio::import('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/MOS Half-time n2644 Questionnaire F29_TO_F42.sav') # Additional phenotypes 
  dates <-  rio::import('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/MOS1_datum_2023-01-24_Pawel_Marju.xls', range='A1:D2500')
  atb   <-  fread('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/MOS_Halftime_n2644_Prescribed_Drug_register_until_2019.tsv')
  family.id     <-  fread('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/MOSlopnr_family_n2644.csv')
  overlapscapis <-  fread('/proj/sens2019512/MOS/Phenotypes/overlap_scapis_mods_ids.csv')
  microb  <-  fread('/proj/sens2019512/MOS/Metagenomics/mos_metagenomics_mgs_relative_abundances_v1.0.tsv')
  cm      <-  fread('/proj/sens2019512/MOS/Metagenomics/mos_metagenomics_technical_variables_v1.0.tsv')
  cci     <- fread("/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_charlsonindex.tsv", na.strings=c("", NA, "NA"))
  polypharmacy <- rio::import('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/Polypharmacy_MOS_n1748_to_GB_AL250417.sav')
  med_28 <- rio::import('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/Medications_MOS_variables_n28_to_GB_250511AL.sav')
  antipsycho <- rio::import('/proj/sens2019512/nobackup/users/baldanzi/atb_gut/Data/MOS/n05a_n1748_selfreported_MOS_toGB250818_AL.sav')
  
  a     <-  grep("HG3A.", colnames(mos1)) # Remove MGS1.0 species
  b     <-  grep("HG3A.", colnames(mos2)) # Remove MGS1.0 species 
  mos1  <-  mos1[,-a]
  mos2  <-  mos2[,-b]
  mos2$BMI_3 <- mos2$BMI_4 <- NULL # Duplicated columns 
  mos1  <-  rename(mos1, physical_activity = ltpa) # Rename the physical activity variable 
  
  # Medications for those extra 39 MOS participants 
  mos2 <- merge(mos2, med_28, by = "lopnrMOS", all.x=T)
  
  
  mos2$plate <- NULL
  mos2$metformin <- NULL
  mos2$Energy_kcal <- NULL
  mos2$Fiber <- NULL
  
  # PPI ####
  mos1$ppi <- mos1$PPI_A02BC
  mos2$ppi <- mos2$PPI_A02BC
  mos2$PPI <- NULL
  mos2$PPI_A02BC <- NULL
  
  
  mos1  <- mos1[,c(colnames(mos2)) ]
  
  setDT(mos1)
  setDT(mos2)
  mos   <- rbind(mos1, mos2, fill=T)
  
  # mos[, polypharmacy :=  rowSums(.SD, na.rm=T), .SDcols = poli]

  
  # Medication names 
  setnames(mos, c("Statins_C10AA", "Metformin_ALL", "Beta_blockers_selective_C07AB", "Antidepressants_SSRI_N06AB"),
           c("statins", "metformin", "betablock", "ssri"))

  rm(mos1)
  rm(mos2)
  rm(a)
  rm(b)
  
  mos <- rename(mos, Data_sample_code = "Data sample code")  # 2262
  
  setDT(mos)
  
 # Data management 
  mos[, sex := factor(sex, c(1,2), c("male","female") )]
  mos[, bornSweden := factor(country_birth, c(0,1), c("no","yes") )]
  mos[, smoking := factor(smoking, 1:4, c("never","current","current","former"))]
  # mos[!is.na(glucose) & diab_diag==0 & glucose >=6.1, diab_diag := 0.5]
  # mos[!is.na(HbA1c) & diab_diag==0 & HbA1c >=42, diab_diag := 0.5]
  mos[!is.na(glucose) & diab_diag ==0 & glucose >7, diab_diag := 1]
  mos[!is.na(HbA1c) & diab_diag == 0 & HbA1c >=48, diab_diag := 1]
  mos[, diab_diag := factor(diab_diag, c(0,1), c("no","yes"))]
  mos[, BP_diag := factor(BP_diag, c(0,1), c("no","yes") )]
  mos[, education:=factor(education, 1:4, c("Compulsory", "Compulsory", "Upper secondary", "University"))]
  
  
  # Visit and fecal sample collection dates -------------------------------------------------------------
 dates <- dates[-which(is.na(dates$lopnr_mos)),]
 setDT(dates)
 dates[V1_date_antro == '2014-05-14' & V2_date_sphyg == '2014-05-22' & faeces_datum == '2015-05-22', faeces_datum:= '2014-05-22' ]
 mos3$faeces_datum <- NULL
 
 mos3 <- merge(mos3, dates[,.(faeces_datum, lopnr_mos)], by.x="lopnrMOS_str9", by.y="lopnr_mos", all.x=T, all.y=F)
 
 
 # Data management mos dates  ----------------------------------------------------------------------------------
  
  # Rewrite Visit 1 and Visit 2 since some participants had V2 before V1 
 setDT(mos3) 
  mos3[, Visit1 := date_antro]
  mos3[, Visit2 := date_sphyg]
  mos3[ date_antro > date_sphyg, Visit1 := date_sphyg]
  mos3[ date_antro > date_sphyg, Visit2 := date_antro]
  mos3[is.na(date_sphyg), Visit2 := faeces_datum]
  mos3[is.na(date_sphyg) & is.na(faeces_datum), Visit2 := Visit1]
  mos3[, Visit1_Visit2 := as.Date(Visit2) - as.Date(Visit1)]
  
  # Cancer 
  mos3[,cancer:="no"]
  mos3[, YearV1 := as.numeric(format(mos3$Visit1, format="%Y"))]
  mos3[, YearCancer := as.numeric(F29_1)]
  mos3[F29 == 1 & (YearCancer==YearV1 | (YearCancer+1)==YearV1), cancer:="1 year before V1"]
  mos3[F29 == 1 & (YearV1-YearCancer)>1 & (YearV1-YearCancer)<=5, cancer:="2-5 years before V1"]
  mos3[F29 == 1 & (YearV1-YearCancer)>6 & (YearV1-YearCancer)<=10, cancer:="6-10 years before V1"]
  
  mos3[, rheumatic := "no"]
  mos3[F30==1, rheumatic := "yes"]
  
  mos3[, IBD := "no"]
  mos3[F37==1 | F38==1, IBD := "yes"]
  
  mos3[, LungDisease := "no"]
  mos3[F40==1 | F41==1 | F42==1, LungDisease := "yes"]
  
  
  mos$sex <- NULL
  stopifnot(all(mos$lopnrMOS %in% mos3$lopnrMOS))
  mos <- merge(mos, mos3, by="lopnrMOS" , all.x=T, all.y=F) 
  
  # Clinical Microbiomics technical data ####
  mos <- merge(mos, cm[,c("lopnrMOS","extraction_plate")] , by = "lopnrMOS", all.x = T, all.y =F)
  
  # Family ID ####
  family.id$V1 <- NULL
  mos <- merge(family.id, mos , by="lopnrMOS", all.y=T)
  # Transform family into a character variable 
  mos$family <- paste0("f",mos$family)
  
  # Coverage / Follow-up ####
  mos[,followup := as.Date(mos$Visit1) - as.Date("2005-07-01")]
  
  
  # Charlson comorbidity index weighted #### -------------------------------------------------------------
  cci <- merge(mos[, .(lopnrMOS_str9, diab_diag, age)], cci, by.x="lopnrMOS_str9", by.y="group", all.x=T )
  
  cci[is.na(cci)] <- 0
  
  # The patient records do not capture uncomplicated diabetes well, so we will enrich that information using self-reported
  cci[diab_diag == 1 & 
        Diabetes_without_chronic_complication == 0 &
        Diabetes_with_chronic_complication == 0, 
      Diabetes_without_chronic_complication := 1]
  

  # # Age was also not part of the inicial CCIw
  # cci[, Age := fcase(
  #   age < 50, 0,
  #   age >= 50 & age < 60, 1,
  #   age >= 60 & age < 70, 2,
  #   age >= 70 & age < 80, 3,
  #   age >= 80, 4
  # )]
  
  # Re-calculate the weighted comorbidity index
  cci$CCIw <- cci$Myocardial_infarction + cci$Congestive_heart_failure + cci$Peripheral_vascular_disease + 
    cci$Cerebrovascular_disease + cci$Chronic_obstructive_pulmonary_disease + cci$Chronic_other_pulmonary_disease + 
    cci$Rheumatic_disease + cci$Dementia + 2*cci$Hemiplegia + cci$Diabetes_without_chronic_complication + 
    2*cci$Diabetes_with_chronic_complication + 2*cci$Renal_disease + cci$Mild_liver_disease + 3*cci$Severe_liver_disease + 
    cci$Peptic_ulcer_disease + 2*cci$Malignancy + 6*cci$Metastatic_solid_tumor + 6*cci$Aids # + cci$Age
  
  
  cci[is.na(CCIw), CCIw := 0]
  mos <- merge(mos, cci[, c("lopnrMOS_str9", "CCIw")], by="lopnrMOS_str9",  all.x=T )
  
  # Polypharmacy ####
  setnames(polypharmacy, c("N_selfrep_Meds_Non_antibiotics_compute_variable"),  c("polypharmacy"))
  
  mos <- merge(mos, polypharmacy, by="lopnrMOS", all.x=T)
  mos[polypharmacy == 0 , polypharmacy_cat := "0"]
  mos[polypharmacy > 0 , polypharmacy_cat := "1-4"]
  mos[polypharmacy > 4 , polypharmacy_cat := "5 or more"]
  
  # Antipsychotics ####
  setnames(antipsycho, "N05A_Antipsychotics_new_0_1", "antipsycho")
  mos <- merge(mos, antipsycho, by = "lopnrMOS", all.x=T)
  mos[, antipsycho := factor(antipsycho, c(0,1), c("no", "yes"))]
  
  
  # Save a temporary mos file without exclusions ####
  fwrite(mos, "/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_noexclusions.tsv")
  # mos <- fread("/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_noexclusions.tsv")
  
  
  ### Exclusion of participants based on the phenotype data ###
  
  
  # Remove participants that were also in SCAPIS 
  mos <- mos[!lopnrMOS_str9 %in% overlapscapis$mosnr, ]
  
  
  # visit1, visit2, faeces_datum DATES ####
  # mos[as.Date(faeces_datum)<= as.Date(Visit2),dd := "1_beforeV2"]
  # mos[as.Date(faeces_datum)<= as.Date(Visit2)+7 & as.Date(faeces_datum)> as.Date(Visit2), dd := "2_beforeV2plus7"]
  # mos[as.Date(faeces_datum)> as.Date(Visit2)+7 ,dd := "3_afterV2plus7"]
  # mos[is.na(faeces_datum) ,dd := "4_NA"]
  
  
  # Restrict to participants with gut microbiota data 
  stopifnot(all(mos$lopnrMOS %in% microb$lopnrMOS))
  
  
  temp_mos <- mos
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
  
  message(paste0("Number of participants left = ",length(unique(atb$lopnrMOS)))) # 2018
  message(paste0("Number of participants that did not have any antibiotic = ",n_noatb)) # 563
  
  # Save the processed prescribed drug register data 
  
  fwrite(atb,"/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_processed_lakemedelregister_revision.tsv")
  # atb = fread("/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_processed_lakemedelregister.tsv")
  
  
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
  downmgs <- fread("/proj/sens2019512/MOS/Metagenomics/mos_metagenomics_mgs_ds_relative_abundances_v1.0.tsv", data.table=F)
  alpha <- fread('/proj/sens2019512/MOS/Metagenomics/mos_metagenomics_alpha_diversity_v1.0.tsv')
  
  richness <- data.frame(lopnrMOS = downmgs$lopnrMOS, richness =  rowSums(downmgs[,-1] > 0))
  
  mos <- merge(mos, alpha, by="lopnrMOS")
  mos <- merge(mos, richness, by="lopnrMOS")
  
  # Save final phenotype MOS data ####
  fwrite(mos, file = "/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_working_dataset_revision_noatbexclusion.tsv")
  
  mos[last_EDATUM >= (Visit1-30), .N]
  mos[is.na(Visit1), .N]
  mos <- mos[is.na(last_EDATUM) | last_EDATUM < (Visit1-30), ]
  fwrite(mos, file = "/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_working_dataset_revision.tsv")
  
  
  ## Calculate species prevalence ----------------------------------------------
  basic.model <- c("age","sex","country_birth","smoking","education","extraction_plate")
  full.model <- c(basic.model,"BMI","CCIw", "ppi", "statins", "metformin", "betablock", "ssri", "polypharmacy_cat", "antipsycho")
  
  cc <- complete.cases(mos[, basic.model, with=F])
  
  
  species <- ifelse(as.matrix(downmgs[downmgs$lopnrMOS %in% mos[cc, lopnrMOS],-1]) == 0, 0, 1)
  prevalence <- colSums(species)/nrow(species)
  
  prevalent.species <- names(prevalence[prevalence >.02]) # 1368
  non.prevalent.species <- names(prevalence[prevalence <=.02])
  
  
  save(list=c("basic.model","prevalent.species","non.prevalent.species", "full.model"),
       file='/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/mos_model_revision.Rdata')
  
  message(paste0("Full model saved = ",paste0(full.model,collapse = ", ")))
  message(paste0("Length of prevalent species = ",length(prevalent.species)))
  
  # Prevalence table 
  prev_table <- data.table(species = names(prevalence), prevalence_mos = prevalence)
  
  fwrite(prev_table, '/proj/sens2019512/nobackup/users/baldanzi/atb_gut/results/prevalence_species_mos.tsv', sep = '\t')
  
  
  
  message("End")
  