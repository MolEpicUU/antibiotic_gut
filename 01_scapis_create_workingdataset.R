# Project Antibiotic Use and the Gut Microbiota

# Version date: 2024-02-13

  # Load packages
  library(data.table)
  library(dplyr)
  library(tidyr)
  
  
  # Import datasets 
  phenotype <- fread('/proj/nobackup/sens2019512/Projects/antibiotic_mgs/dataset/SCAPIS-DATA-PETITION-514-20230310.csv', na.strings = c("NA",NA,"","<NA>")) # nrow = 11249
  patreg    <- fread('/proj/sens2019512/SCAPIS/Gutsy/Phenotypes/Patient_register/Raw/SCAPIS-REGISTER-PETITION-514-20230310-UT_PAR_SV_9771_2018.txt', na.strings=c("NA","",NA)) # patient register 
  drugs     <- fread('/proj/sens2019512/SCAPIS/Gutsy/Phenotypes/Medication/Raw/SCAPIS-REGISTERSU-1715-LMED-20220511-T_R_LMED_27947_2021.txt', na.strings = c("NA",NA,""))
  microb    <- fread('/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_mgs_relative_abundances_v1.0.tsv')  # 9816 participants 
  alpha     <- fread("/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_shannon_diversity_v1.0.tsv",  data.table=F, na.strings = c("NA",NA,""))
  
  # Recoding variables ----------------------------------------------------------------------------
  
  setnames(scapis, 'AgeAtVisitOne', 'age')
  scapis[,Sex := factor(Sex, levels = c("MALE", "FEMALE"), labels = c("male", "female"))]
  setnames(scapis, 'derived_smoke_status', 'smokestatus')
  scapis[,smokestatus := factor(smokestatus, levels = c("NEVER","EX_SMOKER","CURRENT"), labels = c("never", "former", "current"))]
  scapis[,education := factor(cqed001, levels = c(0,1,2,3), labels = c("Compulsory",  "Compulsory", "Upper secondary", "University"))]
  scapis[,leisurePA := factor(cqpa012, levels = c(0,1,2,3),  labels = c("mainly sedentary",   "light intensity",   "regular moderate intensity ",   "regular high intensity"))]
  
  
  # Follow-up /Coverage ####
  scapis[, followup :=  as.Date(Visit1) - as.Date("2005-07-01")]
  
 
  # Restrict to participants with gut microbiota data 
  scapis <- scapis[Subject %in% microb$Subject,] # 9816
  
  
  ### Exclusion of participants based on the phenotype data ###
  
  # Ensure that all participants have a minimum follow up of 8 years. 
  scapis <- scapis[as.Date(scapis$Visit1) > as.Date("2013-07-01"), ]
  
  # Remove participants with interval between visit1 and 2 greater than 60 days. 
  scapis <- scapis[-which(Day1_Day2>60), ]
  
  # Exclude participants with chronic pulmonary disease or IBD
  scapis <-  scapis[-which(scapis$Prev_COPD=="YES"), ]
  scapis <-  scapis[-which(scapis$cqhe066=="YES"), ]
  
  # Participants whose samples were received more than 7 days after Study visit 2
  scapis <- scapis[-which(Day2_DayFecal>7), ] #9250
  
  
  # Drug Register Data ------------------------------------------------------------------------------------------------------------------
  
  # This dataset is in long format, meaning that every row is a prescription. SCAPIS participants are represented in multiple lines
  atb <-  drugs[grep("^J01", ATC),] 
  
    # Change the variable name from ID to Subject, as this is the name in phenotype file
  
  # Restricting to participants included in the phenotype data. 
  atb <- atb[Subject %in% scapis$Subject,] 
  
  # Antibiotic classes ####
  atb[grep("^J0",ATC), class:="other"] 
  atb[grep("^J01D[B,C,D,E]",ATC), class:="Class_cephalosporins"]
  atb[grep("^J01FA",ATC), class:="Class_macrolides"] 
  atb[grep("^J01FF",ATC), class:="Class_lincosamides"]
  atb[grep("^J01MA",ATC), class:="Class_FQs"]
  atb[grep("^J01A",ATC), class:="Class_TCLs"]
  atb[grep("^J01E",ATC), class:="Class_SMZTMP"]
  atb[grep("^J01XE",ATC), class:="Class_NIT"]
  atb[grep("^J01CA",ATC), class:="Class_Peni_Ext"]
  atb[grep("^J01CE",ATC), class:="Class_Peni_BetaS"]
  atb[grep("^J01CF",ATC), class:="Class_Peni_BetaR"]
  atb[grep("^J01CR",ATC), class:="Class_Peni_Comb"]
  
  
  # Merge antibiotic data with dates and intervals from phenotype
  atb <- merge(atb, scapis[,c("Subject","Visit1", "Visit2","DayFecalReceived","Day1_DayFecal","Day2_DayFecal","Day1_Day2")],  by = "Subject", all=T)
  atb[ , c("Visit1", "Visit2","DayFecalReceived","EDATUM") := lapply(.SD,as.Date), .SDcols = c("Visit1", "Visit2","DayFecalReceived","EDATUM") ] 

  
  # Create object with the prescription after 7 days of visit 2 
  atb_after_microb <- atb[which( EDATUM>=(Visit2+7) ), ]
  
  # Ignore prescriptions that occurred after 7 days of visit 2 
  atb[EDATUM>=(Visit2+7), ATC := NA ]
  atb[EDATUM>=(Visit2+7), EDATUM := NA ]
  
  # Transform the phenotype data in data.table
  setDT(scapis)

  
  # Remove participants with antibiotic prescription between visits 
  ids.to.remove <- atb[EDATUM<=(Visit2+7) & EDATUM >= Visit1, unique(Subject)]
  
  atb <- atb[-which(Subject %in% ids.to.remove),]
  scapis <- scapis[-which(Subject %in% ids.to.remove)]
  
  # Remove long-term antibiotic users during fecal sampling ####
  atb[,AtbDay_Day1 := Visit1 - as.Date(EDATUM)] # Calculate intervals between prescription and visit 1
  atb[AtbDay_Day1<0, AtbDay_Day1 := NA]
  
  atb[, forpdddtotal :=  ANTAL * forpddd]
  
  doxy <- atb[ ATC %in% c("J01AA02") & AtbDay_Day1<=90 , .(doxy=sum(forpdddtotal),ATC,EDATUM,forpddd,AtbDay_Day1), by=Subject]
  doxy.users <- doxy[ (AtbDay_Day1<=42 & doxy>=42) | (AtbDay_Day1<=56 & doxy>=56) | (AtbDay_Day1<=84 & doxy>=84)  , unique(Subject)]
  doxy40.users <- atb[ grep("40 mg", lnmn), ] %>% filter(ATC %in% c("J01AA02") & AtbDay_Day1<=90)  %>% pull(Subject)
  
  tetra <- atb[ ATC %in% c("J01AA07") & AtbDay_Day1<=90 , .(tetra=sum(forpdddtotal),ATC,EDATUM,forpddd,AtbDay_Day1), by=Subject]
  tetracycline.users <- tetra[ (AtbDay_Day1<=42 & tetra>=21) | (AtbDay_Day1<=56 & tetra>=28) | (AtbDay_Day1<=84 & tetra>=42), unique(Subject)]
  
  lyme <- atb[ ATC %in% c("J01AA04") & AtbDay_Day1<=90 , .(lyme=sum(forpdddtotal),ATC,EDATUM,forpddd,AtbDay_Day1), by=Subject]
  lymecycline.users <- lyme[ (AtbDay_Day1<=42 & lyme>=21) | (AtbDay_Day1<=56 & lyme>=28) | (AtbDay_Day1<=84 & lyme>=42), unique(Subject)]
  
  methenamine.users <- atb[AtbDay_Day1<=90 & ATC == "J01XX05", unique(Subject)] 
  
  nitro_trime.users <- atb[AtbDay_Day1<=90 & ATC %in% c("J01XE01","J01EA01") & forpdddtotal>=22.5, unique(Subject)]
  
  longterm.users <- c(doxy.users, tetracycline.users, lymecycline.users, methenamine.users, nitro_trime.users) # 28
 
  atb <- atb[-which(Subject %in% longterm.users),] 
  scapis <- scapis[-which(Subject %in% longterm.users),]

  # Participants with prescriptions < 14 days before visit 1 
  ids.to.remove <- atb[AtbDay_Day1<=14 & AtbDay_Day1>=0 , unique(Subject)]
  
  atb <- atb[-which(Subject %in% ids.to.remove),] 
  scapis <- scapis[-which(Subject %in% ids.to.remove),]
  

  # Restrict the atb data to max 8 years before visit
  atb[AtbDay_Day1>8*365.2,ATC:=NA] 
  atb[AtbDay_Day1>8*365.2,EDATUM:=NA]
  atb[is.na(ATC), AtbDay_Day1 := NA]
  
  # Final sample size calculations
  n_final <- length(unique(atb$Subject))
  atb[,atb:=1]
  atb[is.na(ATC),atb:=0]
  temp_atb <- atb[,.(N_atb=sum(atb)),by=Subject]
  n_noatb <- nrow(temp_atb[N_atb==0,])
  n_singledose <- nrow(temp_atb[N_atb==1,])
  
    
  #### Number of participants left ####
  message(paste0("Number of participants left = ",length(unique(atb$Subject))))
  message(paste0("Number of participants that did not have any antibiotic = ",n_noatb))
  
  
   # Categorize the prescriptions according to the period when it was dispensed ####
  atb[  ,t1yr   :=  0 ]
  atb[  ,t1_4yr :=  0 ]
  atb[  ,t4_8yr :=  0 ]
  
  atb[AtbDay_Day1<365, t1yr := 1 ]
  atb[AtbDay_Day1>=365 & AtbDay_Day1 < (4*365.2), t1_4yr := 1 ]
  atb[AtbDay_Day1>=(4*365.2), t4_8yr := 1 ]

  
  # Add the number for atb courses per participant per period
  N_atb <- atb[,.(N1yr = sum(t1yr) , N1_4yr=sum(t1_4yr) , N4_8yr = sum(t4_8yr), N_all = sum(t1yr,t1_4yr,t4_8yr)) , by = Subject]
  
  N_atb[,N_all := sum(N1yr,N1_4yr,N4_8yr), by = Subject] 
  
  
  # Create loop to  count the number of antibiotics for every class ####
  
  classes <- unique(atb$class[!is.na(atb$class)])
  
  for(c in classes){
    
    atb[ !is.na(ATC) ,t1yr   :=  0 ]
    atb[ !is.na(ATC) ,t1_4yr :=  0 ]
    atb[ !is.na(ATC) ,t4_8yr :=  0 ]
    
    atb[class==c & AtbDay_Day1<365, t1yr := 1 ]
    atb[class==c & AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2), t1_4yr := 1 ]
    atb[class==c & AtbDay_Day1>=(4*365.2), t4_8yr := 1 ]
    
    
    # Add the number for atb courses per participant per period
    N_atb_class <- atb[,.(N1yr = sum(t1yr) , N1_4yr=sum(t1_4yr), N4_8yr = sum(t4_8yr), N_all = sum(t1yr,t1_4yr,t4_8yr)) , by = Subject]
    
    names(N_atb_class) <- c("Subject", paste0(c, c("_1yr","_1_4yr","_4_8yr", "_all"))) 
    
    if(all(N_atb$Subject == N_atb_class$Subject)){N_atb <- cbind(N_atb,N_atb_class[,-"Subject"])} else {
      N_atb <- merge(N_atb,N_atb_class, by = "Subject")
    }
    
  }
  
  # More recent antibiotic per participant
  atb_lastdose <- atb[,.SD[which.min(AtbDay_Day1)], by=.(Subject)]
  atb_lastdose <- atb_lastdose[,.(Subject,EDATUM,ATC,class,AtbDay_Day1)]
  
  # Categorize period of the most recent antibiotic 
  atb_lastdose[AtbDay_Day1<365 , last_atb_period:="1yr"]
  atb_lastdose[AtbDay_Day1>=365 & AtbDay_Day1 <(4*365.2) ,last_atb_period:="1_4yr"]
  atb_lastdose[AtbDay_Day1>=(4*365.2) ,last_atb_period:="4_8yr"]
  
  # Rename variables for most recent antibiotic
  
  a <- c("EDATUM","ATC","class","AtbDay_Day1")
  setnames(atb_lastdose, a, paste0("last_",a))
  
  # Merge all antibiotic data (now in wide format)
  N_atb <- merge(N_atb,atb_lastdose, by="Subject", all=T )
  N_atb$last_atb_period[is.na(N_atb$last_ATC)] <- "0"
  
  
  # Single dose ----------------------------------------------------------------
  
  singledoseSubject <-  N_atb[N_all==1,"Subject"] %>% pull(Subject)
  singledose <- atb[Subject %in% singledoseSubject & !is.na(ATC),.(Subject,AtbDay_Day1)]
  setnames(singledose, "AtbDay_Day1","singledoseAtbDay_Day1")
  
  scapis[,singledose := 0]
  scapis[Subject %in% singledoseSubject ,singledose := 1]
  scapis <- merge(scapis, singledose, by = "Subject", all=T)
  
  # Merge all antibiotic data with the phenotype and gut microbiota data (all in wide format)
  scapis <- merge(scapis, N_atb, by="Subject")
  
  
  # Patient register ####
  # Re-name the ID variable 
  setnames(patreg, "ID", "Subject")
  patreg <- patreg[UTDATUM > INDATUM,]
  
  # Limit the patient register data to the subjects in the main dataset
  patreg <- patreg[Subject %in% scapis$Subject,]  # 33010 lines 
  
  # Limit the patient register data to the discharge dates after 01 July 2005
  
  patreg <- patreg[UTDATUM > as.Date("2005-07-01"), ] # 8854
  patreg <- merge(patreg, scapis[,.(Subject, Visit1)], by ="Subject", all.x=T)
  patreg <- patreg[-which(UTDATUM >= as.Date(Visit1)), ]
  
  # Participants who have been hospitalized at least once  
  
  pathospitalized_general <- unique(patreg$Subject)
  
  # Participants who have been hospitalized with a condition that usually needs antibiotics
  
  diagnosis_atb <- c("A39","A40","A41",paste0(c("A"),390:419), "R572", "A327","B377",
                     "A021","A483",
                     "J20","J21",paste0(c("J20","J21"),0:9),
                     "A15","A19",paste0(c("A1"),50:99), "K230","K930",
                     "J440","K112","K113","K122","K140",
                     "J04","J05",paste0(c("J04","J05"),0:9),
                     "H700","H702","H709","H750","J340",
                     "H60", paste0("H60",0:3),
                     "H62", paste0("H62",0:3),"H66","H670","H671",
                     "J36","J369","J390","J391",
                     "A37",paste0("A37",1:9),
                     "B085",
                     "J01","J03","J02",paste0("J0",10:39),
                     "J86",paste0("J86",1:9),
                     paste0("J01",3:8), paste0("J1",30:89), 
                     "A065" , 
                     "J85", "J86", paste0("J8",50:69),
                     "K630","A064","K750","K770","K61", paste0("K61",0:4),
                     "K650","K659",
                     "K67",paste0("K67",0:9), "N300","N308","N309",
                     "N330","N340","N341", "N080","N10","N12","N136","N160","A985",
                     "N431","N74", paste0("N74",c(0:4,8)),
                     "N151","N159","N390","N291","M726","M600","M860","M861","M862","M869",
                     "M900","M901","M902","M00","M000","M001","M002","M008","M009","M01","M010",
                     "M011","M012","M013","M014","M015","M016","M018","M462","M463","M464","M465",
                     "M490","M491","M492","M493","L02","L03","A46","A469","L04","A067","L08",
                     "L050","L00","A481","A32","A31","B59","A020","A021","A066","G060","G07",
                     "G061","G062","A321","A390","G00","G01","A48",paste0("A48",0:9),
                     "A49", paste0("A49",0:9),"B95","B96",paste0("B9",50:69),
                     "A36", paste0("A36",0:9), "A71","A710","A711","A719","B005","H00","H03",
                     "H050", "H061", "H100", "H130","H131","H190","H191","H192","H220","H320","H440",
                     "B999","A32", "A31", "I301","I320","I321","I330","I400","I410","I411", "I412",
                     "I430","I520", "I521","B376","A395","A188","A010","A38","A65",paste0("A6",60:99),
                     paste0("A6",6:9),"D733","E060")
  
  
  pathospitalized_infect <- patreg[,c("Subject",paste0("DIA",1:30))] %>%
    pivot_longer(cols=paste0("DIA",1:30), names_to = "DIA", values_to = "ICD") %>% 
    filter(ICD %in% diagnosis_atb) %>% pull(Subject) %>% unique(.)
  
  
  # Create new variable in the final dataset 
  
  scapis[, hospgeneral:= "no"]
  scapis[Subject %in% pathospitalized_general, hospgeneral:= "yes"]
  
  scapis[, hospinfect:= "no"]
  scapis[Subject %in% pathospitalized_infect, hospinfect:= "yes"]
  
  
  # Antibiotic prescriptions AFTER visit 2 ------------------------------------------------------------------------
  # Categorize the prescriptions according to the period when it was dispensed #
  atb_after_microb[ ,after_1yr :=  0 ]
  atb_after_microb[ ,after_1_4yr :=  0 ]
  
  atb_after_microb[,Day2_Atb := as.Date(EDATUM) - Visit2 ]
  
  atb_after_microb[Day2_Atb<365, after_1yr := 1 ]
  atb_after_microb[Day2_Atb>=365 & Day2_Atb <(4*365.2), after_1_4yr := 1 ]
  
  
  N_atb_after <- atb_after_microb[,.(N_after_1yr = sum(after_1yr) , N_after_1_4yr=sum(after_1_4yr)) , by = Subject]
  
  N_atb_after[,N_after_all := sum(N_after_1yr,N_after_1_4yr), by = Subject]
  
  
  classes <- unique(atb_after_microb$class[!is.na(atb_after_microb$class)])
  
  for(c in classes){
    
    atb_after_microb[, after_1yr := 0 ]
    atb_after_microb[class==c  & Day2_Atb<365, after_1yr := 1 ]
    # atb_after_microb[class==c  & Day2_Atb>=365 & Day2_Atb <(4*365.2), after_1_4yr := 1 ]
    
    
    # Add the number for atb courses per participant per period
    N_atb_class <- atb_after_microb[,.(N1yr=sum(after_1yr)) , by = Subject]
    
    names(N_atb_class) <- c("Subject", paste0(gsub("Class_","",c), c("_after1yr"))) 
    
    if(all(N_atb_after$Subject == N_atb_class$Subject)){N_atb_after <- cbind(N_atb_after,N_atb_class[,-"Subject"])} else {
      N_atb_after <- merge(N_atb_after,N_atb_class, by = "Subject")
    }
    
  }
  
  scapis <- merge(scapis , N_atb_after, by="Subject", all.x=T, all.y=F)
  
  
  # Alpha diversity metrics -------------------------------------------------------------------------------------------------
  
  downmgs <- fread('/proj/sens2019512/SCAPIS/Gutsy/Metagenomics/Processed/scapis_metagenomics_mgs_ds_relative_abundances_v1.0.tsv', data.table=F, na.strings = c("NA",NA,""))
  rownames(downmgs) <- downmgs$scapis_id
  downmgs$scapis_id <- NULL
  
  richness <- data.frame(Subject = rownames(downmgs), richness = rowSums(downmgs > 0))
  
  scapis <- merge(scapis, richness, by="Subject")
  scapis <- merge(scapis, alpha, by.x="Subject", by.y = "scapis_id")
  
  # Save final phenotype data ####
  fwrite(scapis, file = "/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/scapis_working_dataset.tsv")
  
  message("working dataset saved")
  
  # Calculate species prevalence -------------------------------------------------------------------------------------------------
  # Main model covariates 
  
  basic.model <- c("age","Sex","placebirth","smokestatus","education","site_plate")
  full.model <- c(basic.model, "BMI","diabd","rheumatic","cancer", "ppi")
  
  cc <- complete.cases(scapis[,basic.model, with=F])
  
  species <- ifelse(as.matrix(downmgs[scapis[cc,][["Subject"]],]) < 0.01, 0, 1)
  prevalence <- colSums(species)/nrow(species)
  
  prevalent.species <- names(prevalence[prevalence >.01]) 
  non.prevalent.species <- names(prevalence[prevalence <=.01])

  
  save(list=c("basic.model","prevalent.species","non.prevalent.species","full.model"),
       file='/proj/sens2019512/nobackup/users/baldanzi/atb_gut/work/scapis_model.Rdata')
  
  message(paste0("Main model saved = ",paste0(basic.model,collapse = ", ")))
  message(paste0("Length of prevalent species = ",length(prevalent.species)))
  
  
  message("End")
  

  
  