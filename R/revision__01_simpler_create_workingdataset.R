# Project Antibiotic Use and the Gut Microbiota

rm(list=ls())

  # Version date: 2023-11-13
  
  # Load packages
  library(data.table)
  library(dplyr)
  library(tidyr)

# This is the script will:
# 1. Import data on SIMPLER simplers and prepare the data
# 2. Merge phenotype data, gut microbiota data, clinical microbiomics data, place of birth
# 3. Import the prescribed drug register data and prepare the data
# 4. Merge phenotype, gut microbiota, and drug register data 

# Import phenotype datasets 
  simpler   <- fread('/proj/simp2023007/data/phenotypes/FAECES7311infoUPDATED.txt', na.strings = c("NA",NA,""))
  # age, sex, country of birth, site, date fecal sample collection, year fecal sample collection 
  simpler_edu_smok <- fread('/proj/simp2023007/Dataleverans/main_quest_data.csv', na.strings = c("NA",NA,""))
  diabdrug  <- fread('/proj/simp2023007/Dataleverans/diabetes_drugs.csv', na.strings = c("NA",NA,""))
  cancer    <- fread('/proj/simp2023007/Dataleverans/cancer.csv', na.strings=c("NA","",NA))
  patreg    <- fread('/proj/simp2023007/Dataleverans/diseases.csv', na.strings=c("NA","",NA)) # 2567
  hospgeneral <- fread('/proj/simp2023007/Dataleverans/inpat2005.csv', na.strings=c("NA","",NA))
  rheumat_diag <- fread('/proj/simp2023007/Dataleverans/extra_diagnosis_dates.csv')
  polypharmacy_med <- fread('/proj/simp2023007/Dataleverans/medications_for_medscore.csv', na.strings = c("NA", "", NA))
  medication <- fread('/proj/simp2023007/Dataleverans/medscore.csv', na.strings = c("NA", "", NA))
  cci <- fread('/proj/simp2023007/Dataleverans/CCI.csv', na.strings = c("NA", "", NA))
  
# Import other exposure and outcome datasets 
  drugs   <- fread('/proj/simp2023007/data/phenotypes/prescibed_drugs.csv', na.strings = c("NA",NA,"")) # 110331
  
  microb  <- fread('/proj/simp2023007/data/microbiome/processed/simpler_metagenomics_mgs_relative_abundances_v2.0.tsv') # 6150
  downmgs <- fread('/proj/simp2023007/data/microbiome/processed/simpler_metagenomics_mgs_ds_relative_abundances_v2.0.tsv', data.table=F, na.strings = c("NA",NA,""))
  alpha   <- fread('/proj/simp2023007/data/microbiome/processed/simpler_metagenomics_alpha_diversity_v2.0.tsv', na.strings=c("NA","",NA)) 
  techvar <- fread('/proj/simp2023007/data/microbiome/processed/simpler_metagenomics_technical_variables_v2.0.tsv', na.strings=c("NA","",NA)) 
  
  # Data management -------------------------------------------------
  setnames(simpler, "dateofsampling", "Visit1")
  simpler[, Visit1 := as.Date(tolower(Visit1), format='%d%b%Y') ]
  simpler[, placebirth :=  factor(placeofbirth, c("Sverige", "Norden utom Sverige", "EU28 utom Norden", "Europa utom EU28 och Norden", 
                                                    "Sovjetunionen", "Asien","Sydamerika","Afrika"), 
                                    c("Scandinavia","Scandinavia","Europe", "Europe","Europe","Asia","Others","Others"))]
  simpler <- simpler[, center := factor(center, c("Uppsala","Sthlm","V\344ster\345"), c("Uppsala","Stockholm","Vasteras"))]
  simpler <- simpler[SIMPKEY %in% simpler_edu_smok$SIMPKEY, ]
  simpler <- simpler[location %in% "Klinik",]
  
  # Smoking ####
  simpler_edu_smok <- simpler_edu_smok[SIMPKEY %in% simpler$SIMPKEY, ]
  simpler_edu_smok <- simpler_edu_smok[match(simpler$SIMPKEY, SIMPKEY), ]
  simpler_edu_smok$yrvisit <- format(simpler$Visit1, format="%Y")
  
  simpler_edu_smok[D_smok01_u==1, D_smok13:=0]
  simpler_edu_smok[F_smok01_u==1, F_smok13:=0]
  simpler_edu_smok[D_smok13>0 & (D_smok01_u!=3 | is.na(D_smok01_u)), D_smok01_u:=2,]
  simpler_edu_smok[F_smok13>0 & (F_smok01_u!=3 | is.na(F_smok01_u)), F_smok01_u:=2,]
  simpler_edu_smok[, B_der_smok_u := factor(B_der_smok_u, 1:3, c("current","ex","never"))]
  simpler_edu_smok[B_der_smok_u == "never" & is.na(D_smok01_u), D_smok01_u := 1 ]
  simpler_edu_smok[B_der_smok_u == "never" & is.na(F_smok01_u), F_smok01_u := 1 ]
  simpler_edu_smok[is.na(F_smok01_u), F_smok01_u := D_smok01_u]
  
  simpler_edu_smok[, D_smok01_u := factor(D_smok01_u, 1:3, c("never","current","ex"))]
  simpler_edu_smok[, F_smok01_u := factor(F_smok01_u, 1:3, c("never","current","ex"))]
  
  
  simpler_edu_smok[yrvisit<2019, smokestatus:= D_smok01_u]
  simpler_edu_smok[yrvisit<2019 & is.na(smokestatus), smokestatus:= B_der_smok_u]
  simpler_edu_smok[yrvisit<2019 & is.na(smokestatus), smokestatus:= F_smok01_u]
  simpler_edu_smok[yrvisit>=2019, smokestatus:= F_smok01_u]
  simpler_edu_smok[yrvisit>=2019 & is.na(smokestatus), smokestatus:= D_smok01_u]
  simpler_edu_smok[yrvisit>=2019 & is.na(smokestatus), smokestatus:= B_der_smok_u]
  simpler_edu_smok[is.na(smokestatus) & (D_smok13==0 | F_smok13 ==0), smokestatus:= "never"]
  
  simpler_edu_smok[, education := factor(B_der_educ1, 1:3, c("Compulsory", "Upper secondary","University"))]
  
  simpler <- merge(simpler, simpler_edu_smok[,.(SIMPKEY, smokestatus, education) ])
  
  # DNA extraction plate ####
  # techvar <-   techvar[SIMPKEY %in% simpler$SIMPKEY, ]
  # techvar <- techvar[match(simpler$SIMPKEY, SIMPKEY), ]
  simpler <- merge(simpler, techvar[, .(SIMPKEY, aliquoting_plate)], by="SIMPKEY")
  
  
  # Diabetes diagnosis ####
  # Based on the prescribed drug register 
  # At least one prescription for diabetes medication in the 18m before sampling 
  diabdrug <- merge(diabdrug, simpler[, .(Visit1, SIMPKEY)], by="SIMPKEY")
  diabdrug <- diabdrug[as.Date(EDATUM) <= as.Date(Visit1), ]
  diabdrug[, diff := as.Date(Visit1) - as.Date(EDATUM)]
  
  diabdrug_18 <- diabdrug[diff < 1.5*365.2, unique(SIMPKEY), drop=T]
  simpler[, diabd := ifelse(SIMPKEY %in% diabdrug_18, "yes","no")]
  
  
  # Rheumatic, IBD, COPD #### 
  simpler[ , rheumatic  := ifelse(SIMPKEY %in% rheumat_diag[, SIMPKEY], "yes","no")]
  simpler[ , IBD   := ifelse(SIMPKEY %in% patreg[IBD == 1, SIMPKEY]    , "yes","no")] 
  simpler[ , COPD  := ifelse(SIMPKEY %in% patreg[J424344 == 1, SIMPKEY], "yes","no")]
  
  # Cancer
  # Cancer diagnosis 
  cancer <- cancer[-grep("^C44", cancer$ICDO10),]
  cancer <- merge(cancer, simpler[, .(SIMPKEY, Visit1)], by="SIMPKEY")
  cancer[ , diadate:= as.Date(tolower(diadate), format="%d%b%Y")]
  cancer <- cancer[ diadate<Visit1, ]
  cancer[ , diff:= as.numeric((Visit1 - diadate)/365.25)]
  cancer <- cancer[ diff <8 , ]
  cancer <- cancer[, .SD[which.min(diff)], by=.(SIMPKEY)]
  cancer[ diff<1           , cancer := "1 year before V1"]
  cancer[ diff>=1 & diff<4 , cancer := "2-5 year before V1"]
  cancer[ diff>=4 & diff<8 , cancer := "6-8 year before V1"]
  
  simpler <- merge(simpler, cancer[,.(SIMPKEY, cancer)], by = "SIMPKEY", all.x=T)
  simpler[is.na(cancer) , cancer := "no"]
  
  # PPI ####
  ppi_reg <-  copy(drugs[grep("^A02BC", ATC), ])
  ppi_reg  <- merge(ppi_reg, simpler[, .(SIMPKEY, Visit1)], by="SIMPKEY")
  ppi_reg <- ppi_reg[as.Date(Visit1) - as.Date(EDATUM) <= 365.2 & as.Date(Visit1) - as.Date(EDATUM) > 0, ]
  simpler[ , ppi := ifelse(SIMPKEY %in% ppi_reg$SIMPKEY, "yes","no")]
  
  # Statins
  statins_reg <-  medication[C10AA == 1, SIMPKEY]
  simpler[ , statins := ifelse(SIMPKEY %in% statins_reg, "yes","no")]
  
  # Metformin ####
  metformin_reg <-  copy(drugs[grep("^A10BA", ATC), ])
  metformin_reg  <- merge(metformin_reg, simpler[, .(SIMPKEY, Visit1)], by="SIMPKEY")
  metformin_reg <- metformin_reg[as.Date(Visit1) - as.Date(EDATUM) <= 365.2 & as.Date(Visit1) - as.Date(EDATUM) > 0, ]
  simpler[ , metformin := ifelse(SIMPKEY %in% metformin_reg$SIMPKEY, "yes","no")]
  
  # betablock ####
  betablock_reg <-  copy(drugs[grep("^C07AB", ATC), ])
  betablock_reg  <- merge(betablock_reg, simpler[, .(SIMPKEY, Visit1)], by="SIMPKEY")
  betablock_reg <- betablock_reg[as.Date(Visit1) - as.Date(EDATUM) <= 365.2 & as.Date(Visit1) - as.Date(EDATUM) > 0, ]
  simpler[ , betablock := ifelse(SIMPKEY %in% betablock_reg$SIMPKEY, "yes","no")]
  
  # SSRI 
  ssri_reg <-  medication[N06AB == 1, SIMPKEY]
  simpler[ , ssri := ifelse(SIMPKEY %in% ssri_reg, "yes","no")]
  
  # betablock ####
  antipsycho_reg <-  copy(polypharmacy_med[grep("^N05A", ATC), ])
  antipsycho_reg  <- merge(antipsycho_reg, simpler[, .(SIMPKEY, Visit1)], by="SIMPKEY")
  antipsycho_reg <- antipsycho_reg[as.Date(Visit1) - as.Date(EDATUM) <= 365.2 & as.Date(Visit1) - as.Date(EDATUM) > 0, ]
  simpler[ , antipsycho := ifelse(SIMPKEY %in% antipsycho_reg$SIMPKEY, "yes","no")]
  
  # Hospitalization with infection 
  hospinfect  <- patreg[!is.na(UTDATUM), ]
  hospinfect  <- merge(hospinfect, simpler[, .(SIMPKEY, Visit1)], by="SIMPKEY")
  hospinfect  <- hospinfect[as.Date(INDATUM)<Visit1, ]
  hospinfect  <- hospinfect[INDATUM>as.IDate("2005-07-01"), ]
  hospinfect  <- hospinfect[,.(SIMPKEY, sepsis, ear, lung, softmuscle, gastro, urogen, brain, kardit, tuberkmyco, otherinfec) ] %>% 
    mutate(infection = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% filter(infection>=1)
  
  simpler[, hospinfect := ifelse(SIMPKEY %in% hospinfect$SIMPKEY, "yes","no")]
  
  
  # Hospitalization general 
  hospgeneral[, indate := as.Date(tolower(indate), format = "%d%b%Y")]
  hospgeneral   <- hospgeneral[ indate> as.Date("2005-07-01"), ]
  hospgeneral   <- merge(hospgeneral, simpler[, .(SIMPKEY, Visit1)], by="SIMPKEY")
  hospgeneral   <- hospgeneral[ indate<Visit1, ]
  simpler[ , hospgeneral := ifelse(SIMPKEY %in% hospgeneral$SIMPKEY, "yes","no")]
  
  # Follow-up/Coverage
  simpler[ , followup := Visit1 - as.Date("2005-07-01")]
  
  
  # Charlson comorbidity index weighted #### -------------------------------------------------------------
  cci = merge(cci, simpler[, .(SIMPKEY, metformin)], by = "SIMPKEY")
  cci[diabetes == 0 & metformin == "yes", CCI_weighted := CCI_weighted + 1]
  
  setnames(cci, "CCI_weighted", "CCIw")
  simpler <- merge(simpler, cci[, c("SIMPKEY", "CCIw")], by="SIMPKEY", all.x=T )
  
  
  # Polypharmacy ####
  to_keep <- unique(c( grep('^A0[2,3,4,5,6,7,8,9]', polypharmacy_med$ATC, v=T),
                       grep('^A[10,14,15,16]', polypharmacy_med$ATC, v=T),
                       grep('^[C,G,H]', polypharmacy_med$ATC, v=T),
                       grep('^J0[2,4,5,6]', polypharmacy_med$ATC, v=T),
                       grep('^M0[1,3,5,6]', polypharmacy_med$ATC, v=T),
                       grep('^[L, N, P, R03]', polypharmacy_med$ATC, v=T)))
  
  to_keep <- to_keep[!to_keep %in% c('N02BE01', 'M01AE')]
  to_keep <- to_keep[!grepl("^C05", to_keep)]
  
  polypharmacy_med <- polypharmacy_med[ATC %in% to_keep,]
  
  polypharmacy_med <- merge(polypharmacy_med, simpler[,c("SIMPKEY","Visit1")],  by = "SIMPKEY", all.y=T)
  
  polypharmacy_med <- polypharmacy_med[EDATUM < Visit1, ]
  polypharmacy_med <- polypharmacy_med[EDATUM > (Visit1 - 365), ]
  
  poly6m <- polypharmacy_med[EDATUM > (Visit1 - 180), .N , by = c("ATC","SIMPKEY")][N>1,]
  poly6m <- poly6m[, .(polypharmacy6m = .N), by = "SIMPKEY"]
  
  simpler <- merge(simpler, poly6m,  by = "SIMPKEY", all.x=T)
  simpler[is.na(polypharmacy6m),  polypharmacy6m := "0"]
  simpler[polypharmacy6m == 0, polypharmacy6m_cat := "0"]
  simpler[polypharmacy6m > 0 , polypharmacy6m_cat := "1-4"]
  simpler[polypharmacy6m > 4 , polypharmacy6m_cat := "5 or more"]
  
  poly12m <- polypharmacy_med[, .N , by = c("ATC", "SIMPKEY")][N>2,]
  poly12m <- poly12m[, .(polypharmacy12m = .N), by = "SIMPKEY"]
  
  simpler <- merge(simpler, poly12m,  by = "SIMPKEY", all.x=T)
  simpler[is.na(polypharmacy12m), polypharmacy12m := "0"]
  simpler[polypharmacy12m == 0, polypharmacy12m_cat := "0"]
  simpler[polypharmacy12m > 0 , polypharmacy12m_cat := "1-4"]
  simpler[polypharmacy12m > 4 , polypharmacy12m_cat := "5 or more"]
  
  message("Polypharmacy done")
  
  # simpler <- merge(simpler, polypharmacy[, c("SIMPKEY", "polypharmacy12m", "polypharmacy6m")], by ="SIMPKEY", all.x=T)
  # simpler[is.na(polypharmacy6m),  polypharmacy6m := "0"]
  # simpler[polypharmacy6m == 0, polypharmacy6m_cat := "0"]
  # simpler[polypharmacy6m > 0 , polypharmacy6m_cat := "1-4"]
  # simpler[polypharmacy6m > 4 , polypharmacy6m_cat := "5 or more"]
  # simpler[is.na(polypharmacy12m),  polypharmacy12m := "0"]
  # simpler[polypharmacy12m == 0, polypharmacy12m_cat := "0"]
  # simpler[polypharmacy12m > 0 , polypharmacy12m_cat := "1-4"]
  # simpler[polypharmacy12m > 4 , polypharmacy12m_cat := "5 or more"]
  
  
  fwrite(simpler, "/proj/nobackup/simp2023007/users/baldanzi/work/simpler_noexclusions.tsv") #5888
  
  # simpler <- fread("/proj/nobackup/simp2023007/users/baldanzi/work/simpler_noexclusions.tsv", , na.strings = c("NA",NA,""))

  # Restrict to participants with gut microbiota data 
  simpler <- simpler[SIMPKEY %in% microb$SIMPKEY, ] # 5888 # It will increase when using microbiota_v1.1
 
  # Ensure that all participants have a minimum follow up of 8 years. 
  simpler[ as.Date(Visit1)  <= as.Date("2013-07-01"), .N]
  simpler <- simpler[ as.Date(Visit1)  > as.Date("2013-07-01"), ] # 5442
  
  # Exclusion based on the phenotype data 
  simpler <- simpler[ COPD == "no", ]
  simpler <- simpler[ IBD  == "no", ]  
  
  # Antibiotic drug management --------------------------------------------------
  atb <-  copy(drugs[grep("^J01", ATC), ]) # 42 559
  length(unique(atb$SIMPKEY)) # 5637
  
  # Restricting to participants included in the phenotype data. 
  atb <- atb[SIMPKEY %in% simpler$SIMPKEY,] # 33278
  length(unique(atb$SIMPKEY)) # 4673
  
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
  
  # message("\nprescriptons before 14 days")
  # # Participants with prescriptions < 14 days before visit 1 
  # ids.to.remove <- atb[AtbDay_Day1<=14 , unique(SIMPKEY)]
  # 
  # atb <- atb[-which(SIMPKEY %in% ids.to.remove), ]
  # simpler <- simpler[-which(SIMPKEY %in% ids.to.remove), ]
  # 
  # message(paste(length(ids.to.remove), "prescription 14 days before fecal sampling")) # 41
  # length(unique(atb$SIMPKEY)) #  5056 
  
  # Restrict the atb data to max 8 years before visit
  atb[AtbDay_Day1>8*365.2,ATC:=NA] 
  atb[AtbDay_Day1>8*365.2,EDATUM:=NA]
  atb[is.na(ATC), AtbDay_Day1 := NA]
  
  # Final sample size calculations
  atb <- merge(simpler[,"SIMPKEY"], atb, by = "SIMPKEY", all=T)
  n_final <- length(unique(atb$SIMPKEY)) # 5056
  atb[,atb:=1]
  atb[is.na(ATC),atb:=0]
  temp_atb <- atb[,.(N_atb=sum(atb)),by=SIMPKEY]
  n_noatb <- nrow(temp_atb[N_atb==0,]) # 1315
  
  #### Number of participants left ####
  message(paste0("Number of participants left = ",length(unique(atb$SIMPKEY))))
  message(paste0("Number of participants that did not have any antibiotic = ",n_noatb))
  
  
  # Save the processed prescribed drug register data 
  fwrite(atb,"/proj/nobackup/simp2023007/users/baldanzi/work/simpler_processed_lakemedelregister_revision.tsv")
  # atb <- fread("/proj/nobackup/simp2023007/users/baldanzi/work/simpler_processed_lakemedelregister.tsv")
  
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
  simpler <- merge(simpler, N_atb, by="SIMPKEY") # 5108
  
  
  # Alpha diversity -----------------------------------------------------------
  richness <- data.frame(SIMPKEY = downmgs$SIMPKEY, richness = rowSums(downmgs[,-1] > 0))
  
  simpler <- merge(simpler, richness, by="SIMPKEY")
  simpler <- merge(simpler, alpha , by="SIMPKEY")
  
  # Save final simpler data ####
  fwrite(simpler, file = "/proj/nobackup/simp2023007/users/baldanzi/work/simpler_working_dataset_revision_noatbexclusion.tsv", sep="\t")
  
  
  simpler[ last_EDATUM >= (Visit1-30), .N]
  simpler[is.na(Visit1), .N]
  simpler <- simpler[is.na(last_EDATUM) | last_EDATUM < (Visit1-30), ]
  fwrite(simpler, file = "/proj/nobackup/simp2023007/users/baldanzi/work/simpler_working_dataset_revision.tsv", sep="\t")
  
  # Calculate species prevalence -------------------------------------------------------------------------------------------------
  # Main model covariates 
  basic.model <- c("age", "sex", "placebirth", "smokestatus", "education", "center","aliquoting_plate")
  #full.model <- c(main.model, "BMI","CCIw", "diabd","rheumatic","cancer", "ppi")
  
  full.model <- c(basic.model, "BMI", "CCIw", "ppi", "statins", "metformin", "betablock", "ssri", "polypharmacy12m_cat", "antipsycho")
  
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
  
  # Prevalence table 
  prev_table <- data.table(species = names(prevalence), prevalence_simpler = prevalence)
  
  fwrite(prev_table, '/proj/simp2023007/nobackup/users/baldanzi/results/prevalence_species_simpler.tsv', sep = '\t')
  
  
  
  
  message("End")
  
  
  