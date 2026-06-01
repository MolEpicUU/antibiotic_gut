# Project Antibiotic Use and the Gut Microbiota

# Table 2

  rm(list=ls())
  
  setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed')
  
  library(data.table)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(patchwork)

  scapis <- fread('revision__work/scapis_atb_deidentif_data.tsv')
  setnames(scapis, "Sex", "sex")
  simpler <- fread('revision__work/simpler_atb_deidentif_data.tsv')
  mos <- fread('revision__work/mos_atb_deidentif_data.tsv')


# Table 2 ---------------------------------------------------------------------

  scapis_all  <- scapis[, grep("^Class.*yr", names(scapis), v=T), with=F]
  simpler_all <- simpler[, grep("Class.*yr", names(simpler), v=T), with=F]
  mos_all     <- mos[, grep("Class.*yr", names(mos)), v=T, with=F]
  
  unique_atb <- unique(str_extract(grep("^Class.*yr", names(scapis), v=T), "Class_.*_"))
  
  
    list_cohorts <- lapply(list(scapis_all, simpler_all, mos_all), function(cohort){
    
    for(atb in unique_atb){
      setnames(cohort, paste0(atb, "1yr"), paste0(atb, "period1"))
      
      var <- paste0(atb, c(2:4), "yr")
      cohort$temp <- rowSums(cohort[, var , with=F])
      setnames(cohort, "temp", paste0(atb, "period2"))
      cohort <- cohort[, !c(var), with=F]
      
      var <- paste0(atb, c(5:8), "yr")
      cohort$temp <- rowSums(cohort[, var , with=F])
      setnames(cohort, "temp", paste0(atb, "period3"))
      cohort <- cohort[, !c(var), with=F]
    }
      return(cohort)
    })
  
    
    scapis_all <- list_cohorts[[1]]
    simpler_all <- list_cohorts[[2]]
    mos_all <- list_cohorts[[3]]
    
  
    
  scapis_all <- scapis_all[, .(Antibiotic = names(.SD), 
                               SCAPIS = paste0(lapply(.SD, function(x) sum(x>0)), " (",
                                          lapply(.SD, function(x) round(sum(x>0)*100/.N,1)), "%)"))]
  
  simpler_all <- simpler_all[, .(Antibiotic = names(.SD), 
                               SIMPLER = paste0(lapply(.SD, function(x) sum(x>0)), " (",
                                               lapply(.SD, function(x) round(sum(x>0)*100/.N,1)), "%)"))]
  
  mos_all <- mos_all[, .(Antibiotic = names(.SD), 
                                 MOS = paste0(lapply(.SD, function(x) sum(x>0)), " (",
                                                  lapply(.SD, function(x) round(sum(x>0)*100/.N,1)), "%)"))]
  
  
  
  tab2 <- merge(scapis_all, simpler_all, by = 'Antibiotic')
  tab2 <- merge(tab2, mos_all, by = 'Antibiotic')
    
  tab2$Antibiotic  
  
  tab2[, period := str_extract(Antibiotic, "period.")]   
  tab2[, period := factor(period, c("period1", "period2", "period3"), c("<1yr", "1-4yr", "4-8yr"))]
  tab2[, Antibiotic := gsub("Class_", "", gsub("_period.", "", Antibiotic))]
  tab2 <- tab2[Antibiotic != "other", ]
  tab2[, Antibiotic := factor(Antibiotic, c("Nall","Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb","cephalosporins","macrolides",
                               "lincosamides","TCLs","FQs","SMZTMP","NIT"),
                             c("Any","Penicillin V","Flucloxacillin","Penicillins with extended spectrum","Amoxicillin-clavulanic acid",
                               "Cephalosporins","Macrolides","Clindamycin","Tetracyclines","Fluoroquinolones","Sulfamethoxazole-trimethoprim",
                               "Nitrofurantoin"))]

 tab2 <- tab2[order(Antibiotic), ]
 tab2 <- pivot_wider(tab2, id_cols = Antibiotic, names_from = period, values_from = c("SCAPIS", "SIMPLER", "MOS"), names_vary = "fastest")
 tab2
 
 fwrite(tab2, 'revision__tables/revision_table_2.tsv', sep = '\t')

