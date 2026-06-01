# Project Antibiotic Use and the Gut Microbiota

# This script will classify the antibiotics most used overall and by period 

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


# Figure 2 - Number of antibiotic users ---------------------------------------------------------------------

  scapis_all  <- scapis[, grep("^Class.*yr", names(scapis), v=T), with=F]
  simpler_all <- simpler[, grep("Class.*yr", names(simpler), v=T), with=F]
  mos_all     <- mos[, grep("Class.*yr", names(mos)), v=T, with=F]
  cohorts_all <- rbind(scapis_all, simpler_all, mos_all)
  
  Ntotal <- nrow(cohorts_all)
  
  unique_atb <- unique(str_extract(grep("^Class.*yr", names(cohorts_all), v=T), "Class_.*_"))


  for(atb in unique_atb){
    
    var <- paste0(atb, c(1:8), "yr")
    cohorts_all$temp <- rowSums(cohorts_all[, var , with=F])
    setnames(cohorts_all, "temp", paste0(atb, "period4"))
    
    setnames(cohorts_all, paste0(atb, "1yr"), paste0(atb, "period1"))
    
    var <- paste0(atb, c(2:4), "yr")
    cohorts_all$temp <- rowSums(cohorts_all[, var , with=F])
    setnames(cohorts_all, "temp", paste0(atb, "period2"))
    
    var <- paste0(atb, c(5:8), "yr")
    cohorts_all$temp <- rowSums(cohorts_all[, var , with=F])
    setnames(cohorts_all, "temp", paste0(atb, "period3"))
    
    var <- paste0(atb, c(2:8), "yr")
    cohorts_all <- cohorts_all[, !c(var), with=F]
  }
  
  cohorts_all <- cohorts_all[, .(antibiotic = names(.SD), number_users = as.numeric(lapply(.SD, function(x) sum(x>0))))]
  
  cohorts_all[, c("antibiotic", "period") := tstrsplit(antibiotic, "_period")]
  cohorts_all[, antibiotic := gsub("Class_", "", antibiotic)]
  cohorts_all[, period := factor(period, c("1","2","3", "4"), c("<1y", "1-4y", "4-8y", "total_8years"))]
  
  cohorts_all[, proportion_users := number_users/Ntotal]

  antibiotics_order <- cohorts_all[period == "total_8years", .(antibiotic, antibiotic_order = rank(-number_users))]
  
  cohorts_all <- merge(cohorts_all, antibiotics_order, by = "antibiotic")
  
  
  abx_names <- c("Penicillin V", "Tetracyclines", "Penicillin extended spectrum", "Flucloxacillin", 
                 "Fluoroquinolones", "Nitrofurantoin", "Clindamycin", "SMZ-TMP", "Cephalosporins",
                 "Macrolides", "Amox-clav")
  
  abx_abbrev <- antibiotics_order[order(antibiotic_order), antibiotic]
  
  cohorts_all[, antibiotic := factor(antibiotic, abx_abbrev, abx_names)]
  cohorts_all <- cohorts_all[order(antibiotic_order, period)]
  setcolorder(cohorts_all, c("antibiotic_order", "antibiotic", "period", "number_users", "proportion_users"))

  cohorts_all
    
  fwrite(cohorts_all, file = "revision__work/antibiotics_number_of_users.tsv", sep = "\t")
    
    
    
  
  