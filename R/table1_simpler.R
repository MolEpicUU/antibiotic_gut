# Project Antibiotic and the gut microbiota 

# Gabriel Baldanzi 

# 2025-05-31

# Table 1. Population descriptive table - SIMPLER

rm(list=ls())

suppressMessages({library(compareGroups)
library(data.table)})
library(Hmisc)
library(sjmisc)

  

  # Import data set
  setwd('/proj/nobackup/simp2023007/users/baldanzi/')
  simpler <- fread("work/simpler_working_dataset_revision.tsv", na.strings = c("NA", NA, ""))


  # Import models 
  load('work/simpler_model_revision.Rdata')
  
  cc <- complete.cases(simpler[,basic.model,with=F])
  simpler <- simpler[cc, ]
  
  simpler[, atb_cat:=ifelse(N_all >=2 , 2, N_all )]
  simpler[, atb_cat:=factor(atb_cat, c(0,1,2), c("none", "1 course",">=2 courses"))]

  simpler[, sex := factor(sex, c("M","F"), c("male","female") )]
  simpler[, placebirth := factor(placebirth, c("Scandinavia", "Europe","Asia","Others"), c("Scandinavia","no","no","no") )]
  simpler[, smokestatus := factor(smokestatus, c("never","ex","current"), c("Never","Former","Current"))]
  simpler[, diabd := factor(diabd,  c("no","yes") )]
  simpler[, education:=factor(education, c("Compulsory", "Upper secondary", "University"))]
  simpler[, polypharmacy_cat := fifelse(polypharmacy12m_cat == "5 or more", "yes", "no")]
  simpler[, CCI2 := fifelse(CCIw >= 2, "yes", "no")]

  myvars <- c("age" , "sex", "center",
            "smokestatus", "BMI", "education",
            "placebirth", "diabd", 
            "atb_cat", "N_all","rheumatic","cancer", "ppi",
            "CCI2", "metformin", "betablock", "ssri", "statins", "antipsycho", 
            "polypharmacy_cat",
            "hospinfect","hospgeneral", "aliquoting_plate")

  mylabels <- c("Age (years)", "Sex", "Study site", "Smoking",  
              "BMI (kg/m2)",  "Education", 
              "Born in Scandinavia", "Diabetes",
              "Number of ATB courses", "Number of antibiotic courses",
              "Rheumatic disease","Previous cancer", "PPI",
              "Charlson index >=2", "Metformin", "Beta-blocker", "SSRI", "Statins", "Antipsychotics", 
              "Polypharmacy",
              "Hospitalization for infection", "Any hospitalization", "Aliquoting plate")

  temp.data <- simpler[,myvars,with=F]

  temp.data[,Female:=ifelse(sex=="female","yes","no")]
  temp.data[!is.na(cancer),cancer:=ifelse(cancer=="no","no","yes")]


  message("Assign labels")

# Variables for table 1 


  
  j=1
  for(i in myvars){
    label(temp.data[[i]]) <- mylabels[j]
    j=j+1
  }

# Table 1 part 1 ####

  t1 <- compareGroups(atb_cat ~  age + Female + smokestatus + education + placebirth,
                      include.miss = FALSE, chisq.test.perm = TRUE,
                      data = temp.data, method=2)
  
  
  t = createTable(t1, hide.no = "no",
                  digits = 1, show.n = T)
  
  
  export2csv(t, file = '/proj/nobackup/simp2023007/wharf/baldanzi/baldanzi-simp2023007/atbgut/results/Table1_part1_simpler.csv', 
             which.table="descr", sep="\t", nmax = TRUE)
  

# Table 1 part 2 ####
  full.model <- full.model[-which(full.model %in% c("CCIw", "polypharmacy12m_cat"))]
  cc <- complete.cases(temp.data[,c(full.model, "CCI2", "polypharmacy_cat"),with=F])
  temp.data <- temp.data[cc, ]

t1 <- compareGroups(atb_cat ~ 
                      BMI + CCI2 + metformin + betablock + ssri + statins + ppi + antipsycho + 
                      polypharmacy_cat +
                      hospinfect + hospgeneral,
                    include.miss = FALSE, chisq.test.perm = TRUE,
                    data = temp.data, method=2)


t = createTable(t1, hide.no = "no", digits = 1, show.n = T)


export2csv(t, file = '/proj/nobackup/simp2023007/wharf/baldanzi/baldanzi-simp2023007/atbgut/results/Table1_part2_simpler.csv', 
           which.table="descr", sep="\t", nmax = TRUE)

