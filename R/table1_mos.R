# Project Antibiotic and the gut microbiota 

# Gabriel Baldanzi 

# 2025-04-10

# Table 1. Population descriptive table - MOS

library(compareGroups)
library(data.table)
library(Hmisc)
library(sjmisc)

  setwd('/proj/nobackup/sens2019512/users/baldanzi/atb_gut/')

  mos <- fread("work/mos_working_dataset_revision.tsv", na.strings = c("","NA",NA))

  # Import models 
  load('work/mos_model_revision.Rdata')
  full.model <- full.model[-which(full.model =="extraction_plate")]
  
  cc <- complete.cases(mos[,basic.model,with=F])
  mos <- mos[cc, ]
  
  mos <- mos[last_EDATUM<Visit1-30 | is.na(last_EDATUM), ]
  
  mos[, atb_cat:=ifelse(N_all >=2 , 2, N_all )]
  mos[, atb_cat:=factor(atb_cat, c(0,1,2), c("none", "1 course",">=2 courses"))]

  mos[, sex := factor(sex, c(1,2), c("male","female") )]
  mos[, country_birth := factor(country_birth, c(0,1), c("no","Sweden") )]
  mos[, smoking := factor(smoking, c("never","former","current"), c("Never","Former","Current"))]
  mos[, diab_diag := factor(diab_diag,  c("no","yes") )]
  mos[, education:=factor(education, c("Compulsory", "Upper secondary", "University"))]
  mos[, polypharmacy_cat := fifelse(polypharmacy_cat == "5 or more", "yes", "no")]
  mos[, CCI2 := fifelse(CCIw >= 2, "yes", "no")]
  
  myvars <- c("age" , "sex", "smoking", 
              "BMI", "education",
              "country_birth", "diab_diag", "BP_diag",
              "atb_cat", "N_all",
              "rheumatic","cancer", "ppi",
              "CCI2", "metformin", "betablock", "ssri", "statins", "antipsycho",
              "polypharmacy_cat")
  
  mylabels <- c("Age (years)", "Sex", "Smoking",  
                "BMI (kg/m2)",  "Education", 
                "Born in Swedn", "Diabetes", "Hypertension", 
                "Number of ATB courses", "Number of antibiotic courses",
                "Rheumatic disease","Previous cancer", "PPI",
                "Charlson index >=2", "Metformin", "Beta-blocker", "SSRI", "Statins", "Antipsychotics", 
                "Polypharmacy")


  temp.data <- mos[,myvars,with=F]

  temp.data[!is.na(cancer),cancer:=ifelse(cancer=="no","no","yes")]
  temp.data[,Female:=ifelse(sex=="female","yes","no")]
  temp.data[,ppi:=ifelse(ppi=="1","yes","no")]
  temp.data[is.na(ppi), ppi:="no"]
  temp.data[,metformin:=ifelse(metformin=="1","yes","no")]
  temp.data[is.na(metformin), metformin:="no"]
  temp.data[,betablock:=ifelse(betablock=="1","yes","no")]
  temp.data[is.na(betablock), betablock:="no"]
  temp.data[,ssri:=ifelse(ssri=="1","yes","no")]
  temp.data[is.na(ssri), ssri:="no"]
  temp.data[,statins:=ifelse(statins=="1","yes","no")]
  temp.data[is.na(statins), statins:="no"]


  message("Assign labels")

# Variables for table 1 



j=1
for(i in myvars){
  label(temp.data[[i]]) <- mylabels[j]
  j=j+1
}

# Table 1 part 1 ####

t1 <- compareGroups(atb_cat ~  age + Female + smoking + education + country_birth,
                    include.miss = FALSE, chisq.test.perm = TRUE,
                    data = temp.data, method=2)


t = createTable(t1, hide.no = "no", 
                digits = 1, show.n = T)

export2csv(t, file = 'results/Table1_part1_mos.csv', 
           which.table="descr", sep="\t", nmax = TRUE)

# Table 1 part 2 ####
  full.model <- full.model[-which(full.model == "CCIw")]
  cc <- complete.cases(temp.data[, c(full.model, "CCI2"),with=F])
  temp.data <- temp.data[cc, ]
  
  t1 <- compareGroups(atb_cat ~ 
                        BMI + CCI2 + metformin + betablock + ssri + statins + ppi + antipsycho + 
                        polypharmacy_cat,
                      include.miss = FALSE, chisq.test.perm = TRUE,
                      data = temp.data, method=2)
  
  
  t = createTable(t1, hide.no = "no", 
                  digits = 1, show.n = T)
  
  export2csv(t, file = 'results/Table1_part2_mos.csv', 
             which.table="descr", sep="\t", nmax = TRUE)
  

