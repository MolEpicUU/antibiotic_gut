# Project Antibiotic and the gut microbiota 

# Gabriel Baldanzi 

# 2025-04-10

# Table 1. Population descriptive table - SCAPIS

  library(compareGroups)
  library(data.table)
  library(Hmisc)
  library(sjmisc)

  setwd('/proj/nobackup/sens2019512/users/baldanzi/atb_gut/')

  scapis <- fread("work/scapis_working_dataset_revision.tsv", na.strings = c("","NA",NA))

  # Import models 
  load('work/scapis_model_revision.Rdata')
  full.model <- full.model[-which(full.model =="site_plate")]
  
  cc <- complete.cases(scapis[,c(basic.model,"Site"),with=F])
  scapis <- scapis[cc, ]
  
  
  scapis[, atb_cat:=ifelse(N_all >=2 , 2, N_all )]
  scapis[, atb_cat:=factor(atb_cat, c(0,1,2), c("none", "1 course",">=2 courses"))]

  scapis[, Sex := factor(Sex, c("male","female") )]
  scapis[, placebirth := ifelse(placebirth=="Scandinavia", "yes", "no")]
  scapis[, smokestatus := factor(smokestatus, c("never","former","current"), c("Never","Former","Current"))]
  scapis[, diabd := factor(diabd,  c("no","yes") )]
  scapis[, education:=factor(education, c("Compulsory", "Upper secondary", "University"))]
  scapis[, polypharmacy12m_cat := fifelse(polypharmacy12m_cat == "5 or more", "yes", "no")]
  scapis[, CCI2 := fifelse(CCIw >= 2, "yes", "no")]

  myvars <- c("age" , "Sex",
            "smokestatus", "BMI", "education",
            "placebirth", "diabd", "Site",
            "atb_cat", "N_all","rheumatic","cancer",
            "hospgeneral", "hospinfect", "ppi",
            "CCI2", "metformin", "betablock", "ssri", "statins", "antipsycho", 
            "polypharmacy12m_cat")

  mylabels <- c("Age (years)", "Sex", "Smoking",  
              "BMI (kg/m2)",  "Education", 
              "Born in Scandinavia", "Diabetes", "Site",
              "Number of ATB courses", "Number of antibiotic courses",
              "Rheumatic disease","Previous cancer", 
              "Any hospitalization", "Hospitalization for infection", "PPI",
              "Charlson index >=2", "Metformin", "Beta-blocker", "SSRI", "Statins", "Antipsychotics", 
              "Polypharmacy")

  temp.data <- scapis[,myvars,with=F]

  temp.data[!is.na(cancer),cancer:=ifelse(cancer=="no","no","yes")]
  temp.data[,Female:=ifelse(Sex=="female","yes","no")]



  message("Assign labels")

# Variables for table 1 

j=1
for(i in myvars){
  label(temp.data[[i]]) <- mylabels[j]
  j=j+1
}

# Table 1 part 1 ####

t1 <- compareGroups(atb_cat ~ age + Female + smokestatus + education + placebirth,
                    include.miss = FALSE, chisq.test.perm = TRUE,
                    data = temp.data, method=2)


t = createTable(t1, hide.no = "no", 
                digits = 1, show.n = T)

export2csv(t, file = 'results/Table1_part1_scapis.csv', 
           which.table="descr", sep="\t", nmax = TRUE)

# Table 1 part 2 ####
  full.model <- full.model[-which(full.model == "CCIw")]
  cc <- complete.cases(temp.data[, c(full.model, "CCI2"),with=F])
  temp.data <- temp.data[cc, ]
  
  t1 <- compareGroups(atb_cat ~ 
                        BMI + CCI2 + metformin + betablock + ssri + statins + ppi + antipsycho + 
                      polypharmacy12m_cat + hospinfect + hospgeneral ,
                      include.miss = FALSE, chisq.test.perm = TRUE,
                      data = temp.data, method=2)
  
  
  t = createTable(t1, hide.no = "no", 
                  digits = 1, show.n = T)
  
  export2csv(t, file = 'results/Table1_part2_scapis.csv', 
             which.table="descr", sep="\t", nmax = TRUE)
  

