# Project antibiotic and gut microbiome 

# wrapper functions for stratified 

setwd('users/baldanzi/')

classes <- c("Class_Peni_Ext",  
             "Class_Peni_BetaS",    
            "Class_lincosamides", 
            "Class_SMZTMP",     
            "Class_FQs" , 
            "Class_Peni_Comb",     
            "Class_NIT",
            "Class_Peni_BetaR",    
            "Class_TCLs", 
            "Class_macrolides",     
            "Class_cephalosporins")



# Sex ####

interaction <- "sex"

for(abx in classes){
  
  jobname <- paste0("scapis_interaction_", interaction, "_", abx)
  out_file <- paste0("out/scapis_interaction_", interaction, "_", abx, ".out")
  sub <- system(paste("sbatch -A projname -n 16 -t 2:00:00 -J", jobname, 
               "-o" , out_file, 
               "--wrap='ml R_packages; Rscript revision__10.8_scapis_species_interaction.R ", 
               abx, interaction, "'", 
             "--mail-user=gabriel.baldanzi@medsci.uu.se  --mail-type=FAIL"), intern=T)
  cat(sub)
  
}


# Age  ####
interaction <- "age"

for(abx in classes){
  
  jobname <- paste0("scapis_interaction_", interaction, "_", abx)
  out_file <- paste0("out/scapis_interaction_", interaction, "_", abx, ".out")
  sub <- system(paste("sbatch -A projname -n 16 -t 2:00:00 -J", jobname, 
                      "-o" , out_file, 
                      "--wrap='ml R_packages; Rscript revision__10.8_scapis_species_interaction.R ", 
                      abx, interaction, "'", 
                      "--mail-user=gabriel.baldanzi@medsci.uu.se  --mail-type=FAIL"), intern=T)
  
}

