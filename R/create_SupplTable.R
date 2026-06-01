# Project antibiotic use and gut microbiota 

# Supplementary tables 

# This is a long script that combines .tsv files into a single large xlsx file. Some very large .tsv were later added manually.

  
  library(xlsx)
  library(tidyverse)
  library(data.table)
  library(stringr)
  message("finished loading pckgs")
  rm(list=ls())
  t0 <- Sys.time()
  print(t0)
  
  setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Error_SuppTable_20260504/')

# Exposures 
  
  exp <- c( 'Class_Peni_BetaS_1yr', 'Class_Peni_BetaS_1_4yr', 'Class_Peni_BetaS_4_8yr',
            'Class_Peni_BetaR_1yr', 'Class_Peni_BetaR_1_4yr', 'Class_Peni_BetaR_4_8yr', 
            'Class_Peni_Ext_1yr', 'Class_Peni_Ext_1_4yr', 'Class_Peni_Ext_4_8yr',
            'Class_Peni_Comb_1yr', 'Class_Peni_Comb_1_4yr', 'Class_Peni_Comb_4_8yr',
            'Class_cephalosporins_1yr', 'Class_cephalosporins_1_4yr', 'Class_cephalosporins_4_8yr',
            'Class_macrolides_1yr', 'Class_macrolides_1_4yr', 'Class_macrolides_4_8yr',
            'Class_lincosamides_1yr', 'Class_lincosamides_1_4yr', 'Class_lincosamides_4_8yr',
            'Class_TCLs_1yr', 'Class_TCLs_1_4yr', 'Class_TCLs_4_8yr',
            'Class_FQs_1yr', 'Class_FQs_1_4yr', 'Class_FQs_4_8yr',
            'Class_SMZTMP_1yr', 'Class_SMZTMP_1_4yr', 'Class_SMZTMP_4_8yr', 
            'Class_NIT_1yr', 'Class_NIT_1_4yr', 'Class_NIT_4_8yr')
  
  expafter <- c('Class_Peni_BetaS_after1yr',
           'Class_Peni_BetaR_after1yr', 
           'Class_Peni_Ext_after1yr',
           'Class_Peni_Comb_after1yr',
           'Class_cephalosporins_after1yr',
           'Class_macrolides_after1yr',
           'Class_lincosamides_after1yr',
           'Class_TCLs_after1yr',
           'Class_FQs_after1yr',
           'Class_SMZTMP_after1yr', 
           'Class_NIT_after1yr')
  
  exp_single <- c('classClass_Peni_BetaS1_4yr', 'classClass_Peni_BetaS4_8yr',
           'classClass_Peni_BetaR1_4yr', 'classClass_Peni_BetaR4_8yr', 
           'classClass_Peni_Ext1_4yr', 'classClass_Peni_Ext4_8yr',
            'classClass_Peni_Comb1_4yr', 'classClass_Peni_Comb4_8yr',
           'classClass_cephalosporins1_4yr', 'classClass_cephalosporins4_8yr',
           'classClass_macrolides1_4yr', 'classClass_macrolides4_8yr',
           'classClass_lincosamides1_4yr', 'classClass_lincosamides4_8yr',
           'classClass_TCLs1_4yr', 'classClass_TCLs4_8yr',
           'classClass_FQs1_4yr', 'classClass_FQs4_8yr',
            'classClass_SMZTMP1_4yr', 'classClass_SMZTMP4_8yr', 
           'classClass_NIT1_4yr', 'classClass_NIT4_8yr')
  
  atbnames_single <- rep(c('Penicillin V', 
                    'Flucloxacillin',
                    'Penicillins extended spectrum',
                    'Amoxicillin-clavulanic acid',
                    'Cephalosporins',
                    'Macrolides',
                    'Clindamycin',
                    'Tetracyclines',
                    'Fluoroquinolones',
                    'Sulfamethoxazole-trimethoprim',
                    'Nitrofurantoin'), each = 2)
  
  atbperiods_single <- rep(c("<4 years", "4-8 years"),11 )
  
  exp_doseresponse <- paste0(rep(paste0(c('Class_Peni_BetaS_1_4yr', 'Class_Peni_BetaS_4_8yr',
           'Class_Peni_BetaR_1_4yr', 'Class_Peni_BetaR_4_8yr', 
           'Class_lincosamides_1_4yr', 'Class_lincosamides_4_8yr',
           'Class_TCLs_1_4yr', 'Class_TCLs_4_8yr',
           'Class_FQs_1_4yr', 'Class_FQs_4_8yr'), "_cat"),each = 2),c("1","2"))
  
  # Antibiotics 
  atbnames <- rep(c('Penicillin V', 
                'Flucloxacillin',
                'Penicillins extended spectrum',
                'Amoxicillin-clavulanic acid',
                'Cephalosporins',
                'Macrolides',
                'Clindamycin',
                'Tetracyclines',
                'Fluoroquinolones',
                'Sulfamethoxazole-trimethoprim',
                'Nitrofurantoin'), each = 3)
  


# periods 
  atbperiods <- rep(c("<1 year", "1-4 years", "4-8 years"),11 )
  keyatb <- data.frame(antibiotic = rep(atbnames,2), period = rep(atbperiods,2), model = rep(c("basic model","full model"), each=33))

# Function to clean exposure, antibiotics, and periods 
  
  clean_exposure_fun <- function(dat, remove = F, exposure = "exposure"){
    dat$antibiotic <- factor(dat[[exposure]], exp, atbnames)
    dat$period <- factor(dat[[exposure]], exp, atbperiods)
    if(remove) {dat$exposure <- NULL}
    dat$model <- factor(dat$model, c("basic.model","full.model"), 
                        c("basic model","full model"))
    dat <- merge(keyatb, dat, by=c("model","antibiotic","period"))
    if("q.value" %in% names(dat)) {
      setcolorder(dat, c("antibiotic","period","outcome","model","beta","SE","p.value","q.value"))
    } else {setcolorder(dat, c("antibiotic","period","outcome","model","beta","SE","p.value"))}
    return(dat)
  }
  
  
  clean_exposureafter_fun <- function(dat, remove = F, exposure = "exposure"){
    
    dat$antibiotic <- factor(dat[[exposure]], expafter, c('Penicillin V', 
                                                          'Flucloxacillin',
                                                          'Penicillins extended spectrum',
                                                          'Amoxicillin-clavulanic acid',
                                                          'Cephalosporins',
                                                          'Macrolides',
                                                          'Clindamycin',
                                                          'Tetracyclines',
                                                          'Fluoroquinolones',
                                                          'Sulfamethoxazole-trimethoprim',
                                                          'Nitrofurantoin'))
    dat$period <- "1 year after"
    if(remove) {dat$exposure <- NULL}
    dat$model <- factor(dat$model, c("basic.model","full.model"), c("basic model","full model"))
    
    if("q.value" %in% names(dat)) {
      setcolorder(dat, c("antibiotic","period","outcome","model","beta","SE","p.value","q.value"))
    } else {setcolorder(dat, c("antibiotic","period","outcome","model","beta","SE","p.value"))}
    return(dat)
  }
  
  
  clean_exposuredoseresponse_fun <- function(dat, remove = F, exposure = "exposure"){
    dat$antibiotic <- factor(dat[[exposure]], exp_doseresponse, paste(rep(c('Penicillin V', 
                                                                            'Flucloxacillin',
                                                                            'Clindamycin',
                                                                            'Tetracyclines',
                                                                            'Fluoroquinolones'), each = 4), c("1 course", "2 courses")))
    dat$period <- NA
    dat$period[grep("1_4yr", dat[[exposure]] )] <- "1-4 years" 
    dat$period[grep("4_8yr", dat[[exposure]] )] <- "4-8 years"
    
    if(remove) {dat$exposure <- NULL}
    dat$model <- factor(dat$model, c("basic.model","full.model"), c("basic model","full model"))
    #dat <- merge(keyatb, dat, by=c("model","antibiotic","period"))
    if("q.value" %in% names(dat)) {
      setcolorder(dat, c("antibiotic","period","outcome","model","beta","SE","p.value","q.value"))
    } else {setcolorder(dat, c("antibiotic","period","outcome","model","beta","SE","p.value"))}
    return(dat)
  }
  
  clean_exposure_singledose_fun <- function(dat, remove = F, exposure = "exposure"){
    dat$antibiotic <- factor(dat[[exposure]], exp_single, atbnames_single)
    dat$period <- factor(dat[[exposure]], exp_single, atbperiods_single)
    if("N.y" %in% names(dat)) {
      dat <- dplyr::rename(dat, N = N.x, Nexposed = N.y)
    }
    if(remove) {dat$exposure <- NULL}
    dat$model <- factor(dat$model, c("basic.model","full.model"), c("basic model","full model"))
    if("q.value" %in% names(dat)) {
      setcolorder(dat, c("antibiotic","period","outcome","model","beta","SE","p.value","q.value"))
    } else {setcolorder(dat, c("antibiotic","period","outcome","beta","SE","p.value"))}
    return(dat)
  }

  
  clean_estimates <- function(dat){
    dat$beta <- round(dat$beta,3)
    dat$SE <- round(dat$SE,3)
    if("lci" %in% colnames(dat)) dat <- dplyr::rename(dat, LCI = lci, HCI = hci)
    dat$LCI <- NULL
    dat$HCI <- NULL
    dat <- dplyr::rename(dat, "p-value" = p.value)
    dat <- dplyr::rename(dat, "s.e." = SE)
    if("q.value" %in% names(dat)) {dat <- dplyr::rename(dat, "q-value" = q.value)}
    if(unique(dat$cohort) == "Meta-analysis")     {
      dat <- mutate(dat, Q = round(Q,2)) %>% dplyr::rename("Cochran's Q" = Q, "Cochran's Q p-value" = Qpval) 
      } else if ("Qpval" %in% names(dat)){
      dat$Qpval <- dat$Q <- dat$I2 <- NULL
      }
    if(unique(dat$cohort) == "Meta-analysis" & "GVIF" %in% names(dat))     { dat$N <- dat$GVIF <-  NULL}
    if(unique(dat$cohort) == "Meta-analysis" & "gvif" %in% names(dat))     { dat$N <- dat$gvif <-  NULL}
    if("gvif" %in% names(dat))    { dat <- dplyr::rename(dat, "GVIF" = gvif)}
    if("GVIF" %in% names(dat))    { dat$GVIF <- round(dat$GVIF, 2) }
    if("t.value" %in% names(dat)) {dat$t.value <- NULL}
    if("df" %in% names(dat)) {dat$df <- NULL}
    if('dfbeta' %in% colnames(dat)) {
      dat$dfbeta <- round(dat$dfbeta, 3)
      dat$beta_infl <- round(dat$beta_infl, 3)
      dat$se_infl <- round(dat$se_infl, 3)
    }
    return(dat)
  }

  rename_all <- function(dat){
    coh <- unique(dat$cohort)
    dat$cohort <- NULL
    if("message" %in% colnames(dat)) dat$message <- NULL
    nn <- colnames(dat)[-which(colnames(dat) %in% c("antibiotic","exposure","period","model","outcome"))]
    setnames(dat, nn, paste(nn, coh, sep="_"))
    if(paste0("message_", coh) %in% colnames(dat)) dat[[paste0("message_", coh)]] <- NULL
    return(dat)
  }

  merge_df_list <- function(df_list) {
    
    df_list <- lapply(df_list, function(x){
      if("exposure" %in% colnames(x)) x$exposure <- NULL
      return(x)})
    
    i <- length(df_list)
    
    final_df <- merge(df_list[[1]],df_list[[2]], by = c("antibiotic","period","outcome", "model"), all=T)
      for(j in 3:i){
        final_df <- merge(final_df, df_list[[j]], by = c("antibiotic","period","outcome","model"), all.x=T)
      }
    return(final_df)
  }
  
  
  
  rename_all_timewindown <- function(dat){
    coh <- unique(dat$cohort)
    dat$cohort <- NULL
    nn <- colnames(dat)[-which(colnames(dat) %in% c("antibiotic","exposure","period","model","outcome", "time_windown"))]
    setnames(dat, nn, paste(nn, coh, sep="_"))
    return(dat)
  }
  
  merge_df_list_timewindown <- function(df_list) {
    
    df_list <- lapply(df_list, function(x){
      if("exposure" %in% colnames(x)) x$exposure <- NULL
      return(x)})
    
    i <- length(df_list)
    
    final_df <- merge(df_list[[1]],df_list[[2]], by = c("antibiotic","period","outcome", "time_windown", "model"), all.x=T)
    for(j in 3:i){
      final_df <- merge(final_df, df_list[[j]], by = c("antibiotic","period","outcome","time_windown", "model"), all.x=T)
    }
    return(final_df)
  }
  
  
  
  taxa = fread('/Users/gabba126/Documents/PhD_projects/Microbiome/Taxonomy/Taxonomy_CHAMP.tsv')
  
  taxa_merge <- function(dat, full = T){
    if(full){
      dat <- merge(dat, taxa, by.x="outcome", by.y="MGS")
    dat$outcome <- dat$MainTax
    dat$MainTax <- NULL
    dat$Level <- NULL
    dat$species <- NULL
    } else {dat <- merge(dat, taxa[,.(MGS,MainTax)], by.x="outcome", by.y="MGS")

    dat$outcome <- dat$MainTax
    dat$MainTax <- NULL}
    return(dat)
  }
  
  clean_alpha <- function(dat) {
    dat$outcome <- factor(dat$outcome, c("shannon","richness","invsimpson"), c("shannon","richness","inverse simpson"))
    if("message" %in% colnames(dat)) dat$message <- NULL
    return(dat)
  }
  
  separate_sa <- function(dat, strata = c("male", "female")){
    dat[, sa := model]
    dat[, model := "full.model"]
    dat1 <- dat[sa == strata[1], ]
    dat2 <- dat[sa == strata[2], ]
    
    list_dat <- list(dat1, dat2)
    names(list_dat) <- strata
    return(list_dat)
  }
  
  
  
  merge_pairs_sa <- function(dat, hetero = NULL){
    
    
   dat <- lapply(1:length(dat) , function(i) dat[[i]][, strata := names(dat[i])])
    
    
    dat_meta <- rbind(dat[[1]], dat[[2]])
    setDT(dat_meta)
    dat_meta$exposure <- NULL
    value_vars <- setdiff(names(dat_meta), c("antibiotic", "period", "outcome" ,"model", "strata"))
    dat_meta <- dat_meta %>% pivot_wider(id_cols = c(antibiotic, period, outcome, model), names_from = "strata", values_from = all_of(value_vars), names_vary = "slowest")
    setDT(dat_meta)

    
    dat_scapis <- rbind(dat[[3]], dat[[4]])
    setDT(dat_scapis)
    dat_scapis$exposure <- NULL
    value_vars <- setdiff(names(dat_scapis), c("antibiotic", "period", "outcome" ,"model", "strata"))
    dat_scapis <- dat_scapis %>% pivot_wider(id_cols = c(antibiotic, period, outcome, model), names_from = "strata", values_from = all_of(value_vars), names_vary = "slowest")
    setDT(dat_scapis)

    dat_mos <- rbind(dat[[5]], dat[[6]])
    setDT(dat_mos)
    dat_mos$exposure <- NULL
    value_vars <- setdiff(names(dat_mos), c("antibiotic", "period", "outcome" ,"model", "strata"))
    dat_mos <- dat_mos %>% pivot_wider(id_cols = c(antibiotic, period, outcome, model), names_from = "strata", values_from = all_of(value_vars), names_vary = "slowest")
    setDT(dat_mos)

    
    if(length(dat)==7) { dat_simpler <- dat[[7]]
    } else {
    dat_simpler <- rbind(dat[[7]], dat[[8]])
    }
    setDT(dat_simpler)
    dat_simpler$exposure <- NULL
    
    
    value_vars <- setdiff(names(dat_simpler), c("antibiotic", "period", "outcome" ,"model", "strata"))
    dat_simpler <- dat_simpler %>% pivot_wider(id_cols = c(antibiotic, period, outcome, model), names_from = "strata", values_from = all_of(value_vars), names_vary = "slowest")
    setDT(dat_simpler)
    
    return(list(dat_meta, dat_scapis, dat_simpler, dat_mos))
    
  }
  

  
  remove_obj <- function(){
    list_classes <- do.call(c,lapply(ls(), function(x) "data.table" %in% class(get(x)) ))
    rm(list=ls()[list_classes])
  }
  


# Supplemental Table 1 - alpha --------------------------------------------------------

  alpha_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha.tsv')
  alpha_meta <- alpha_meta[!exposure %in% c("N1yr","N1_4yr","N4_8yr"), ]

  df_list <- lapply(c("Meta-analysis","SCAPIS","SIMPLER","MOS"), function(co) alpha_meta[cohort==co,]) 
  
  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)


  
  mer_alpha <- merge_df_list(df_list)

  
  write.xlsx2(mer_alpha[order(mer_alpha$model, mer_alpha$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 1", col.names=T, row.names=F, append=F)

  # Supplemental Table 2 - sex stratified  and alpha  -------------------------------------------------------
  alpha_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha_sa__age_sex.tsv')
  list_meta_sex <- separate_sa(alpha_meta, strata = c("male", "female"))
  

  scapis  <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_scapis_alpha_sensitivityanalyses_age_sex.tsv')
  simpler <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_simpler_alpha_sensitivityanalyses_age_sex.tsv')
  mos     <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_mos_alpha_sensitivityanalyses_age_sex.tsv')
  
  list_scapis_sex <- separate_sa(scapis, strata = c("male", "female"))
  list_simpler_sex <- separate_sa(simpler, strata = c("male", "female"))
  list_mos_sex <- separate_sa(mos, strata = c("male", "female"))

  
  df_list <- c(list_meta_sex, list_scapis_sex, list_mos_sex,  list_simpler_sex)

  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  df_list <- lapply(df_list, function(dat) {
    setDT(dat)
    dat[, model := "full model"]
  })
  
  
  df_list <- merge_pairs_sa(df_list)
  mer_alpha <- merge_df_list(df_list)

  
  # Interaction term 
  alpha_int <- fread("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha__interaction.tsv", na.strings = c("NA", "", NA))
  alpha_int[, c("reference", "exposure") := tstrsplit(exposure, ":", fixed = T)]
  alpha_int_sex <- alpha_int[model == "sex"]
  alpha_int_age <- alpha_int[model == "age"]
  
  alpha_int_sex[, model := "full.model"]
  alpha_int_sex <- clean_alpha(alpha_int_sex)
  alpha_int_sex <- clean_exposure_fun(alpha_int_sex, remove = T)
  setnames(alpha_int_sex, c("beta","SE", "p.value"), paste0("interaction_", c("beta","SE", "p.value")))
  mer_alpha <- merge(mer_alpha, alpha_int_sex[, c("antibiotic", "period", "outcome", "interaction_beta", "interaction_SE", "interaction_p.value", "p_omnibus", "q_omnibus")], by = c("antibiotic", "period", "outcome"))
  
  setcolorder(mer_alpha, c("antibiotic", "period", "outcome", "model", grep("Meta-analysis", colnames(mer_alpha), v=T), 
                           "interaction_beta", "interaction_SE", "interaction_p.value", "p_omnibus", "q_omnibus"))
  
  to_remove <- grep("sa_" ,colnames(mer_alpha), v=T)
  mer_alpha <- mer_alpha[, !c(to_remove), with=F]
  
  write.xlsx2(mer_alpha[order(mer_alpha$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 2", col.names=T, row.names=F, append=T)
  
  
  # Supplemental Table 3 - age stratified  and alpha  -------------------------------------------------------
  alpha_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha_sa__age_sex.tsv')
  list_meta_age <- separate_sa(alpha_meta, strata = c("age <= 55 years", "age > 55 years"))
  
  
  scapis  <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_scapis_alpha_sensitivityanalyses_age_sex.tsv')
  simpler <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_simpler_alpha_sensitivityanalyses_age_sex.tsv')
  mos     <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_mos_alpha_sensitivityanalyses_age_sex.tsv')
  
  list_scapis_age <- separate_sa(scapis, strata = c("age <= 55 years", "age > 55 years"))
  list_simpler_age <- separate_sa(simpler, strata = c("age > 55 years"))
  list_simpler_age[[2]] <- NULL
  list_mos_age <- separate_sa(mos, strata = c("age <= 55 years", "age > 55 years"))
  
  
  df_list <- c(list_meta_age, list_scapis_age, list_mos_age, list_simpler_age)

  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  df_list <- lapply(df_list, function(dat) {
    setDT(dat)
    dat[, model := "full model"]
  })
  
  

  df_list <- merge_pairs_sa(df_list)
  mer_alpha <- merge_df_list(df_list)
  
  
  # Interaction term 
  alpha_int_age[, model := "full.model"]
  alpha_int_age <- clean_alpha(alpha_int_age)
  alpha_int_age <- clean_exposure_fun(alpha_int_age, remove = T)
  setnames(alpha_int_age, c("beta","SE", "p.value"), paste0("interaction_", c("beta","SE", "p.value")))
  mer_alpha <- merge(mer_alpha, alpha_int_age[, c("antibiotic", "period", "outcome", "interaction_beta", "interaction_SE", "interaction_p.value", "p_omnibus", "q_omnibus")], by = c("antibiotic", "period", "outcome"))
  
  setcolorder(mer_alpha, c("antibiotic", "period", "outcome", "model", grep("Meta-analysis", colnames(mer_alpha), v=T), 
                           "interaction_beta", "interaction_SE", "interaction_p.value", "p_omnibus", "q_omnibus"))
  
  to_remove <- grep("sa_" ,colnames(mer_alpha), v=T)
  mer_alpha <- mer_alpha[, !c(to_remove), with=F]
  
  write.xlsx2(mer_alpha[order(mer_alpha$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 3", col.names=T, row.names=F, append=T)
  
  
  # Supplemental Table 4 - time windowns and alpha  -------------------------------------------------------
  
  meta <- fread("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha_sa_timewindows.tsv")
  meta[, time_windown := model]
  meta[, model := "full.model"]
  scapis <- fread("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_scapis_alpha_sensitivityanalyses_timewindows.tsv")
  simpler <- fread("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_simpler_alpha_sensitivityanalyses_timewindows.tsv")
  mos <- fread("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/res_mos_alpha_sensitivityanalyses_timewindows.tsv")
  
  
  df_list <- list(meta, scapis, simpler, mos)
  names(df_list) <- c("Meta-analysis","SCAPIS","SIMPLER","MOS")
  
  df_list <- lapply(df_list, function(dat) dat[, model := "full.model"])

  
  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all_timewindown)
  
  mer_alpha <- merge_df_list_timewindown(df_list)
  
  
  write.xlsx2(mer_alpha[order(mer_alpha$model, mer_alpha$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 4", col.names=T, row.names=F, append=T)
  
  # Supplemental Table 5 - negative exposure  and alpha  -------------------------------------------------------
  
  alpha_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha_negexposure.tsv')
  alpha_meta <- alpha_meta[!exposure %in% c("N1yr","N1_4yr","N4_8yr","Nafter1yr"), ]
  alpha_meta <- alpha_meta[grep("after",exposure),]
  alpha_meta[, q.value := p.adjust(p.value), by = c("outcome", "model", "cohort")]
  df_list <- lapply(c("Meta-analysis","SCAPIS","SIMPLER","MOS"), function(co) alpha_meta[cohort==co,]) 
  
  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_exposureafter_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  
  mer_alpha <- merge_df_list(df_list)
  mer_alpha$`N_Meta-analysis` <- NULL

  
  write.xlsx2(mer_alpha[order(mer_alpha$model, mer_alpha$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 5", col.names=T, row.names=F, append=T)
  
  # Supplemental Table 6 - negative exposure2  and alpha  -------------------------------------------------------
  
  alpha_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha_negexposure_preatbzero.tsv')

  alpha_meta <- alpha_meta[grep("after",exposure),]
  alpha_meta[, q.value := p.adjust(p.value), by = c("outcome", "model", "cohort")]
  

  df_list <- lapply(c("Meta-analysis","SCAPIS","SIMPLER","MOS"), function(co) alpha_meta[cohort==co,]) 
  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_exposureafter_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  
  mer_alpha <- merge_df_list(df_list)
  mer_alpha$`N_Meta-analysis` <- NULL
  
  
  write.xlsx2(mer_alpha[order(mer_alpha$model, mer_alpha$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 6", col.names=T, row.names=F, append=T)
  
  
  # Supplemental Table 7 - antibiotic type and alpha sensitivity analysis  --------------------------------------------------------
  
  res_alpha_sa <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha_sa.tsv')
  
  res_alpha_sa <- res_alpha_sa[!exposure %in% c("N1yr","N1_4yr", "N4_8yr"),]
  
  meta_sa1 <- res_alpha_sa[model == "full.model_hospitalizedinfect" & cohort == "Meta-analysis",]
  meta_sa1$model <- "full.model"
  
  meta_sa2 <- res_alpha_sa[model == "full.model_hospitalizedgeneral" & cohort == "Meta-analysis",]
  meta_sa2$model <- "full.model"
  
  scapis_sa1 <- res_alpha_sa[model == "full.model_hospitalizedinfect" & cohort == "SCAPIS",]
  scapis_sa1$model <- "full.model"
  scapis_sa1$cohort <- "SCAPIS_hospinfec"
  
  scapis_sa2 <- res_alpha_sa[model == "full.model_hospitalizedgeneral" & cohort == "SCAPIS",]
  scapis_sa2$model <- "full.model"
  scapis_sa2$cohort <- "SCAPIS_hospgen"
  
  simpler_sa1 <- res_alpha_sa[model == "full.model_hospitalizedinfect" & cohort == "SIMPLER",]
  simpler_sa1$model <- "full.model"
  simpler_sa1$cohort <- "SIMPLER_hospinfec"
  
  simpler_sa2 <- res_alpha_sa[model == "full.model_hospitalizedgeneral" & cohort == "SIMPLER",]
  simpler_sa2$model <- "full.model"
  simpler_sa2$cohort <- "SIMPLER_hospgen"
  
  df_list <- list(meta_sa1, meta_sa2, scapis_sa1, scapis_sa2, simpler_sa1, simpler_sa2)
  
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  
  mer_sa <- merge_df_list(df_list)
  mer_sa$GVIF_hospgen <- mer_sa$GVIF_hospinfec <- NULL
  mer_sa$GVIF_SCAPIS_hospinfec <- mer_sa$GVIF_SCAPIS_hospgen <- NULL
  mer_sa$GVIF_SIMPLER_hospinfec <- mer_sa$GVIF_SIMPLER_hospgen <- NULL
  mer_sa$N_hospgen <- mer_sa$N_hospinfec <- NULL
  
  write.xlsx2(mer_sa[order(mer_sa$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 7", col.names=T, row.names=F, append=T)
  
  
  # Supplemental Table 8 - singledose  and alpha  -------------------------------------------------------
  
  alpha_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_alpha_singledose.tsv')
  alpha_meta <- alpha_meta[exposure!="classother1_4yr",]
  scapis <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/scapis_alpha_singledose.tsv')
  simpler <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/simpler_alpha_singledose.tsv')
  simpler <- simpler[exposure!="classother1_4yr",]
  mos <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/mos_alpha_singledose.tsv')

  df_list <- list(alpha_meta, scapis, simpler, mos)  
  df_list <- lapply(df_list, clean_alpha)
  df_list <- lapply(df_list, clean_exposure_singledose_fun, T)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)

  mer_alpha <- merge_df_list(df_list)
  mer_alpha <- mer_alpha[order(mer_alpha$model),]
  mer_alpha <- mer_alpha[order(mer_alpha$model, mer_alpha$outcome),]
  mer_alpha$`war_Meta-analysis` <- NULL
  
  write.xlsx2(mer_alpha[order(mer_alpha$model, mer_alpha$outcome),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 8", col.names=T, row.names=F, append=T)
  
  remove_obj()
  

  
  # Supplemental Table 9 - antibiotic and species --------------------------------------------------------
  
  taxa = fread('/Users/gabba126/Documents/PhD_projects/Microbiome/Taxonomy/Taxonomy_CHAMP.tsv')
  
  # Species prevalence
  species_prevalence_scapis  <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/prevalence_species_scapis.tsv')
  species_prevalence_mos     <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/prevalence_species_mos.tsv')
  species_prevalence_simpler <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/prevalence_species_simpler.tsv')
  species_prevalence <- merge(species_prevalence_scapis, species_prevalence_simpler, by = "species")
  species_prevalence <- merge(species_prevalence, species_prevalence_mos, by = "species")
  species_prevalence[, Median := round(median(c(prevalence_scapis, prevalence_mos, prevalence_simpler)),1), by = species]
  species_prevalence[, Min := round(min(c(prevalence_scapis, prevalence_mos, prevalence_simpler)),1), by = species]
  species_prevalence[, Max := round(max(c(prevalence_scapis, prevalence_mos, prevalence_simpler)),1), by = species]
  species_prevalence[, prevalence := paste0(Median, " [",Min,"-",Max,"]")]
  
  
  N_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_species_class.tsv')
  metainflobs <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/metainflobs_byclass.tsv')
  metainflobs$cohort <- "Inf_obs"
  metainflobs$Q <- NULL
  metainflobs$Qpval <- NULL
  metainflobs$I2 <- NULL
  
  N_scapis <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/scapis_species_class.tsv')

  N_scapis[outcome %in% N_meta$outcome, q.value := p.adjust(p.value, method = "BH"), by = model]
  N_scapis$cohort <- "SCAPIS"
  
  N_simpler <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/simpler_species_class.tsv')

  N_simpler[outcome %in% N_meta$outcome, q.value := p.adjust(p.value, method = "BH"), by = model]
  N_simpler$cohort <- "SIMPLER"
  
  N_mos <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/mos_species_class.tsv')

  N_mos[outcome %in% N_meta$outcome, q.value := p.adjust(p.value, method = "BH"), by = model]
  N_mos$cohort <- "MOS"
  
  N_scapis_io <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/scapis_infobs_byclass.tsv')
  N_simpler_io <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/simpler_infobs_byclass.tsv')
  N_mos_io <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/mos_infobs_byclass.tsv')
  
  res_hetero <- N_meta[Qpval<0.05 & q.value<0.05 & model == "full.model",.(outcome,exposure,model)]
  metainflobs <- merge(res_hetero, metainflobs, by=c("outcome","exposure","model"))
  N_scapis_io <- merge(res_hetero, N_scapis_io[, .(outcome,exposure, model, beta_infl, se_infl, p.value_inf, dfbeta)], by=c("outcome","exposure","model"))
  N_simpler_io <- merge(res_hetero, N_simpler_io[, .(outcome,exposure, model, beta_infl, se_infl, p.value_inf, dfbeta)], by=c("outcome","exposure","model"))
  N_mos_io <- merge(res_hetero, N_mos_io[, .(outcome,exposure, model, beta_infl, se_infl, p.value_inf, dfbeta)], by=c("outcome","exposure","model"))
  
  N_scapis <- merge(N_scapis, N_scapis_io, by=c("outcome","exposure","model"), all.x=T)
  N_simpler <- merge(N_simpler, N_simpler_io, by=c("outcome","exposure","model"), all.x=T)
  N_mos <- merge(N_mos, N_mos_io, by=c("outcome","exposure","model"), all.x=T)
  
  df_list <- list(N_meta, metainflobs, N_scapis, N_simpler, N_mos)
  
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  mer_species <- merge_df_list(df_list)
  
  mer_species <- merge(mer_species, species_prevalence[, .(species, prevalence)], by.x="outcome", "species", all.x = T)
  
  mer_species <- taxa_merge(mer_species , full=T)
  
  mer_species$df_MOS <- NULL
  setcolorder(mer_species, c("outcome","prevalence", "antibiotic", "period","model"))
  
  length(unique(mer_species$outcome))
  
  fwrite(mer_species[order(mer_species$model),], file="suppl.table9.csv")
  remove_obj()
  
  
  # Supplemental Table 10 - sensitivity analysis species --------------------------------------------------------
  
  meta_clean <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_species_class_clean.tsv')
  meta_sa <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_byclass_sa12.tsv')
  meta_sa1 <- meta_sa[model == "class_sa1",]
  meta_sa1$model <- "full.model"
  meta_sa1$cohort <- "Meta-analysis"
  meta_sa2 <- meta_sa[model == "class_sa2",]
  meta_sa2$model <- "full.model"
  meta_sa2$cohort <- "Meta-analysis"
  
  scapis_sa1 <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/scapis_class_sa1.tsv')
  scapis_sa1$cohort <- "SCAPIS_sa1"
  scapis_sa2<- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/scapis_class_sa2.tsv')
  scapis_sa2$cohort <- "SCAPIS_sa2"
  simpler_sa1 <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/simpler_class_sa1.tsv')
  simpler_sa1$cohort <- "SIMPLER_sa1"
  simpler_sa2 <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/simpler_class_sa2.tsv')
  simpler_sa2$cohort <- "SIMPLER_sa2"
  
  df_list <- list(meta_sa1, meta_sa2,scapis_sa1, scapis_sa2, simpler_sa1, simpler_sa2)
  
  clean_q <- function(dat) {
    dat <- merge(meta_clean[model=="full.model" & q.value<.05,.(exposure, outcome)], dat, by=c("exposure","outcome"))
    
    return(dat)
  }
  
  df_list <- lapply(df_list, clean_q)
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  
  mer_sa <- merge_df_list(df_list)
  mer_sa$model <- NULL
  mer_sa <- taxa_merge(mer_sa , full=F)
  
  write.xlsx2(mer_sa, "Supp.tables_temp.xlsx", sheetName="Suppl. Table 10", col.names=T, row.names=F, append=T)
  
  remove_obj()
  
  # Supplemental Table 11 - Single dose and species  --------------------------------------------------------

  
  N_meta <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/meta_species_singledose.tsv')
  
  N_scapis <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/scapis_species_singledose.tsv')
  N_scapis <-N_scapis[outcome %in% N_meta$outcome]
  N_scapis[, q.value := p.adjust(p.value, method = "BH"), by = model]
  
  N_simpler <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/simpler_species_singledose.tsv')
  N_simpler <-N_simpler[outcome %in% N_meta$outcome]
  N_simpler[, q.value := p.adjust(p.value, method = "BH"), by = model]
  
  N_mos <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/mos_species_singledose.tsv')
  N_mos <- N_mos[outcome %in% N_meta$outcome]
  N_mos[, q.value := p.adjust(p.value, method = "BH"), by = model]
  
  df_list <- list(N_meta, N_scapis, N_simpler, N_mos)
  df_list <- lapply(df_list, clean_exposure_singledose_fun, T)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  mer_species <- merge_df_list(df_list)
  mer_species <- taxa_merge(mer_species , full=F)
  
  mer_species$`war_Meta-analysis` <- NULL
  
  supptable13 <- mer_species[order(mer_species$model),]
  
  fwrite(mer_species[order(mer_species$model),], file="suppltable_11.tsv", sep = '\t')
  # data.table file too large for R to handle direct into xlsx.
  #write.xlsx2(mer_species[order(mer_species$model),], "Supp.tables_temp.xlsx", sheetName="Suppl. Table 9", col.names=T, row.names=F, append=T)
  #remove_obj()

  
  remove_obj()
  
  # Supplemental Table 12 - sex-stratified  --------------------------------------------------------
  setwd("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/")
  meta_male <- fread('meta_species__sa_stratified_allantibiotics_male.tsv', na.strings = c("NA", "", NA))
  meta_female <- fread('meta_species__sa_stratified_allantibiotics_female.tsv', na.strings = c("NA", "", NA))
  stopifnot(all(meta_male$outcome %in% meta_female$outcome))
  
  list_meta_sex <- list(male =  meta_male, female = meta_female)
  
  scapis_male   <- fread('scapis_species_class_sa_stratified_allantibiotics__male.tsv', na.strings = c("NA", "", NA)  ) %>% mutate(cohort = "SCAPIS") %>% filter(outcome %in% meta_male$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  scapis_female   <- fread('scapis_species_class_sa_stratified_allantibiotics__female.tsv', na.strings = c("NA", "", NA)  ) %>% mutate(cohort = "SCAPIS") %>% filter(outcome %in% meta_male$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  simpler_male <- fread('simpler_species_class_sa_stratified_allantibiotics__male.tsv', na.strings = c("NA", "", NA) ) %>% mutate(cohort = "SIMPLER") %>% filter(outcome %in% meta_male$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  simpler_female <- fread('simpler_species_class_sa_stratified_allantibiotics__female.tsv', na.strings = c("NA", "", NA) )  %>% mutate(cohort = "SIMPLER") %>% filter(outcome %in% meta_male$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  mos_male <- fread('mos_species_class_sa_stratified_allantibiotics__male.tsv', na.strings = c("NA", "", NA) )  %>% mutate(cohort = "MOS") %>% filter(outcome %in% meta_male$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  mos_female <- fread('mos_species_class_sa_stratified_allantibiotics__female.tsv', na.strings = c("NA", "", NA) )  %>% mutate(cohort = "MOS") %>% filter(outcome %in% meta_male$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  list_scapis_sex <- list(male =  scapis_male, female = scapis_female)
  list_simpler_sex <- list(male =  simpler_male, female = simpler_female)
  list_mos_sex <- list(male =  mos_male, female = mos_female)
  
  
  df_list <- c(list_meta_sex, list_scapis_sex, list_mos_sex,  list_simpler_sex)
  df_list <- lapply(df_list, function(dat) {
    setDT(dat)
    dat[, model := "full.model"]
    setDF(dat)
  })

  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  df_list <- lapply(df_list, function(dat) {
    setDT(dat)
    dat[, model := "full model"]
  })
  
  
  df_list <- merge_pairs_sa(df_list)
  mer_sa <- merge_df_list(df_list)
  
  mer_sa <- taxa_merge(mer_sa , full=T)
  
  setDT(mer_sa)
  
  
  # Intercation term 
  spint <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/meta_species__sa_interaction.tsv', na.strings = c("NA", NA, ""))
  setnames(spint, "class", "antibiotic")
  main <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/meta_species_class_clean.tsv', na.strings = c("NA", NA, ""))
  main <- main[model=="full.model"]
  
  spint[, c("reference", "exposure") := tstrsplit(exposure, ":", fixed = T)]
  spint <- merge(main[, .(outcome, exposure, main = beta)], spint)
  
  # Remove results with influential observations in main analysis 
  to_keep <- spint[, .N, .(outcome, antibiotic, model )][N==3,]
  spint <- merge(spint, to_keep[, !c("N"), with=F], by = c("outcome", "antibiotic", "model"))
  
  # Recalculate q_omnibus 
  q_omnibus <- unique(spint[, .(antibiotic, outcome , model , p_omnibus)])
  q_omnibus[, q_omnibus := p.adjust(p_omnibus, "BH"), by = "model"]
  spint <- spint[, q_omnibus := NULL]
  spint <- merge(spint, q_omnibus, by = c("outcome", "antibiotic", "model","p_omnibus"))
  
  spint_sex <- spint[model == "sex"]
  spint_age <- spint[model == "age"]
  
  spint_sex[, model := "full.model"]
  spint_sex <- clean_exposure_fun(spint_sex, remove = T)
  spint_sex <- taxa_merge(spint_sex, full = F)
  setnames(spint_sex, c("beta","SE", "p.value"), paste0("interaction_", c("beta","SE", "p.value")))
  mer_sa <- merge(mer_sa, spint_sex[, c("antibiotic", "period", "outcome", "interaction_beta", "interaction_SE", "interaction_p.value",  "p_omnibus", "q_omnibus")], by = c("antibiotic", "period", "outcome"), all.x=T)
  
  setcolorder(mer_sa , neworder = c("outcome", "antibiotic", "period", "model" ,  "beta_Meta-analysis_male", "s.e._Meta-analysis_male",  "p-value_Meta-analysis_male", "q-value_Meta-analysis_male", "Cochran's Q_Meta-analysis_male", "Cochran's Q p-value_Meta-analysis_male", "I2_Meta-analysis_male",  
                                    "beta_Meta-analysis_female","s.e._Meta-analysis_female", "p-value_Meta-analysis_female","q-value_Meta-analysis_female",  "Cochran's Q_Meta-analysis_female" ,"Cochran's Q p-value_Meta-analysis_female" , "I2_Meta-analysis_female", 
                                    "interaction_beta", "interaction_SE", "interaction_p.value", "p_omnibus", "q_omnibus"))
  
  
  to_remove <- grep("Nexp_", colnames(mer_sa), value = T)
  to_remove
  mer_sa <- mer_sa[, !c(to_remove), with=F]
  to_remove <- grep("N_.*Meta-analysis_.*male", colnames(mer_sa), value = T)
  to_remove
  mer_sa <- mer_sa[, !c(to_remove), with=F]
  
  mer_sa$cohort <- NULL

  setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Error_SuppTable_20260504/')
  fwrite(mer_sa[order(mer_sa$model),], file="suppl.table12.csv")
  # data.table file too large for R to handle direct into xlsx.
  remove_obj()
  
  # Supplemental Table 13 - age-stratified  --------------------------------------------------------
  setwd("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results")
  
  
  meta_agebelow55 <- fread('meta_species__sa_stratified_allantibiotics_agebelow55.tsv', na.strings = c("NA", "", NA)) %>% mutate(model = "full.model")
  meta_ageabove55 <- fread('meta_species__sa_stratified_allantibiotics_ageabove55.tsv', na.strings = c("NA", "", NA)) %>% mutate(model = "full.model")
  
  list_meta_age <- list("age <= 55 years" =  meta_agebelow55, "age > 55 years" = meta_ageabove55)
  
  scapis_agebelow55   <- fread('scapis_species_class_sa_stratified_allantibiotics__agebelow55.tsv', na.strings = c("NA", "", NA)  ) %>% mutate(cohort = "SCAPIS", model = "full.model") %>% filter(outcome %in% meta_agebelow55$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  scapis_ageabove55   <- fread('scapis_species_class_sa_stratified_allantibiotics__ageabove55.tsv', na.strings = c("NA", "", NA)  ) %>% mutate(cohort = "SCAPIS", model = "full.model") %>% filter(outcome %in% meta_agebelow55$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))

  simpler_ageabove55 <- fread('simpler_species_class_sa_stratified_allantibiotics__ageabove55.tsv', na.strings = c("NA", "", NA) )  %>% mutate(cohort = "SIMPLER", model = "full.model") %>% filter(outcome %in% meta_agebelow55$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  mos_agebelow55 <- fread('mos_species_class_sa_stratified_allantibiotics__agebelow55.tsv', na.strings = c("NA", "", NA) )  %>% mutate(cohort = "MOS", model = "full.model") %>% filter(outcome %in% meta_agebelow55$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  mos_ageabove55 <- fread('mos_species_class_sa_stratified_allantibiotics__ageabove55.tsv', na.strings = c("NA", "", NA) )  %>% mutate(cohort = "MOS", model = "full.model") %>% filter(outcome %in% meta_agebelow55$outcome) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
  
  list_scapis_age <- list("age <= 55 years" =  scapis_agebelow55, "age > 55 years" = scapis_ageabove55)
  list_simpler_age <- list("age > 55 years" = simpler_ageabove55)
  list_mos_age <- list("age <= 55 years" =  mos_agebelow55, "age > 55 years" = mos_ageabove55)
  
  
  df_list <- c(list_meta_age, list_scapis_age, list_mos_age,  list_simpler_age)
  
  
  df_list <- lapply(df_list, clean_exposure_fun)
  df_list <- lapply(df_list, clean_estimates)
  df_list <- lapply(df_list, rename_all)
  df_list <- lapply(df_list, function(dat) {
    setDT(dat)
    dat[, model := "full model"]
  })
  
  
  df_list <- merge_pairs_sa(df_list)
  mer_sa <- merge_df_list(df_list)
  
  mer_sa <- taxa_merge(mer_sa , full=T)
  
  spint_age[, model := "full.model"]
  spint_age <- clean_exposure_fun(spint_age, remove = T)
  spint_age <- taxa_merge(spint_age, full = F)
  setnames(spint_age, c("beta","SE", "p.value"), paste0("interaction_", c("beta","SE", "p.value")))
  mer_sa <- merge(mer_sa, spint_age[, c("antibiotic", "period", "outcome", "interaction_beta", "interaction_SE", "interaction_p.value",  "p_omnibus", "q_omnibus")], by = c("antibiotic", "period", "outcome"), all.x=T)
  
  setcolorder(mer_sa , neworder = c("outcome", "antibiotic", "period", "model" ,  "beta_Meta-analysis_age <= 55 years", "s.e._Meta-analysis_age <= 55 years",  "p-value_Meta-analysis_age <= 55 years", "q-value_Meta-analysis_age <= 55 years", "Cochran's Q_Meta-analysis_age <= 55 years", "Cochran's Q p-value_Meta-analysis_age <= 55 years", "I2_Meta-analysis_age <= 55 years",  
                                    "beta_Meta-analysis_age > 55 years","s.e._Meta-analysis_age > 55 years", "p-value_Meta-analysis_age > 55 years","q-value_Meta-analysis_age > 55 years",  "Cochran's Q_Meta-analysis_age > 55 years" ,"Cochran's Q p-value_Meta-analysis_age > 55 years" , "I2_Meta-analysis_age > 55 years", 
                                    "interaction_beta", "interaction_SE", "interaction_p.value",  "p_omnibus", "q_omnibus"))
  
  
 
  
  to_remove <- grep("Nexp_", colnames(mer_sa), value = T)
  to_remove
  mer_sa <- mer_sa[, !c(to_remove), with=F]
  
  to_remove <- grep("N_.*Meta-analysis_.*age", colnames(mer_sa), value = T)
  to_remove
  mer_sa <- mer_sa[, !c(to_remove), with=F]
  
  
  setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__tables')
  fwrite(mer_sa[order(mer_sa$model),], file="suppl.table13.csv")
  # data.table file too large for R to handle direct into xlsx.
  remove_obj()
  
  
  
  
  # Supplemental Table 14 - Species  cardiometabolic markers --------------------------------------------------------
  
  
  scapis <- fread('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results/scapis_healthoutcomes_spearman.tsv')
  scapis <- scapis[outcome %in% c("BMI", "WHR", "SBP", "TG", "TC", "non-HDL", "HbA1c", "CRP"), ]
  scapis[, outcome := factor(outcome, c("BMI", "WHR", "SBP", "TG", "TC", "non-HDL", "HbA1c", "CRP"))]
  scapis <- scapis[order(outcome), ]
  
  scapis[, method := NULL]
  scapis[, covariates := NULL]
  
  setnames(scapis, c("species", "outcome"), c("outcome", "marker"))
  scapis <- taxa_merge(scapis , full=F)
  setnames(scapis, c("outcome", "p.value", "q.value"), c("species", "p-value", "q-value"))
  setcolorder(scapis, c("species", "marker", "rho", "p-value", "q-value", "N"))
  
  value_vars <- setdiff(colnames(scapis), c("species", "marker"))

  scapis <- scapis %>% tidyr::pivot_wider(id_cols = species, names_from = marker, values_from = all_of(value_vars), names_vary = "slowest")
  
  print(dim(scapis))
  setDT(scapis)
  write.xlsx2(scapis, "Supp.tables_temp.xlsx", sheetName="Suppl. Table 11", col.names=T, row.names=F, append=T)
  
  
  # Supplemental Table 15 --------------------------------------------------------
  
  # The suppl. table 15 is the list of hospital discharge diagnosis treated with 
  # antibiotics and that was manually pasted in the final document using the 
  # table submitted for the ethical application. 
  