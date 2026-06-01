# Project - Antibiotic use and the gut microbiota


# This script will plot the results from the full model analysis and create correlation plots
# with the results from sensitivity analysis 

  rm(list=ls())
  library(data.table)
  library(patchwork)
  library(ggpubr)
  library(tidyverse)
  
  setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results')
  
  res_main <- fread("meta_species_class_clean.tsv")
  res_main <- res_main[model== "full.model" & q.value<0.05,]

  res_scapissimpler <- res_main
  
  res_sa <- fread("meta_byclass_sa12.tsv")
  res_sa1 <- res_sa[model == "class_sa1",]
  res_sa2 <- res_sa[model == "class_sa2",]
  
  
  res_main <- res_main %>% separate(exposure, into =c("antibiotic", "period"), sep = "_(?=\\d)", extra = "merge", remove = F)
  res_scapissimpler <- res_scapissimpler %>% separate(exposure, into =c("antibiotic", "period"), sep = "_(?=\\d)", extra = "merge", remove = F)
  res_sa1 <- res_sa1 %>% separate(exposure, into =c("antibiotic", "period"), sep = "_(?=\\d)", extra = "merge", remove = F)
  res_sa2 <- res_sa2 %>% separate(exposure, into =c("antibiotic", "period"), sep = "_(?=\\d)", extra = "merge", remove = F)
  
  setDT(res_main)
  setDT(res_scapissimpler)
  setDT(res_sa1)
  setDT(res_sa2)
  
  fix.exp.fun <- function(var){
    factor(var, levels = c("Class_Peni_BetaS","Class_Peni_BetaR","Class_Peni_Ext",
                           "Class_Peni_Comb","Class_cephalosporins","Class_macrolides",
                           "Class_lincosamides","Class_TCLs","Class_FQs","Class_SMZTMP",
                           "Class_NIT") ,
           labels = c("Penicillin V","Flucloxacillin", "Penicillin ES",
                      "Amox-clav","Cephalosporins","Macrolides",
                      "Clindamycin","Tetracyclines","Fluoroquinolones","SMZ-TMP",
                      "Nitrofurantoin"))
  }
  
  res_main[, antibiotic := fix.exp.fun(antibiotic)]
  res_scapissimpler[, antibiotic := fix.exp.fun(antibiotic)]
  res_sa1[, antibiotic := fix.exp.fun(antibiotic)]
  res_sa2[, antibiotic := fix.exp.fun(antibiotic)]
  
  # Set order 
  res_main <- res_main[order(exposure, outcome), ]
  
  # Create a key for joining
  res_main[, key_id := .I]
  
  # Merge the key_id into the other tables using exposure and outcome
  res_scapissimpler <- merge(res_scapissimpler, res_main[, .(exposure, outcome, key_id)], by = c("exposure", "outcome"))
  res_sa1 <- merge(res_sa1, res_main[, .(exposure, outcome, key_id)], by = c("exposure", "outcome"))
  res_sa2 <- merge(res_sa2, res_main[, .(exposure, outcome, key_id)], by = c("exposure", "outcome"))
  
  # Order them by key_id
  setorder(res_scapissimpler, key_id)
  setorder(res_sa1, key_id)
  setorder(res_sa2, key_id)
  
  
  # Cor plot function 
  corplot_fun <- function(abx , res){
    
    res_abx <- res[antibiotic == abx,][["beta"]]
    res_main_abx <- res_scapissimpler[antibiotic == abx, ][["beta"]]
    
    spcor <- round(cor(res_abx, res_main_abx, method = "spearman"),2)
    spcor <- paste("Sp.cor=", spcor)
    
    min_lim <- min(c(res_abx, res_main_abx))
    max_lim <- max(c(res_abx, res_main_abx))
    
    min_lim <- min_lim + 0.2*min_lim
    max_lim <- max_lim + 0.2*max_lim
    
    min_x <- min_lim - 0.15*min_lim
    max_y <- max_lim - 0.17*max_lim

    
    data.frame(main = res_main_abx, sa = res_abx) %>%
      ggplot(aes(x=main, y= sa)) + 
      annotate("text", x = -Inf, y = Inf, label = spcor, size = 2.3, hjust = -0.1, vjust=1.5) +
      geom_abline(slope=1,intercept = 0, linetype="dashed", col = "gray70") +
      geom_point(col = "gray50", size = 1, alpha = .4) + 
      stat_smooth(method = "lm", col = "black", linewidth=1) + 
      scale_x_continuous(expand = c(0,0), limits = c(min_lim, max_lim)) + 
      scale_y_continuous(expand = c(0,0), limits = c(min_lim, max_lim)) +
      theme_classic() + 
      coord_cartesian() +
      ggtitle(abx) +
      xlab("full model") + ylab("sens. analysis") +
      theme(plot.title = element_text(size=7, hjust=.5,margin = margin(b = 0.1,0,0,0, unit= "mm")),
            axis.title.x = element_text(size=7, margin = margin(t = -0.3, b = 0.3 ,0,0, unit = "mm")), 
            axis.title.y = element_text(size=7, margin = margin(0,0,0,0)), 
            axis.text = element_text(size=6, margin = margin(b = -0.3,0,0,0, unit = "mm")), 
            plot.margin = margin(t=1,r = 1,0,0, unit = "mm"), 
            panel.border = element_rect( colour = "black", fill = NA, linewidth = 0.5), 
            axis.line = element_blank(), 
            axis.ticks.length = unit(.2 , units = "mm"), 
            aspect.ratio = 1)
    
  }
  
  antibiotics <- unique(res_main$antibiotic)
  
  antibiotics <- c("Clindamycin", "Flucloxacillin", "Fluoroquinolones", "Tetracyclines", "Cephalosporins", "Macrolides", "Penicillin V", "Penicillin ES", "Amox-clav",
                   "SMZ-TMP", "Nitrofurantoin")
  
  list_plots <- lapply(antibiotics, corplot_fun, res = res_sa1)
  p1 <- wrap_plots(list_plots, ncol = 4, nrow=3)
  
  list_plots <- lapply(antibiotics, corplot_fun, res = res_sa2)
  p2 <- wrap_plots(list_plots, ncol = 4, nrow=3)
  
  pfinal <- cowplot::plot_grid(p1, p2, ncol = 1, labels = c("a","b"), label_size = 7,
                               label_x = 0.12,
                               label_y = 0.98)


  ggsave(filename = "../../Revision_3/SupplFigures/SupplFig_10.png",plot=pfinal, width = 180, height = 200, units = "mm", dpi = 400)


  
  
