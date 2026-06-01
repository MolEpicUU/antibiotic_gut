# Project Antibiotic use and the gut microbiota

# Forest plot with the results for alpha diversity 

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

rm(list=ls())

setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__results')

load("../revision__work/antibiotic_order_figures")

  # Import data set 
  resmeta <- fread('meta_alpha.tsv')
  resmeta <- resmeta[model == "full.model",]
  resmeta <- resmeta[,.(exposure, outcome, beta, LCI, HCI, q.value, cohort, model)]
  resmeta[, cohort := factor(cohort, rev(c("SCAPIS","SIMPLER", "MOS","Meta-analysis")))]
  
  x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                         "cephalosporins","macrolides",
              "lincosamides","TCLs","FQs","SMZTMP","NIT"))
  
  lev = do.call(c,lapply(x, function(w) paste(w, c("4_8yr","1_4yr","1yr"), sep="_")))
  
  lab = c("Penicillin V", "Flucloxacillin", "Penicillin ES",  "Amox-clav", 
          "Cephalosporins",
          "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")
  
  lab = do.call(c,lapply(lab, function(w) paste(w, c("4-8yr","1-4yr","<1yr"))))
  
  
  res <- rbind(resmeta)
  res[, exp_atb := gsub("[1-9].*", "", exposure)]
  res[, period := stringr::str_extract(exposure, "[1-9].*")]
  
  
  res[, exposure := factor(exposure, lev, lab)]
  res <- res[order(model, cohort, outcome ,exposure), ]
  res[, exposure := factor(exposure, rev(lab), rev(lab))]
  res[, fill := period]
  res[, fill := ifelse(q.value>.05, "No", fill )]
  res[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]
  

  
  lab.plot <- gsub(".*1-4","1-4", unique(res$exposure[res$exp_atb!="N"]))
  lab.plot <- gsub("4-8","  4-8", lab.plot)
  lab.plot <- gsub(".*<","<",lab.plot)

  
  
  order_plot <- order_atb

  
  order_plot <- do.call(c, lapply(order_plot, function(x) paste(x,c("4-8yr", "1-4yr", "<1yr"))))
  
  
  # Limit the CI that is shown 
  res[outcome == "Shannon", x_min := -0.4]
  res[outcome == "Shannon", x_max :=  0.7]
  res[outcome == "Richness", x_min := -100]
  res[outcome == "Richness", x_max :=  100]
  res[outcome == "Inv. Simpson", x_min := -10]
  res[outcome == "Inv. Simpson", x_max :=  25]
  
  # Create capped LCI and HCI columns for plotting
  res[, capped_LCI := LCI]
  res[, capped_HCI := HCI]
  res[outcome == "Shannon", capped_LCI := pmax(LCI, x_min)]
  res[outcome == "Shannon", capped_HCI := pmin(HCI,  x_max)]
  res[outcome == "Richness", capped_LCI := pmax(LCI, x_min)]
  res[outcome == "Richness", capped_HCI := pmin(HCI,  x_max)]
  res[outcome == "Inv. Simpson", capped_LCI := pmax(LCI, x_min)]
  res[outcome == "Inv. Simpson", capped_HCI := pmin(HCI,  x_max)]
  
  # Create indicators for whether CI was cropped
  res[, LCI_cropped := LCI < x_min]
  res[, HCI_cropped := HCI > x_max]
  res[LCI_cropped == F, x_min := NA]
  res[HCI_cropped == F, x_max := NA]
  
  
  min_y <- c(0,seq(6.5, 33.5,by = 6))
  max_y <- seq(3.5, 33.5,by = 6)
  n_list <- length(min_y)
  
  ret_obj <- lapply(1:n_list, function(i) {
    geom_rect(color = "white",fill="gray96", xmin = -Inf, xmax = Inf, ymin = min_y[i], ymax = max_y[i], alpha = 0.3 )
  })
  
  forestplot_alpha <- res[exp_atb!="N", ] %>% mutate(exposure = factor(exposure, rev(order_plot))) %>%
    
    ggplot(aes(y = exposure, fill = fill, col = period, group = cohort, shape = cohort)) +
    ret_obj[1] + ret_obj[2] + ret_obj[3] + ret_obj[4] + ret_obj[5] + ret_obj[6] +
    geom_vline(xintercept=0, lty=2, linewidth=.3, color ="gray50") +
    
    # Capped CI lines
    geom_linerange(aes(xmin = capped_LCI, xmax = capped_HCI),
                   linewidth=0.4, orientation = "y", position=position_dodge(width = .8)) +
    
    # Arrows for LCI cropped
    geom_segment(aes(x = x_min, xend = x_min + 0.01),
                 arrow = arrow(length = unit(0.1, "cm"), ends = "first"),
                 inherit.aes = T, position=position_dodge(width = .8)) +
    
    # Arrows for HCI cropped
    geom_segment(data = res,
                 aes(x = x_max, xend = x_max - 0.01),
                 arrow = arrow(length = unit(0.1, "cm"), ends = "first"),
                 inherit.aes = T, position=position_dodge(width = .8)) +
    
    # Points
    geom_point(aes(x=beta), stroke=.4, size=1.5, position=position_dodge(width = .8), col = "gray10") +
    ylab("") + xlab("regression coefficient") +
    theme(legend.title = element_blank()) +
    scale_fill_manual(breaks = c("1yr","1_4yr","4_8yr","No"),
                      values = c("#D2A0D2","#74DC97","#DCC2A7", "white"), guide = "none") +
    scale_color_manual(breaks = c("1yr","1_4yr","4_8yr"),
                       values = c("#D2A0D2","#74DC97","#DCC2A7"), guide = "none") +
    scale_shape_manual(breaks = c("SCAPIS","SIMPLER","MOS","Meta-analysis"),
                       values = c(21, 22, 24, 23)) +
    scale_y_discrete(breaks = unique(res[exp_atb!="N", exposure]), labels = lab.plot) +
    scale_x_continuous(labels = scales::label_number(drop0trailing = T)) + 
    theme_classic() +
    facet_wrap(~outcome, scales = "free_x") +
    theme(legend.title = element_blank(), 
          legend.margin = margin(0,0,0,0),
          strip.background = element_blank(), 
          legend.position = "top", 
          strip.text.x = element_text(size=10), 
          axis.text.x = element_text(size = 6), 
          axis.title.x = element_text(size = 8),
          panel.spacing = unit(.3, "mm"), 
          plot.margin = margin(t=0, r = 0, b = 0, l=0.1))

  forestplot_alpha
  
  ggsave(plot= forestplot_alpha, "figuresuppl_alpha_allcohorts.png", width = 8, height = 10, units = "in",
         path="/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/revision__figures", dpi = 400)
  
  