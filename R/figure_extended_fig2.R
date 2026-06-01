# Project Antibiotic use and the gut microbiota

# Forest plot of the negative control exposure 

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)


rm(list=ls())

setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/')
load("../revision__work/antibiotic_order_figures")

  # Import data set 

  var <- c('exposure', 'outcome', 'beta', 'LCI', 'HCI', 'q.value', 'cohort', 'model')
  resmeta             <- fread('meta_alpha_negexposure.tsv')
  resmeta_preatbzero  <- fread('meta_alpha_negexposure_preatbzero.tsv')
  
  resmeta <- resmeta[,.(exposure, outcome, beta, LCI, HCI, q.value, cohort, model)]
  resmeta_preatbzero <- resmeta_preatbzero[,.(exposure, outcome, beta, LCI, HCI, q.value, cohort, model)]
  
  
  lev <- paste0("Class_",c("lincosamides", "Peni_BetaR", "FQs","cephalosporins", 
                         "macrolides", "TCLs",
                         "Peni_BetaS","Peni_Comb","SMZTMP","NIT", "Peni_Ext"),"_after1yr")
  
  lab <-  c("Clindamycin","Flucloxacillin","Fluoroquinolones", "Cephalosporins",
          "Macrolides", "Tetracyclines", "Penicillin V", "Amox-clav", 
              "SMZ-TMP","Nitrofurantoin","Penicillin ES")
  
  lab <-  paste(lab, "1yr after")
  
  # fix names 
  resmeta[, exp_atb := gsub("[1-9].*", "", exposure)]
  resmeta[, exp_atb := gsub("after", "", exp_atb)]
  resmeta[, exposure := factor(exposure, rev(lev), rev(lab))]
  resmeta[, model := factor(model, c("full.model"), c("Full model"))]
  resmeta[, fill := as.character(cohort)]
  resmeta[, fill := ifelse(q.value>.05, "No", fill )]
  resmeta[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]
  
  
  resmeta_preatbzero[, exp_atb := gsub("[1-9].*", "", exposure)]
  resmeta_preatbzero[, exp_atb := gsub("after", "", exp_atb)]
  resmeta_preatbzero[, exposure := factor(exposure, rev(lev), rev(lab))]
  resmeta_preatbzero[, model := factor(model, c("full.model"), c("Full model"))]
  resmeta_preatbzero[, fill := as.character(cohort)]
  resmeta_preatbzero[, fill := ifelse(q.value>.05, "No", fill )]
  resmeta_preatbzero[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]

  
  # Order exposures 
  order_plot <- do.call(c, lapply(order_atb, function(x) paste(x,c("1yr after"))))
  
  min_y <- c(-Inf,seq(2.5,16.5,by = 2))
  max_y <- seq(1.5,16.5,by = 2)
  n_list <- length(min_y)
  
  ret_obj <- lapply(1:n_list, function(i) {
    geom_rect(color = "white",fill="gray97", xmin = -Inf, xmax = Inf, ymin = min_y[i], ymax = max_y[i], alpha = 0.3 )
  })
  
  
  pneg1 <- resmeta %>% filter(model=="Full model" & cohort =="Meta-analysis") %>% mutate(exposure = factor(exposure, rev(order_plot))) %>% 

    ggplot(aes(x=beta, y = exposure, xmin = LCI, xmax = HCI, color = cohort, fill = fill)) +
    ret_obj[1] + ret_obj[2] + ret_obj[3] + ret_obj[4] + ret_obj[5] + ret_obj[6] +
    geom_vline(xintercept=0, lty=2, linewidth=.3, color ="gray50") +
    geom_linerange(linewidth=0.4, orientation = "y",position=position_dodge(width = .7)) +
    ylab("") + xlab("regression coefficient") +
    geom_point(position=position_dodge(width = .7), shape = 23, size = 2) + 
    theme(legend.title = element_blank()) + # "#D2A0D2","#74DC97","#DCC2A7"
    scale_color_manual(breaks = c("Meta-analysis"), values = c("gray20")) + 
    scale_fill_manual(breaks = c("Meta-analysis", "No","Shannon","Richness","Inv. Simpson"), 
                      values = c("gray10", "white", "white","gray98","white"), 
                      guide = "none") + 
    theme_classic() + 
    facet_wrap("outcome", scales = "free_x") + 
    scale_x_continuous(expand = c(0.2,0), labels = label_number(drop0trailing = TRUE)) +
    theme(legend.title = element_blank(), strip.background = element_blank(), 
          strip.text = element_text(size=7), 
          legend.position = "right", 
          strip.text.x = element_text(size=7), 
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 7), 
          panel.spacing = unit(.5, "mm"),
          plot.title = element_text(hjust = .5)) +
    guides(color=guide_legend(override.aes=list(fill=NA)))
  
  
  pneg2 <- resmeta_preatbzero %>% filter(model=="Full model" & cohort =="Meta-analysis") %>% mutate(exposure = factor(exposure, rev(order_plot))) %>% 
    ggplot(aes(x=beta, y = exposure, xmin = LCI, xmax = HCI, color = cohort, fill = fill)) +
    ret_obj[1] + ret_obj[2] + ret_obj[3] + ret_obj[4] + ret_obj[5] + ret_obj[6] +
    geom_vline(xintercept=0, lty=2, linewidth=.3, color ="gray50") +
    geom_linerange(linewidth=0.4, orientation = "y",position=position_dodge(width = .7)) +
    ylab("") + xlab("regression coefficient") +
    geom_point(position=position_dodge(width = .7), shape = 23, size = 2) + 
    theme(legend.title = element_blank()) + # "#D2A0D2","#74DC97","#DCC2A7"
    scale_color_manual(breaks = c("Meta-analysis"), values = c("gray20")) + #, guide = "none") +
    scale_fill_manual(breaks = c("Meta-analysis", "No","Shannon","Richness","Inv. Simpson"), 
                      values = c("gray10", "white", "white","gray98","white"), 
                      guide = "none") + 
    theme_classic() + 
    facet_wrap(~outcome, scales = "free_x") + 
    scale_x_continuous(expand = c(0.1,0.1), labels = label_number(drop0trailing = TRUE)) +
    theme(legend.title = element_blank(), strip.background = element_blank(), 
          strip.text = element_text(size=7), 
          legend.position = "right", 
          strip.text.x = element_text(size=7), 
          axis.text.x = element_text(size = 6), 
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 7),
          plot.margin = margin(t=.1,r=0,b=.1,l=1, unit = "mm"), 
          panel.spacing = unit(.5, "mm"),
          plot.title = element_text(hjust = .5)) +
    guides(color=guide_legend(override.aes=list(fill=NA)))
  
  pfinal <- pneg1 + pneg2 + plot_layout(ncol = 2) + plot_annotation(tag_levels = "a") &
    theme(legend.position = "none", 
          plot.tag = element_text(size=7, face = "bold"), plot.tag.position = c(x=0.05,y=0.97), plot.margin = margin(t=0,b=0,r=0,l=0)) 
  
  pfinal
  
  
  ggsave("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_3/Extended_Data_Figure2.tiff", 
         pfinal, width = 180, height = 150, units = "mm", dpi = 600)
  
  message("Saved plot - END ")
  