# Project Antibiotic use and the gut microbiota


# alpha-diversity associations stratified by sex and age 

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

rm(list=ls())

setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results')

load("../revision__work/antibiotic_order_figures")

  # Import data set 
  resmeta <- fread('meta_alpha_sa__age_sex.tsv')
  resmeta <- resmeta[,.(exposure, outcome, beta, LCI, HCI, q.value, cohort, model)]

  
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
  
  # Clean model names ####
  res[, fill := model]
  res[, model := factor(model, c("female", "male", "age > 55 years", "age <= 55 years"), 
                        c("Females", "Males", "≤55 years", ">55 years"))]
  res[, fill := factor(fill, c("female", "male", "age > 55 years", "age <= 55 years", "No"), 
                        c("Females", "Males", "≤55 years", ">55 years", "No"))]
  
  
  res[q.value>=.05 , fill := "No"]
  res[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]
  
  lab.plot <- gsub(".*1-4","1-4", unique(res$exposure[res$exp_atb!="N"]))
  lab.plot <- gsub("4-8","  4-8", lab.plot)
  lab.plot <- gsub(".*<","<",lab.plot)

  
  
  order_plot <- order_atb

  
  order_plot <- do.call(c, lapply(order_plot, function(x) paste(x,c("4-8yr", "1-4yr", "<1yr"))))
  
  
  
  plot_fun <- function(dat, breaks, lab_values){
  
    min_y <- c(0,seq(6.5,33.5,by = 6))
    max_y <- seq(3.5,33.5,by = 6)
    n_list <- length(min_y)
    
    ret_obj <- lapply(1:n_list, function(i) {
      geom_rect(color = "white",fill="gray95", xmin = -Inf, xmax = Inf, ymin = min_y[i], ymax = max_y[i], alpha = 0.3 )
    })
  
  forestplot_alpha <- dat %>% mutate(exposure = factor(exposure, rev(order_plot))) %>% 
    filter(model %in% breaks) %>% 

  ggplot(aes(x=beta, y = exposure, xmin = LCI, xmax = HCI, fill = fill, col = model, group = model)) +
    ret_obj[1] + ret_obj[2] + ret_obj[3] + ret_obj[4] + ret_obj[5] + ret_obj[6] +  
    geom_vline(xintercept=0, lty=2, linewidth=.3, color ="gray50") +
    geom_linerange(linewidth=0.5, orientation = "y",position=position_dodge(width = .7)) +
    ylab("") + xlab("regression coefficient") + 
    geom_point(stroke=.4, size=1.5, position=position_dodge(width = .7), shape = 23) + # , col = "gray20") + 
    theme(legend.title = element_blank()) + #                     "#D2A0D2","#74DC97","#DCC2A7"
    scale_fill_manual(breaks = c(breaks,"No"), values = c("darkblue", "#E69F00", "white") , guide = "none") +
    scale_color_manual(breaks = rev(breaks), values = rev(c("darkblue", "#E69F00")), labels = rev(lab_values)) +
    scale_y_discrete(breaks = unique(res[exp_atb!="N", exposure]), labels = lab.plot) +
    scale_x_continuous(labels = scales::label_number(drop0trailing = T)) +
    theme_classic() + 
    facet_wrap(~outcome, scales = "free_x") + 
    theme(legend.title = element_blank(), 
          plot.title = element_text(size=7, margin = margin(b=0)),
          strip.background = element_blank(), 
          legend.position = "top", 
          strip.text.x = element_text(size=7, color = "black"), 
          axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 7, margin = margin(l=-3, unit = "mm")), 
          axis.title = element_text(size=7),
          legend.text = element_text(size = 7),
          plot.margin = margin(t=0,r=0,b=0,l=0), 
          panel.spacing = unit(.1, "mm")) 
  
  return(forestplot_alpha)
  }
  

  breaks_sex = c("Females", "Males")
  breaks_age = c("≤55 years", ">55 years")
  
  plot_sex <- plot_fun(res, breaks = breaks_sex, lab_values = c("Females n = 7591", "Males n = 7383"))

  
  # Age stratified 
  plot_age <- plot_fun(res, breaks = breaks_age, lab_values = c("≤55 years n = 4259", ">55 years n = 10715"))
  
  
  finalplot <- plot_sex + plot_age + plot_layout(ncol=2) + plot_annotation(tag_levels = "a") &theme(plot.tag = element_text(size=7, face = "bold"), plot.tag.position = c(x=0.1,y=0.97), plot.margin = margin(t=0,b=0,r=0,l=0), legend.margin = margin(b=-5, unit = "mm"))
  
  ggsave(plot= finalplot, "SupplFig_7.png", width = 180, height = 200, units = "mm", dpi = 400, 
         path="/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_3/SupplFigures")
  
  