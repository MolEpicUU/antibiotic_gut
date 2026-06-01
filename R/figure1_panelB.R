# Forest plot with the associations between alpha diversity and previous antibiotic use 

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

rm(list=ls())

setwd('/nobackup/users/baldanzi/atb_gut/results')

  # Import order of antibiotics 
  order_atb <- fread("../work/antibiotics_number_of_users.tsv")
  order_atb <- order_atb[period == "total_8years"]
  order_atb <- order_atb[order(-proportion_users)]
  order_atb[antibiotic == "Penicillin extended spectrum", antibiotic := "Penicillin ES"]
  order_atb <- order_atb$antibiotic
  order_atb <- factor(order_atb, order_atb)
  save(order_atb, file="../work/antibiotic_order_figures")

  load("../work/antibiotic_order_figures")
  
  # Import data set 
  resmeta <- fread('meta_alpha.tsv')
  resmeta <- resmeta[cohort=="Meta-analysis" & model == "full.model",]
  resmeta <- resmeta[,.(exposure, outcome, beta, LCI, HCI, q.value, cohort, model)]
  
  x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                         "cephalosporins","macrolides",
              "lincosamides","TCLs","FQs","SMZTMP","NIT"))
  
  lev = do.call(c,lapply(x, function(w) paste(w, c("4_8yr","1_4yr","1yr"), sep="_")))
  
  lab = c("Penicillin V", "Flucloxacillin", "Penicillin ES",  "Amox-clav", 
          "Cephalosporins",
          "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")
  
  lab = do.call(c,lapply(lab, function(w) paste(w, c("4-8yr","1-4yr","<1yr"))))
  
  
  # keep only the atb with q<0.05 meta-analysis results 
  res <- rbind(resmeta)
  res[, exp_atb := gsub("[1-9].*", "", exposure)]
  res[, period := stringr::str_extract(exposure, "[1-9].*")]
  exp_to_keep <- res[model == "full.model" & q.value<.05, exp_atb]
  
  
  res[, exposure := factor(exposure, lev, lab)]
  res <- res[order(model, cohort, outcome ,exposure), ]
  res[, exposure := factor(exposure, rev(lab), rev(lab))]
  res[, fill := period]
  res[, fill := ifelse(q.value>.05, "No", fill )]
  res[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]
  
  res <- res[exp_atb %in% c(exp_to_keep), ]

  
  lab.plot <- gsub(".*1-4","1-4", unique(res$exposure))
  lab.plot <- gsub("4-8","  4-8", lab.plot)
  lab.plot <- gsub(".*<","<",lab.plot)

  
  order_plot <- order_atb

  
  order_plot <- do.call(c, lapply(order_plot, function(x) paste(x,c("4-8yr", "1-4yr", "<1yr"))))
  
  
  min_y <- c(0,seq(6.5,21.5,by = 6))
  max_y <- seq(3.5,21.5,by = 6)
  n_list <- length(min_y)
  
  ret_obj <- lapply(1:n_list, function(i) {
    geom_rect(color = "white",fill="gray98", xmin = -Inf, xmax = Inf, ymin = min_y[i], ymax = max_y[i], alpha = 0.3 )
  })
  
  forestplot_alpha <- res[exp_atb!="N", ] %>% mutate(exposure = factor(exposure, rev(order_plot))) %>% 
  ggplot(aes(x=beta, y = exposure, xmin = LCI, xmax = HCI, fill = fill, col = period)) +
    ret_obj[1] + ret_obj[2] + ret_obj[3] + ret_obj[4] +
    geom_vline(xintercept=0, lty=2, linewidth=.3, color ="gray50") +
    geom_linerange(linewidth=0.4, orientation = "y",position=position_dodge(width = .7)) +
    ylab("") + xlab("regression coefficient") +
    geom_point(stroke=.3, size=2, position=position_dodge(width = .7), shape = 23, col = "gray10") + 
    theme(legend.title = element_blank()) +                    #       "#D2A0D2","#74DC97","#DCC2A7"
    scale_fill_manual(breaks = c("1yr","1_4yr","4_8yr","No"), values = c("#D2A0D2","#74DC97","#DCC2A7", "white") , guide = "none") +
    scale_color_manual(breaks = c("1yr","1_4yr","4_8yr"), values = c("#D2A0D2","#74DC97","#DCC2A7") , guide = "none") +
    scale_y_discrete(breaks = unique(res[exp_atb!="N", exposure]), labels = lab.plot) +
    scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    theme_classic() + 
    facet_wrap(~outcome, scales = "free_x") + 
    theme(legend.title = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "right", 
          strip.text.x = element_text(size=7), 
          axis.text.x = element_text(size = 6), 
          axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size=7),
          panel.spacing = unit(1.3, "mm")) 
  
  forestplot_alpha
  
  ggsave(plot= forestplot_alpha, "figure1_panelB_foreplot_alpha.pdf", 
         width = 180, height = 110, units = "mm",
         path="wharf/baldanzi/baldanzi-sens2019512/atbgut")
  
  
  