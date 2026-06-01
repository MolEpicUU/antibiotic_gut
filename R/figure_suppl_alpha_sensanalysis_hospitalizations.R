# Project Antibiotic use and the gut microbiota

# Forest plot of sensitivity analysis with alpha diversity removing hospitalized participants 

  
  rm(list=ls())
  
  library(data.table)
  library(ggplot2)
  library(ggh4x)
  library(dplyr)
  

  setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results')  
  load("../revision__work/antibiotic_order_figures") # order_atb object
    
  resmeta_sa <- fread('meta_alpha_sa.tsv')
  resmeta_sa[, atb := gsub("Class_","",gsub("_[1,4].*r","",exposure)),]
  resmeta_sa[, period := stringr::str_extract(exposure, "[1,4].*"),]
  resmeta_sa <- resmeta_sa[cohort == "Meta-analysis",]
  
  resmeta <- fread('meta_alpha.tsv')
  resmeta[, atb := gsub("Class_","",gsub("_[1,4].*r","",exposure)),]
  resmeta[, period := stringr::str_extract(exposure, "[1,4].*"),]
  resmeta <- resmeta[cohort=="Meta-analysis" & model == "full.model",]
  
  
  
  # Divide the two sensitivity analyses results 
  resmeta_sa[, type_sa := model]
  resmeta_sa[, type_sa := factor(type_sa, c("full.model_hospitalizedinfect","full.model_hospitalizedgeneral", "full.model"), 
                                 c("Exclusion hospitalized\nfor infection", "Exclusion hospitalized\nany reason", "Full model")) ]
  
  resmeta[, type_sa := "Full model"]
  
  resmeta_sa <- rbindlist(list(resmeta_sa, resmeta), fill = T)
  
  # Fix labels 
  x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                         "cephalosporins","macrolides",
                         "lincosamides","TCLs","FQs","SMZTMP","NIT"))
  
  lev = do.call(c,lapply(x, function(w) paste(w, c("4_8yr","1_4yr","1yr"), sep="_")))
  
  lab = c("Penicillin V", "Flucloxacillin", "Penicillin ES",  "Amox-clav", 
          "Cephalosporins",
          "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")
  
  lab = do.call(c,lapply(lab, function(w) paste(w, c("4-8yr","1-4yr","<1yr"))))
  
  
  resmeta_sa[, exp_atb := gsub("[1-9].*", "", exposure)]
  resmeta_sa[, period := stringr::str_extract(exposure, "[1-9].*")]
  resmeta_sa[, exposure := factor(exposure, lev, lab)]
  resmeta_sa <- resmeta_sa[order(model, cohort, outcome ,exposure), ]
  resmeta_sa[, exposure := factor(exposure, rev(lab), rev(lab))]
  resmeta_sa[, fill := type_sa]
  resmeta_sa[q.value>0.05, fill := "No"]
  resmeta_sa[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]
  
  lab.plot <- gsub(".*1-4","1-4", unique(resmeta_sa$exposure))
  lab.plot <- gsub("4-8","  4-8", lab.plot)
  lab.plot <- gsub(".*<","<",lab.plot)
  
  
  
  # Order exposures 
  order_plot <- do.call(c, lapply(order_atb, function(x) paste(x,c("4-8yr", "1-4yr", "<1yr"))))
  
  breaks_names <- c("Full model", "Exclusion hospitalized\nfor infection", "Exclusion hospitalized\nany reason")
  lab_names <- c("Full model n = 14974", 
                 "Exclusion hospitalized for infection\nn included = 12727", 
                 "Exclusion hospitalized any reason\nn = 8135")
  
  resmeta_sa[, type_sa := factor(type_sa, rev(breaks_names))]
  
  
  min_y <- c(0,seq(6.5,33.5,by = 6))
  max_y <- seq(3.5,33.5,by = 6)
  n_list <- length(min_y)
  
  ret_obj <- lapply(1:n_list, function(i) {
    geom_rect(color = "white",fill="gray97", xmin = -Inf, xmax = Inf, ymin = min_y[i], ymax = max_y[i], alpha = 0.3 )
  })
  
  forestplot_alpha <- resmeta_sa  %>% mutate(exposure = factor(exposure, rev(order_plot))) %>% 
  
    ggplot(aes(x=beta, y = exposure, xmin = LCI, xmax = HCI, fill = fill, col = type_sa, group = type_sa)) +
    ret_obj[1] + ret_obj[2] + ret_obj[3] + ret_obj[4] + ret_obj[5] + ret_obj[6] +
    geom_vline(xintercept=0, lty=2, linewidth=.2, color ="gray60") +
    geom_linerange(linewidth=0.5, orientation = "y",position=position_dodge(width = .8)) +
    ylab("") + xlab("regression coefficient") +
    geom_point(stroke=.4, size=1.5, position=position_dodge(width = .8), shape = 23) + #, color ="gray50") + 
    theme(legend.title = element_blank()) + #                     "#D2A0D2","#74DC97","#DCC2A7"
    scale_fill_manual(breaks = c(breaks_names, "No"), values = c("black","#E69F00","gray50", "white") , guide = "none") +
    scale_color_manual(breaks = c(breaks_names), values = c("black","#E69F00","gray50") , labels = lab_names) +
    scale_y_discrete(breaks = unique(resmeta_sa$exposure), labels = lab.plot) +
    theme_classic() + 
    facet_nested_wrap(~outcome, scales = "free_x", ncol = 6, nrow = 1,solo_line = T, nest_line = element_line(color="black")) + 
    theme(legend.title = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "left", 
          legend.margin = margin(r=0, l = 0, b = 0, t = 0),
          legend.text = element_text(size=7),
          plot.margin = unit(c(l=0, t = 0, b = 0, r = -4), "mm"),
          strip.text.x = element_text(size=7),
          axis.title.x = element_text(size=7),
          axis.text.y = element_text(size=7),
          axis.text.x = element_text(size = 6), 
          panel.spacing = unit(1, "mm")) 
  
  
  ggsave(plot= forestplot_alpha, "SupplFig_9.png", width = 180, height = 200, units = "mm",
          path="/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_3/SupplFigures", dpi = 400)
  
  
