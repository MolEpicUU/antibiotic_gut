# Project Antibiotic use and the gut microbiota

# Sensivity analysis with different time windows for recent antibiotic use

# Forest plot 

library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)

setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/')

load("../revision__work/antibiotic_order_figures")

# Import data set 
 res <- fread('meta_alpha_sa_timewindows.tsv')
 res <- res[cohort == "Meta-analysis",]
 setnames(res, "model", "time_windown")

 res <- res[,.(exposure, outcome, beta, LCI, HCI, q.value, cohort, time_windown)]
 
 x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                        "cephalosporins","macrolides",
                        "lincosamides","TCLs","FQs","SMZTMP","NIT"))
 
 lev = do.call(c,lapply(x, function(w) paste(w, c("4_8yr","1_4yr","1yr"), sep="_")))
 
 lab = c("Penicillin V", "Flucloxacillin", "Penicillin ES",  "Amox-clav", 
         "Cephalosporins",
         "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")
 
 lab = do.call(c,lapply(lab, function(w) paste(w, c("4-8yr","1-4yr","<1yr"))))
 
 
 res[, exp_atb := gsub("[1-9].*", "", exposure)]
 res[, period := stringr::str_extract(exposure, "[1-9].*")]
 
 res[, exposure := factor(exposure, lev, lab)]
 res <- res[order(cohort, outcome ,exposure), ]
 res[, exposure := factor(exposure, rev(lab), rev(lab))]
 res[, fill := period]
 res[, fill := ifelse(q.value>.05, "No", fill )]
 res[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]
 
 lab.plot <- gsub(".*1-4","1-4", unique(res$exposure))
 lab.plot <- gsub("*.4-8","4-8", lab.plot)
 lab.plot <- gsub(".*<","<",lab.plot)
 
 
 order_plot_atb <- order_atb
 
 
 order_plot <- do.call(c, lapply(order_plot_atb, function(x) paste(x,c("4-8yr", "1-4yr", "<1yr"))))
 
   
   p1 <- res %>% filter(outcome == "Shannon") %>% 
     mutate(atb = gsub(" $", "",str_extract(exposure, "^.* "))) %>% 
     mutate(atb = factor(atb, order_plot_atb)) %>% 
     mutate(exposure = factor(exposure, rev(order_plot))) %>%
     mutate(period = factor(period, c("4_8yr", "1_4yr", "1yr"), c("4-8yr", "1-4yr", "<1yr"))) %>%
   ggplot(aes(x=beta, y = time_windown, xmin = LCI, xmax = HCI, fill = fill, col = period, group = rev(period))) +   
     geom_vline(xintercept=0, lty=2, linewidth=.3, color ="darkred") +
     geom_linerange(linewidth=0.6, orientation = "y",position=position_dodge(width = .7)) +
     ylab("") + xlab("regression coefficient") + 
     ggtitle("Shannon") +
     geom_point(stroke=.4, size=1.8, position=position_dodge(width = .7), shape = 23, col = "gray20") +
     scale_fill_manual(breaks = c("1yr","1_4yr","4_8yr","No"), values = c("#D2A0D2","#74DC97","#DCC2A7", "white") , guide = "none") +
     scale_color_manual(breaks = rev(c("<1yr", "1-4yr","4-8yr" )), values = rev(c("#D2A0D2","#74DC97","#DCC2A7"))) +
     scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE)) + 
     theme_classic() + 
     theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
           plot.title = element_text(size=7, hjust=0, face = "bold"),
           strip.text = element_text(size = 7, margin = margin(.01, .1, .02, .1)),
           strip.background = element_blank(), 
           legend.title = element_blank(),
           axis.text.x = element_text(size = 6), 
           axis.text.y = element_text(size = 7), 
           axis.title.x = element_text(size=7),
           panel.spacing.y =  unit(0, "lines"),
           panel.spacing.x =  unit(0.3, "lines"), 
           plot.margin = margin(2, 2, 2, 2), 
           legend.position = c(1,0), 
           legend.margin = margin(unit = "mm", r = 5, b = 2), 
           legend.justification = c(.95,0), 
           legend.key.height = unit(3, "mm")) + 
     facet_wrap(~atb, scales = "free_x", ncol = 4, nrow= 3)  
   
   p2 <- res %>% filter(outcome == "Richness") %>% 
     mutate(atb = gsub(" $", "",str_extract(exposure, "^.* "))) %>% 
     mutate(atb = factor(atb, order_plot_atb)) %>% 
     mutate(exposure = factor(exposure, rev(order_plot))) %>%
     mutate(period = factor(period, c("4_8yr", "1_4yr", "1yr"), c("4-8yr", "1-4yr", "<1yr"))) %>%
     ggplot(aes(x=beta, y = time_windown, xmin = LCI, xmax = HCI, fill = fill, col = period, group = rev(period))) +   
     geom_vline(xintercept=0, lty=2, linewidth=.3, color ="darkred") +
     geom_linerange(linewidth=0.6, orientation = "y",position=position_dodge(width = .7)) +
     ylab("") + xlab("regression coefficient") +
     ggtitle("Richness") +
     geom_point(stroke=.4, size=1.8, position=position_dodge(width = .7), shape = 23, col = "gray20") +
     scale_fill_manual(breaks = c("1yr","1_4yr","4_8yr","No"), values = c("#D2A0D2","#74DC97","#DCC2A7", "white") , guide = "none") +
     scale_color_manual(breaks = rev(c("<1yr", "1-4yr","4-8yr" )), values = rev(c("#D2A0D2","#74DC97","#DCC2A7"))) +
     scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE)) + 
     theme_classic() + 
     theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
           plot.title = element_text(size=7, hjust=0, face = "bold"),
           strip.text = element_text(size = 7, margin = margin(.01, .1, .02, .1)),
           strip.background = element_blank(), 
           legend.title = element_blank(),
           axis.text.x = element_text(size = 6), 
           axis.text.y = element_text(size = 7), 
           axis.title.x = element_text(size=7),
           panel.spacing.y =  unit(0, "lines"),
           panel.spacing.x =  unit(0.3, "lines"), 
           plot.margin = margin(2, 2, 2, 2),
           legend.position = c(1,0),            
           legend.margin = margin(unit = "mm", r = 5, b = 2), 
           legend.justification = c(.95,0), 
           legend.key.height = unit(3, "mm")) +
     facet_wrap(~atb, scales = "free_x", ncol = 4, nrow= 3)  

   p3 <- res %>% filter(outcome == "Inv. Simpson") %>% 
     mutate(atb = gsub(" $", "",str_extract(exposure, "^.* "))) %>% 
     mutate(atb = factor(atb, order_plot_atb)) %>% 
     mutate(exposure = factor(exposure, rev(order_plot))) %>%
     mutate(period = factor(period, c("4_8yr", "1_4yr", "1yr"), c("4-8yr", "1-4yr", "<1yr"))) %>%
     ggplot(aes(x=beta, y = time_windown, xmin = LCI, xmax = HCI, fill = fill, col = period, group = rev(period))) +   
     geom_vline(xintercept=0, lty=2, linewidth=.3, color ="darkred") +
     geom_linerange(linewidth=0.6, orientation = "y",position=position_dodge(width = .7)) +
     ylab("") + xlab("regression coefficient") +
     ggtitle("Inv. Simpson") +
     geom_point(stroke=.4, size=1.8, position=position_dodge(width = .7), shape = 23, col = "gray20") +
     scale_fill_manual(breaks = c("1yr","1_4yr","4_8yr","No"), values = c("#D2A0D2","#74DC97","#DCC2A7", "white") , guide = "none") +
     scale_color_manual(breaks = rev(c("<1yr", "1-4yr","4-8yr" )), values = rev(c("#D2A0D2","#74DC97","#DCC2A7"))) +
     scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE)) + 
     theme_classic() + 
     theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
           plot.title = element_text(size=7, hjust=0, face = "bold"),
           strip.text = element_text(size = 7, margin = margin(.01, .1, .02, .1)),
           strip.background = element_blank(), 
           legend.title = element_blank(),
           axis.text.x = element_text(size = 6), 
           axis.text.y = element_text(size = 7), 
           axis.title.x = element_text(size=7),
           panel.spacing.y =  unit(0, "lines"),
           panel.spacing.x =  unit(0.3, "lines"), 
           plot.margin = margin(2, 2, 2, 2),            
           legend.position = c(1,0), 
           legend.margin = margin(unit = "mm", r = 5, b = 2), 
           legend.justification = c(.95,0), 
           legend.key.height = unit(3, "mm")) + 
     facet_wrap(~atb, scales = "free_x", ncol = 4, nrow= 3)  

   
   
   
 final_plot = (p1 + plot_spacer() + p2 + plot_spacer() + p3) +
   plot_layout(ncol = 1, heights = c(1, -0.05, 1, -0.05, 1)) +
   plot_annotation(tag_levels = "a")&
   theme(plot.margin = margin(t=.1, b=.1, l=-2, r=.1, unit = "mm"), plot.tag = element_text(size=7, face = "bold"), plot.tag.position = c(x=0.05,y=0.99),
         legend.text = element_text(size=7))
   

 #final_plot
 
ggsave(filename = "../../Revision_3/SupplFigures/SupplFig_8.png", 
        plot = final_plot, dpi = 400,
       width = 180, height = 210, units = "mm")      


     
     
     
     
     
     
     
     
     
     
     
     
     