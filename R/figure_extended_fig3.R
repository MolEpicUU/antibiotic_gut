# Project Antibiotic use and the gut microbiota

# Extended Data Figure 

# Bar plots with the number of individuals exposed to a single antibiotic course in the past 8 years

library(data.table)
library(ggplot2)
library(patchwork)
library(readxl)
library(dplyr)
library(tidyr)

setwd('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/')
load("../revision__work/antibiotic_order_figures") # order_atb object

table_nexposed <- readxl::read_excel('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__tables/Supp.tables_temp.xlsx', sheet = "Suppl. Table 8",)

setDT(table_nexposed)
table_nexposed <- table_nexposed[outcome == "shannon" & model == "full model", .(antibiotic, period, Nexposed_SCAPIS, Nexposed_SIMPLER, Nexposed_MOS)]
table_nexposed[, Total := sum(Nexposed_SCAPIS, Nexposed_SIMPLER, Nexposed_MOS), by = c("antibiotic", "period")]
table_nexposed[, antibiotic := factor(antibiotic,
                               c("Clindamycin","Flucloxacillin","Fluoroquinolones","Tetracyclines",
                                 "Cephalosporins","Macrolides","Penicillin V","Penicillins extended spectrum",
                                 "Amoxicillin-clavulanic acid","Sulfamethoxazole-trimethoprim", "Nitrofurantoin"),
                               c("Clindamycin","Flucloxacillin","Fluoroquinolones","Tetracyclines",
                                 "Cephalosporins", "Macrolides","Penicillin V", "Penicillin ES",
                                 "Amox-clav", "SMZ-TMP","Nitrofurantoin"))]



colnames(table_nexposed) <- gsub("Nexposed_", "", colnames(table_nexposed), )

table_nexposed_wide <- table_nexposed %>% pivot_longer(cols = c("SCAPIS", "SIMPLER", "MOS", "Total"), names_to = "cohort", values_to = "N")
table_nexposed_wide$cohort <-  factor(table_nexposed_wide$cohort, c("Total","SCAPIS", "SIMPLER", "MOS"))
table_nexposed_wide$period <-  factor(table_nexposed_wide$period, c("4-8 years", "<4 years"))


N_anyatb <- sum(table_nexposed_wide$N[table_nexposed_wide$cohort=="Total"])

p_barplot <- table_nexposed_wide %>% filter(cohort == "Total") %>%  
  mutate(antibiotic = factor(antibiotic, order_atb)) %>% 
  ggplot(aes(x = antibiotic, y = N)) +
  geom_bar(position = position_dodge(width = 1), stat = "identity", fill = "gray60") +
  facet_grid(~period) +
  scale_y_continuous(breaks = c(0, 50 , 100, 150, 200, 400, 600, 800), expand = c(0, 0 ), limits = c(0,805)) +
  ylab("Number of individuals") +
  theme_classic() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7, angle = 25, hjust = 1) ,
        axis.text.y = element_text(size=6),
        axis.title = element_text(size = 7),
        strip.background=element_blank(),
        strip.text.x = element_text(size=7),
        axis.line.x = element_blank(),
        panel.border = element_rect(colour = "gray50", fill = NA),
  )


# Panel B ####


# Import data set 

resmeta_alpha <- fread('meta_alpha_singledose.tsv')
scapis <- fread('scapis_alpha_singledose.tsv')
mos <- fread('mos_alpha_singledose.tsv')
simpler <- fread('simpler_alpha_singledose.tsv')
Nexp <- rbind(scapis[, .(cohort,model, exposure,N.y)], simpler[, .(cohort,model, exposure,N.y)], mos[, .(cohort, model, exposure,N.y)])
Nexp <- unique(Nexp[model == "full.model", ])
Nexp <- Nexp[, N:= sum(N.y), by = exposure]
resmeta_alpha <- merge(resmeta_alpha[model=="full.model",], unique(Nexp[, .(exposure, N)]), by = "exposure" , all.x=T, all.y=F)

setnames(resmeta_alpha, c("LCI","HCI"), c("lci","hci"))
resmeta_alpha <- resmeta_alpha[,.(exposure, outcome, beta, lci, hci, p.value, cohort, model,N)]

resmeta_alpha[, q.value := p.adjust(p.value, method = "BH"), by = c("model","cohort")]

resmeta_alpha[, exposure := gsub("class","", exposure)]


x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                       "cephalosporins","macrolides",
                       "lincosamides","TCLs","FQs","SMZTMP","NIT"))

lev = do.call(c,lapply(x, function(w) paste(w, c("4_8yr","1_4yr"), sep="")))

lab = c("Penicillin V", "Flucloxacillin", "Penicillin ES",  "Amox-clav", 
        "Cephalosporins",
        "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")

lab = do.call(c,lapply(lab, function(w) paste(w, c("4-8yr","<4yr"))))


res <- copy(resmeta_alpha)
res[, exp_atb := gsub("[1-9].*", "", exposure)]
res[, period := stringr::str_extract(exposure, "[1-9].*")]
res[, exp_atb := gsub("class", "", exp_atb)]

exp_to_keep <- res[cohort == "Meta-analysis" & model == "full.model" & q.value<.05, exp_atb]


res[, exposure := factor(exposure, lev, lab)]
res <- res[order(model, cohort, outcome ,exposure), ]
res[, exposure := factor(exposure, rev(lab), rev(lab))]
res[, model := factor(model, rev(c("full.model","basic.model")), rev(c("Full model", "Basic model")))]
res[, fill := as.character(period)]
res[, fill := ifelse(q.value>.05, "No", fill )]
res[, outcome := factor(outcome, c("shannon","richness","invsimpson"), c("Shannon","Richness","Inv. Simpson"))]

res <- res[ model == "Full model", ]


lab.plot <- gsub("4-8"," 4-8", unique(res$exposure[res$exp_atb!="N"]))
lab.plot <- gsub(".*<4","<4", lab.plot)



# Order exposures 
order_plot <- do.call(c, lapply(order_atb, function(x) paste(x,c("4-8yr","<4yr"))))

min_y <- c(0,seq(4.5,26.5,by = 4))
max_y <- seq(2.5,26.5,by = 4)
n_list <- length(min_y)

ret_obj <- lapply(1:n_list, function(i) {
  geom_rect(color = "white",fill="gray97", xmin = -Inf, xmax = Inf, ymin = min_y[i], ymax = max_y[i], alpha = 0.3 )
})


psingle <- res %>% mutate(exposure = factor(exposure, rev(order_plot))) %>% 
  filter(!is.na(exposure)) %>% 

  ggplot(aes(x=beta, y = exposure, xmin = lci, xmax = hci, fill = fill, col = period)) +
  ret_obj[1] + ret_obj[2] + ret_obj[3] + ret_obj[4] + ret_obj[5] + ret_obj[6] + 
  geom_vline(xintercept=0, lty=2, linewidth=.3, color ="gray50") +
  geom_linerange(linewidth=0.4, orientation = "y",position=position_dodge(width = .7)) +
  ylab("") + xlab("") +
  geom_point(stroke=.5, size=2, position=position_dodge(width = .7), shape = 23, size = 2, color ="gray50") + 
  theme(legend.title = element_blank()) +                 #  "#D2A0D2","#74DC97","#DCC2A7"
  scale_fill_manual(breaks = c("1_4yr","4_8yr","No"), values = c("#D2A0D2","#DCC2A7", "white") , guide = "none") +
  scale_color_manual(breaks = c("1_4yr","4_8yr"), values = c("#D2A0D2","#DCC2A7") , guide = "none") +
  scale_y_discrete(breaks = unique(res[exp_atb!="N", exposure]), labels = lab.plot) +
  theme_classic() + 
  xlab("regresison coefficient") + 
  facet_wrap(~outcome, scales = "free_x") + 
  theme(legend.title = element_blank(), 
        strip.background = element_blank(), 
        legend.position = "right", 
        strip.text.x = element_text(size=7), 
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7), 
        panel.spacing = unit(.3, "mm")) 


 final_plot <- cowplot::plot_grid(p_barplot, psingle, ncol = 1, labels = c("a", "b"), label_size = 7, 
                                  rel_heights = c(1, 1.2))


ggsave("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_3/Extended_Data_Figure3.tiff", final_plot,
       width = 180, height = 200 , units = "mm", dpi = 600)

message("Saved plot - END ")

