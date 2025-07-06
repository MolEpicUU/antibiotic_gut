# Marginal means of alpha diversity by antibiotic exposure (all antibiotics, not divided by class)


  library(data.table)
  library(Hmisc)
  library(rms)
  library(emmeans)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)

  scapis  <- fread('nobackup/users/baldanzi/atb_gut/work/scapis_working_dataset_revision.tsv', na.strings = c("", NA, "NA"))
  mos     <- fread('nobackup/users/baldanzi/atb_gut/work/mos_working_dataset_revision.tsv', na.strings = c("", NA, "NA"))
  simpler <- fread('nobackup/users/baldanzi/atb_gut/work/simpler_working_dataset_revision.tsv', na.strings = c("", NA, "NA"))
  
  load('nobackup/users/baldanzi/atb_gut/work/scapis_model_revision.Rdata')


  simpler[!is.na(center) & !is.na(dna_extraction_plate), site_plate := paste0(center, dna_extraction_plate)]
  
  simpler <- simpler[, c("age","sex","placebirth","smokestatus","education","site_plate","BMI","CCIw","ppi","statins","metformin","betablock","ssri","polypharmacy12m_cat", "N1yr","N1_4yr","N4_8yr","shannon","richness","invsimpson")]
  
  scapis <- scapis[,   c("age","Sex","placebirth","smokestatus","education","site_plate","BMI","CCIw","ppi","statins","metformin","betablock","ssri","polypharmacy12m_cat", "N1yr","N1_4yr","N4_8yr","shannon","richness","invsimpson")]
  
  mos <- mos[,         c("age","sex","country_birth","smoking", "education","extraction_plate","BMI", "CCIw","ppi","statins","metformin","betablock","ssri","polypharmacy_cat" ,"N1yr","N1_4yr","N4_8yr","shannon","richness","invsimpson")]
  
  colnames(mos) <- colnames(scapis) <- c("age","sex","placebirth","smokestatus","education","site_plate","BMI", "CCIw", "ppi","statins","metformin","betablock","ssri","polypharmacy12m_cat" ,"N1yr","N1_4yr","N4_8yr","shannon","richness","invsimpson")
  
  polldat <- rbind(scapis, simpler, mos)
  polldat <- polldat[complete.cases(polldat), ]
  
  # Data harmonization 
  polldat[, ppi := factor(ppi, c("0","no","1","yes"), c("no","no","yes","yes"))]
  polldat[, statins := factor(statins, c("0","no","1","yes"), c("no","no","yes","yes"))]
  polldat[, metformin := factor(metformin, c("0","no","1","yes"), c("no","no","yes","yes"))]
  polldat[, betablock := factor(betablock, c("0","no","1","yes"), c("no","no","yes","yes"))]
  polldat[, ssri := factor(ssri, c("0","no","1","yes"), c("no","no","yes","yes"))]
  polldat[, sex := factor(sex, c("1","male","M","2", "female","F"), c("male","male","male","female","female","female"))]
  polldat[placebirth == "1", placebirth := "Scandinavia"]
  polldat[placebirth != "Scandinavia", placebirth := "other"]
  polldat[smokestatus == "ex", smokestatus := "former"]
  
  
  # linear reg
  
  cov <- polldat[, c("age","sex","placebirth","smokestatus","education","site_plate","BMI","CCIw", "ppi","statins","metformin","betablock","ssri","polypharmacy12m_cat" )]
  mat <- as.matrix(polldat[, c("shannon","richness","invsimpson")])
  N1yr <- polldat$N1yr
  N1_4yr <- polldat$N1_4yr
  N4_8yr <- polldat$N4_8yr

  shannon_fit <- lm(mat[,"shannon"]  ~ . + rcs(N1yr, c(1,2,4)) + rcs(N1_4yr, 3) + rcs(N4_8yr, 3) , data = cov)
 
  richness_fit <- lm(mat[,"richness"]  ~ . + rcs(N1yr, c(1,2,4)) + rcs(N1_4yr, 3) + rcs(N4_8yr, 3) , data = cov)
 
  invsimpson_fit <- lm(mat[,"invsimpson"]  ~ . + rcs(N1yr, c(1,2,4)) + rcs(N1_4yr, 3) + rcs(N4_8yr, 3) , data = cov)
  
  calc_emmeans <- function(fit){
  rg <- ref_grid(fit, at = list(N1yr = c(0,1,2,3,4), N1_4yr = c(0) , N4_8yr = c(0)), nuisance = c("site_plate", "sex"),  
                 parms = c(1,2,4), weights = "proportional")
  datN1yr <- emmip(rg, ~ N1yr , CIs=TRUE,plotit=FALSE)
  datN1yr$period <- "<1yr"
  matN1yr <- emmeans(rg, ~ N1yr)
  matN1yr <- data.frame(pairs(matN1yr, adjust = "none"))
  matN1yr$q.value <- p.adjust(matN1yr$p.value, method = "BH")
  matN1yr$period <- "<1yr"
  
  rg <- ref_grid(fit, at = list(N1yr = c(0), N1_4yr = c(0,1,2,3,4) , N4_8yr = c(0)), nuisance = c("site_plate", "sex"),
                 parms = c(1,2,4), weights = "proportional")
  datN1_4yr <- emmip(rg, ~ N1_4yr , CIs=TRUE,plotit=FALSE)
  datN1_4yr$period <- "1-4yr"
  matN1_4yr <- emmeans(rg, ~ N1_4yr)
  matN1_4yr <- data.frame(pairs(matN1_4yr, adjust = "none"))
  matN1_4yr$q.value <- p.adjust(matN1_4yr$p.value, method = "BH")
  matN1_4yr$period <- "1-4yr"
  
  rg <- ref_grid(fit, at = list(N1yr = c(0), N1_4yr = c(0) , N4_8yr = c(0,1,2,3,4)), nuisance = c("site_plate", "sex"), 
                 parms = c(1,2,4), , weights = "proportional")
  datN4_8yr <- emmip(rg, ~ N4_8yr , CIs=TRUE,plotit=FALSE)
  datN4_8yr$period <- "4-8yr"
  matN4_8yr <- emmeans(rg, ~ N4_8yr)
  matN4_8yr <- data.frame(pairs(matN4_8yr, adjust = "none"))
  matN4_8yr$q.value <- p.adjust(matN4_8yr$p.value, method = "BH")
  matN4_8yr$period <- "4-8yr"
  
  dat <- rbind(datN1yr[,c("period","xvar","yvar","SE", "LCL", "UCL")], datN1_4yr[,c("period","xvar","yvar","SE", "LCL", "UCL")], datN4_8yr[,c("period","xvar","yvar","SE", "LCL", "UCL")])
  dat$period <- factor(dat$period, c("4-8yr", "1-4yr","<1yr" ))
  
  mat <- rbind(matN1yr[,c("period","contrast","estimate", "p.value","q.value")], matN1_4yr[,c("period","contrast","estimate","p.value","q.value")], matN4_8yr[,c("period","contrast","estimate","p.value","q.value")])
  mat$period <- factor(mat$period, c("4-8yr", "1-4yr","<1yr" ))
  
  return(list(dat=dat, mat = mat))
  
  }
 
  
  dat_mat_shannon <- calc_emmeans(shannon_fit)
  dat_mat_richness <- calc_emmeans(richness_fit)
  dat_mat_invsimpson <- calc_emmeans(invsimpson_fit)
  
  dat_shannon <- dat_mat_shannon[[1]]
  dat_richness <- dat_mat_richness[[1]]
  dat_invsimpson <- dat_mat_invsimpson[[1]]
  
  
  mat_shannon <- dat_mat_shannon[[2]]
  mat_richness <- dat_mat_richness[[2]]
  mat_invsimpson <- dat_mat_invsimpson[[2]]
  
  mat_shannon <- mat_shannon %>% separate(contrast, into = c("group1", "group2"), sep = " - ") %>% 
    mutate(group1 = gsub(".*yr","", group1), group2 = gsub(".*yr","", group2), p.value = -log(p.value)) 
  mat_richness <- mat_richness %>% separate(contrast, into = c("group1", "group2"), sep = " - ") %>% 
    mutate(group1 = gsub(".*yr","", group1), group2 = gsub(".*yr","", group2), p.value = -log(p.value)) 
  mat_invsimpson <- mat_invsimpson %>% separate(contrast, into = c("group1", "group2"), sep = " - ") %>% 
    mutate(group1 = gsub(".*yr","", group1), group2 = gsub(".*yr","", group2), p.value = -log(p.value)) 
    
  
     
  plot_fun <- function(dat, title) {
         ggplot(dat, aes(x=as.factor(xvar),y=yvar, group = period, fill = period))+
           
             geom_line(linewidth = .2, color = "gray30") +
             geom_errorbar(aes(ymin=LCL,ymax=UCL), width=0, color = "gray75", linewidth = 1) + 
             geom_point(shape=22, size = 2, color = "gray20") + 
             
             facet_wrap("period", ncol = 3) + theme_classic() +
             ylab(title) + xlab("Number of antibiotic courses") +
             scale_fill_manual(breaks = c("<1yr","1-4yr","4-8yr"), values = c("#D2A0D2","#74DC97","#DCC2A7") , guide = "none") +
             theme(plot.title = element_text(hjust = .5, size = 9),
                   axis.title.x = element_text(size=8, margin = margin(b=3)), 
                   plot.margin = margin(r=.1, l=.1, b = 1))
  }
  
  
  levels_d <- c("0","1","2","3","4")
  
  
  heatmap_fun <- function(dat, title) {
    dat %>%
      mutate(
        group1 = factor(group1, levels = levels_d),
        group2 = factor(group2, levels = rev(levels_d)),
        star = ifelse(q.value < 0.05, "*", "")  # create a column for stars
      ) %>%
      ggplot(aes(x = group1, y = group2, fill = estimate)) +
      geom_tile(color = "white") +
      ggtitle(title) + 
      geom_text(aes(label = star), color = "black", size = 4) +  # add stars
      scale_fill_gradient2(mid = "gray98", high = "red", midpoint = 0, name = "difference in EMM") +
      scale_y_discrete(expand=c(0,0)) + 
      theme_minimal() +
      ylab("N courses") + xlab("N courses") +
      theme(plot.title = element_text(hjust = .5, size = 9), 
            plot.margin = margin(r=.1, l=.1, t = 1),
            axis.title = element_text(size=8), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(margin = margin(r=3), size=6),
            axis.title.x = element_text(size=6),
            axis.text = element_text(size=6), 
            legend.position = "bottom", 
            strip.text.x = element_text(size=8),
            strip.placement = "outside", 
            strip.text.y.left = element_text(angle = 0),        
            strip.text = element_text(margin = margin(b = 4))) +
      facet_wrap(~period, nrow = 1)
  }
     
       
      plot_shannon      <- plot_fun(dat_shannon, "Shannon") 
      plot_richness     <- plot_fun(dat_richness, "Richness")
      plot_invsimpson   <- plot_fun(dat_invsimpson, "Inv. Simpson")
      
      hm_shannon        <- heatmap_fun(mat_shannon, "Shannon") 
      hm_richness       <- heatmap_fun(mat_richness, "Richness")
      hm_invsimpson     <- heatmap_fun(mat_invsimpson, "Inv. Simpson")
 
 pfinal <- plot_shannon + plot_richness + plot_invsimpson + 
   hm_shannon + hm_richness + hm_invsimpson +  plot_layout(ncol = 3, nrow = 2, 
                                                           heights = c(5,1)) &theme(
     legend.justification = c(1, 0),
     legend.background = element_rect(fill = "white", color = NA),
     legend.margin = margin(t = -5),
     legend.key.height = unit(0.2, "cm"),
     legend.text = element_text(size=6), 
     legend.title = element_text(angle = 0, hjust = 0.5, size = 8)
   )  &
   guides(
     fill = guide_colourbar(
       title.position = "left",
       title.hjust = 0,
       label.position = "bottom",
       direction = "horizontal"
     ))
 
  ggsave(pfinal, file = "figure1_panelA.pdf", width=8, height = 4,
        path = 'wharf/baldanzi/baldanzi-sens2019512/atbgut')
  
  