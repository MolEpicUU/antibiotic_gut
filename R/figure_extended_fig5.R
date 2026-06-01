# Project antibiotic and gut microbiome 


# The aim of this script is to create a heatmap showing the associations of 
# species previously connected to CRC and IBD  with antibiotic use
  
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(stringr)
  
  setwd("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2")
  
  taxa <- fread('/Users/gabba126/Documents/PhD_projects/Microbiome/Taxonomy/Taxonomy_CHAMP.tsv')
  atbres <- fread('revision__results/meta_species_class_clean.tsv')
  load("revision__work/antibiotic_order_figures")
  
  atbres <- merge(taxa[, .(MGS, MainTax)], atbres, by.x = "MGS", by.y = "outcome")
  atbres <- atbres[model == "full.model"]
  
  atb_linked_species <- unique(atbres$MainTax[atbres$q.value < 0.05])
  
  
  # CRC species 
  crc <- xlsx::read.xlsx("revision__externaldata_clinicaloutcomes/CRC_Piccinno_41591_2025_3693_MOESM2_ESM.xlsx", 
                         sheetName = "Table 6", startRow = 2)
  setDT(crc)
  crc <- crc[Feature.type == "SGBs" & Comparison == "control vs CRC"]
  crc <- crc[, c("Microbial.feature", "Effect.size.1", "P.1", "FDR.1")]
  crc[, Microbial.feature := gsub("\\|t__", "", str_extract(Microbial.feature, "^.*\\|t__"))]
  crc[, Microbial.feature := gsub("_", " ", Microbial.feature)]
  crc[, Microbial.feature := gsub("SGB\\d+", "", Microbial.feature)]
  crc <- crc[!grepl("^GGB| bacterium| sp |Incertae", Microbial.feature)]
  crc[, Microbial.feature := trimws(Microbial.feature)]
  crc <- crc[lengths(strsplit(Microbial.feature, "\\s+"))>=2]
  
  crc[, species_str := gsub(" ", ".*", Microbial.feature)]
  setnames(crc, c("Microbial.feature",  "Effect.size.1", "P.1", "FDR.1"), c("species", "effect_size", "pvalue", "qvalue"))
  crc <- crc[qvalue<0.05]
  crc[effect_size>0, .N]
  crc_species_pos_total <- crc %>% filter(effect_size>0) %>% pull(species_str) %>% unique(.)
  length(crc_species_pos_total)
  crc_species_pos_str <- paste0(crc_species_pos_total, collapse = "|") 
  crc_species_pos <- unique(grep(crc_species_pos_str, atb_linked_species, value=T))
  length(crc_species_pos)

  crc[effect_size<0, .N]
  crc_species_neg_total <- crc %>% filter(effect_size<0) %>% pull(species_str) %>% unique(.)
  length(crc_species_neg_total)
  crc_species_neg_str <- paste0(crc_species_neg_total, collapse = "|") 
  crc_species_neg <- unique(grep(crc_species_neg_str, atb_linked_species, value=T))
  length(crc_species_neg)
  
  # IBD #### ===================================================================
  
  ibd <- xlsx::read.xlsx("revision__externaldata_clinicaloutcomes/Nagata_SupplData_41467_2024_54797_MOESM6_ESM.xlsx", sheetIndex = 1, startRow = 2)

   setDT(ibd)
   setnames(ibd, c("Feature", "Coefficient..IBD.Japanese.4D.", "P.value..IBD.Japanese.4D.", "Q.value..IBD.Japanese.4D."),
            c("species", "coef_IBD", "pval_IBD", "qval_IBD"))

  ibd <- ibd[qval_IBD<0.05]
  ibd <- ibd[, c("species", "coef_IBD", "pval_IBD", "qval_IBD")]
  ibd[, species := gsub(" \\[ID\\:\\d+\\]", "", species)]
  ibd[, species := gsub("\\[|\\]", "", species)]
  ibd <- ibd[!grepl("sp\\.", species)]
  ibd <- ibd[!grepl("Unassigned species| bacterium|^bacterium|species incertae sedis|\\/", species)]
  ibd[, species := trimws(species)]
  ibd <- ibd[length(unlist(strsplit(species, "\\s+")))>=2]
  ibd[, species_str := gsub(" ", ".*", species)]


  ibd_species_pos_total <- ibd[coef_IBD > 0 ] %>%  pull(species_str) %>% unique(.)
  ibd_species_pos_str <- paste0(ibd_species_pos_total, collapse = "|")
  ibd_species_pos <- unique(grep(ibd_species_pos_str, atb_linked_species, value=T))

  ibd_species_neg_total <- ibd[coef_IBD < 0] %>% pull(species_str) %>% unique(.)
  ibd_species_neg_str <- paste0(ibd_species_neg_total, collapse = "|")
  ibd_species_neg <- unique(grep(ibd_species_neg_str, atb_linked_species, value=T))
  
  print(length(crc_species_pos_total))
  print(length(crc_species_neg_total))
  cat("Total species-level associations with CRC")
  print(length(c(crc_species_pos_total,crc_species_neg_total)))
  print(length(ibd_species_pos_total))
  print(length(ibd_species_neg_total))
  cat("Total species-level associations with IBC")
  print(length(c(ibd_species_pos_total,ibd_species_neg_total)))
  
  print(length(crc_species_pos))
  print(length(crc_species_neg))
  cat("CRC mapped to CHAMP")
  print(length(c(crc_species_pos, crc_species_neg)))
  print(length(ibd_species_pos))
  print(length(ibd_species_neg))
  cat("IBD mapped to CHAMP")
  print(length(c(ibd_species_pos, ibd_species_neg)))
  
  print(length(unique(c(crc_species_pos, crc_species_neg, ibd_species_pos, ibd_species_neg))))
  
  
  # Heatmap - prepare data ####
  
  x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                         "cephalosporins","macrolides",
                         "lincosamides","TCLs","FQs","SMZTMP","NIT"))
  
  lev = do.call(c,lapply(x, function(w) paste(w, c("4_8yr","1_4yr","1yr"), sep="_")))
  
  lab = c("Penicillin V", "Flucloxacillin", "Penicillin ES",  "Amox-clav", 
          "Cephalosporins",
          "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")
  
  lab = do.call(c,lapply(lab, function(w) paste(w, c("4-8yr","1-4yr","<1yr"))))
  
  
  order_plot <- order_atb
  
  
  order_plot <- do.call(c, lapply(order_plot, function(x) paste(x,c("<1yr", "1-4yr", "4-8yr"))))
  
  atbres[, exposure := factor(exposure, lev, lab)]
  
  atbres <- atbres[MainTax %in% c(crc_species_pos, crc_species_neg, ibd_species_pos, ibd_species_neg)] 
  
  mat_function <- function(dat, row, col){
    
    temp_dat <- dat %>% pivot_wider(id_cols = MainTax, values_from = "beta", names_from = "exposure") %>% as.data.frame(.)
    rownames(temp_dat) <- temp_dat$MainTax
    mat_beta <- as.matrix(temp_dat[, -1])
    
    temp_dat <- dat %>% pivot_wider(id_cols = MainTax, values_from = "q.value", names_from = "exposure") %>% as.data.frame(.)
    rownames(temp_dat) <- temp_dat$MainTax
    mat_q <- as.matrix(temp_dat[, -1])
    
    return(list(mat_beta = mat_beta, mat_q = mat_q))
    
    
  }
  
  abx_matrix <- mat_function(atbres)
  abx_matrix[[1]][is.na(abx_matrix[[1]])] <- 0
  abx_matrix[[2]][is.na(abx_matrix[[2]])] <- 1
  
  abx_matrix[[1]] <- abx_matrix[[1]][, order_plot] 
  abx_matrix[[2]] <- abx_matrix[[2]][, order_plot] 
  
  
  # associations with health outcomes 
  
  hm_species <- rownames(abx_matrix[[1]])
  
  mat_health <- data.frame(CRC = ifelse(hm_species %in% crc_species_pos, "increased", 
                           ifelse(hm_species %in% crc_species_neg, "decreased", NA)), 
                           IBD = ifelse(hm_species %in% ibd_species_pos, "increased", 
                                        ifelse(hm_species %in% ibd_species_neg, "decreased", NA)))
  
  rownames(mat_health) <- hm_species
  
  mat_health$CRC <- factor(mat_health$CRC, c("increased", "decreased", "0"))
  mat_health$IBD <- factor(mat_health$IBD, c("increased", "decreased", "0"))
  mat_health <- mat_health[order(mat_health$CRC, mat_health$IBD, rownames(mat_health)), ]
  
  mat_health <- as.matrix(mat_health)
  
  
  abx_matrix[[1]] <- abx_matrix[[1]][rownames(mat_health), ] 
  abx_matrix[[2]] <- abx_matrix[[2]][rownames(mat_health), ] 
  
  # Heatmap ####
  
  heatmap_colors <- colorRampPalette(c("dodgerblue3", "white", "mediumvioletred"))(5)
  
  HM1 <- Heatmap(mat_health, 
                 row_names_side = "left", 
                 row_names_gp = gpar(fontsize = 7), 
                 col = c("increased" = "mediumvioletred", "decreased" = "dodgerblue3"),
                 na_col = "gray95", 
                 column_names_side = "top", 
                 name="Health outcome", 
                 # column_names_centered = T,
                 width = unit(5, "mm"), 
                 column_names_gp = gpar(fontsize=7), 
                 column_title_rot = 0, 
                 column_title_gp = gpar(fontsize=7),
                 heatmap_legend_param = list(
                   direction = "horizontal",
                   title_position = "lefttop",
                   title_gp = gpar(fontsize = 7),  # title font size
                   labels_gp = gpar(fontsize = 5),
                   grid_height = unit(2, "mm")) )
  
  
  HM <-  Heatmap(abx_matrix[[1]], 
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(abx_matrix[[2]][i, j] < 0.05 & abx_matrix[[2]][i, j] >= 0.01) {
                     grid.text('*', x, y = y-unit(0.9, "mm"))
                   } 
                   if(abx_matrix[[2]][i, j] < 0.01) {
                     grid.text('**', x, y = y-unit(0.9, "mm"))
                   }},
                 show_row_dend = F, show_column_dend = F, cluster_rows = T,cluster_columns = F, 
                 row_names_side = "left", 
                 row_names_gp = gpar(fontsize = 7),  
                 column_labels = rep(c("<1yr", "1-4yr", "4-8yr"), 11), 
                 column_names_rot = 90,
                 column_names_centered = T,
                 column_names_gp = gpar(fontsize=7), 
                 column_title_rot = 10, 
                 column_title_gp = gpar(fontsize=7),
                 column_split = rep(order_atb, each = 3), 
                 column_names_side = "top", 
                 col = circlize::colorRamp2(c(-.8, 0, .8), heatmap_colors[c(1, 3, 5)]), 
                 name="Antibiotics", 
                 width = unit(ncol(abx_matrix[[1]]) * 4.2, "mm"), 
                 heatmap_legend_param = list(
                   direction = "horizontal",
                   title_position = "lefttop",
                   title_gp = gpar(fontsize = 7),  # title font size
                   labels_gp = gpar(fontsize = 6),
                   grid_height = unit(2, "mm")) ) 

  # draw(HM1+HM, heatmap_legend_side = "bottom")
  
  tryCatch(
    {dev.off()
      tiff(file="../Revision_3/Extended_Data_Figure5.tiff", width = 200, height=200 , units = "mm", res = 600)
      draw(HM1+HM, heatmap_legend_side = "bottom", padding = unit(c(0.1, .1, 0.1, 0.1), "mm"), ht_gap = unit(-3, "mm"))
      dev.off()}, 
    error = function(e) {
      tiff(file="../Revision_3/Extended_Data_Figure5.tiff", width = 200, height=200, units = "mm", res = 600 )
      draw(HM1+HM, heatmap_legend_side = "bottom", padding = unit(c(0.1, .1, 0.1, 0.1), "mm"), ht_gap = unit(-3, "mm"))
      dev.off()
    })
  
  