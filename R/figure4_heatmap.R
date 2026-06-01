# Project Antibiotic Use and the Gut Microbiota

# Heatmap - health outcome


library(data.table)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

setwd("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results")


# import 
  res_heatlh <- fread("scapis_healthoutcomes_spearman.tsv")
  res_meta <- fread("meta_species_class_clean.tsv")
  
  taxonomy <- fread("/Users/gabba126/Documents/PhD_projects/Microbiome/Taxonomy/Taxonomy_CHAMP.tsv")

  
  res_heatlh <- merge(res_heatlh, taxonomy[, .(MGS, MainTax)], by.x="species", by.y="MGS", all.x=T)
  setcolorder(res_heatlh, c("MainTax"))
  
  
  # Clean antibiotic associations
  
  column_order <- c("Class_lincosamides_1yr","Class_lincosamides_1_4yr","Class_lincosamides_4_8yr",
                   "Class_Peni_BetaR_1yr","Class_Peni_BetaR_1_4yr","Class_Peni_BetaR_4_8yr",
                   "Class_FQs_1yr","Class_FQs_1_4yr","Class_FQs_4_8yr")
  
  res_meta <- res_meta[model == "full.model", ]
  res_meta <- res_meta[exposure %in% column_order, ]
  res_meta <- res_meta[outcome %in% res_heatlh$species, ]
  res_meta <- merge(res_meta, taxonomy[, .(MGS, MainTax)], by.x="outcome", by.y="MGS", all.x=T)
  setnames(res_meta, c("outcome"), c("species"))
  setnames(res_meta, c("exposure"), c("outcome"))
  
  
  mat_function <- function(dat, row, col){
    
    temp_dat <- dat %>% pivot_wider(id_cols = MainTax, values_from = "beta", names_from = "outcome") %>% as.data.frame(.)
    rownames(temp_dat) <- temp_dat$MainTax
    mat_beta <- as.matrix(temp_dat[, -1])
    
    temp_dat <- dat %>% pivot_wider(id_cols = MainTax, values_from = "q.value", names_from = "outcome") %>% as.data.frame(.)
    rownames(temp_dat) <- temp_dat$MainTax
    mat_q <- as.matrix(temp_dat[, -1])
    
    return(list(mat_beta = mat_beta, mat_q = mat_q))
    
    
  }
  
  setnames(res_heatlh, "rho", "beta")
  health_matrix <- mat_function(res_heatlh)
  abx_matrix <- mat_function(res_meta)
  
  # Max-scaling
  col_max <- apply(abs(health_matrix[[1]]), 2, max, na.rm = TRUE)
  
  health_matrix[[1]] <- sweep(health_matrix[[1]], 2, col_max, FUN = "/")
  
  column_order_health <- c("BMI", "WHR", "SBP", "TG", "non-HDL", "HbA1c", "CRP")
  health_matrix[[1]] <- health_matrix[[1]][, column_order_health]
  health_matrix[[2]] <- health_matrix[[2]][, column_order_health]
  
  colnames(health_matrix[[1]]) <- c("BMI", "WHR", "SBP", "TG", "non-HDL", "HbA1c", "CRP")
  
  col_titles <- c("Adiposity", "Adiposity", "Blood\npressure", "Lipids", "Lipids", "Glucose", "Inflam.")
  

  abx_matrix[[1]] <- abx_matrix[[1]][, column_order] 
  abx_matrix[[2]] <- abx_matrix[[2]][, column_order] 
  
  abx_matrix[[1]][is.na(abx_matrix[[1]])] <- 0
  abx_matrix[[2]][is.na(abx_matrix[[2]])] <- 1
  
  heatmap_colors <- colorRampPalette(c("dodgerblue3", "white", "mediumvioletred"))(5)
  
  
 HM <-  Heatmap(abx_matrix[[1]], 
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(abx_matrix[[2]][i, j] < 0.05 & abx_matrix[[2]][i, j] >= 0.01) {
              grid.text('*', x, y = y-unit(0.9, "mm"))
            } 
            if(abx_matrix[[2]][i, j] < 0.01) {
              grid.text('**', x, y = y-unit(0.9, "mm"))
            }},
          show_row_dend = F, show_column_dend = F, 
          cluster_rows = T, cluster_columns = F, 
          row_names_side = "left", 
          row_names_gp = gpar(fontsize = 6.5),  
          column_labels = rep(c("<1yr", "1-4yr", "4-8yr"), 3), 
          column_names_rot = 0,
          column_names_centered = T,
          column_names_gp = gpar(fontsize=6.3), 
          column_title_rot = 0, 
          column_title_gp = gpar(fontsize=7),
          column_split = c(rep("Clindamycin",3), rep("Flucloxacillin", 3), rep("Fluoroquinolones", 3)), 
          rect_gp = gpar(col = "darkgray", lwd = 0.3), 
          column_gap = unit(1.5, "mm"),
          cluster_column_slices = F,
          column_names_side = "top", 
          col = circlize::colorRamp2(c(-1.5, 0, 1.5), heatmap_colors[c(1, 3, 5)]), 
          name="Antibiotics", 
          width = unit(ncol(abx_matrix[[1]])*6.5, "mm"),
          heatmap_legend_param = list(
            direction = "horizontal",
            title_position = "lefttop",
            title_gp = gpar(fontsize = 6),  # title font size
            labels_gp = gpar(fontsize = 5),
            grid_height = unit(2, "mm")) ) +
  
  
  Heatmap(health_matrix[[1]], 
          cell_fun = function(j, i, x, y, w, h, fill) {
            if(health_matrix[[2]][i, j] < 0.05 & health_matrix[[2]][i, j] >= 0.01) {
              grid.text('*', x, y = y-unit(0.9, "mm"))
            } 
            if(health_matrix[[2]][i, j] < 0.01) {
              grid.text('**', x, y = y-unit(0.9, "mm"))
            }},
          show_row_dend = F, show_column_dend = F, 
          cluster_rows = F, 
          cluster_columns = F, 
          row_names_side = "left", 
          column_names_side = "top", 
          col = circlize::colorRamp2(c(-1, 0, 1), heatmap_colors[c(1, 3, 5)]), 
          column_names_gp = gpar(fontsize=6.3), 
          column_names_centered = T,
          name = "Cardiometabolic markers", 
          column_title_gp = gpar(fontsize=7),
          column_order = colnames(health_matrix[[1]]),
          column_split = factor(col_titles,  c("Adiposity", "Blood\npressure", "Lipids",  "Glucose", "Inflam.")), 
          cluster_column_slices = F,
          column_gap = unit(2.5, "mm"),
          rect_gp = gpar(col = "darkgray", lwd = 0.3), 
          column_names_rot = 0,
          width = unit(ncol(health_matrix[[1]])*8.5, "mm"),
          heatmap_legend_param = list(
            grid_height = unit(2, "mm"),
            direction = "horizontal",
            title_position = "lefttop",
            title_gp = gpar(fontsize = 6),  # title font size
            labels_gp = gpar(fontsize = 5)
            ) 
          )

 
# Save Heatmap 
    pdf(file="../../Revision_3/Figure4.pdf", width = 7.09, height=8.8 )
    draw(HM, heatmap_legend_side = "bottom", padding = unit(c(1, 0, 1, 0), "mm"))
    dev.off()
  
  
 
 
 
 
 
 
