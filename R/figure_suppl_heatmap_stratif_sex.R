
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(patchwork)

rm(list=ls())

source('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__atb_scripts/order_species_by_taxa.R')
load("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__work/antibiotic_order_figures")
order_atb <- gsub("Penicillin ES", "Penicillins ES", order_atb)

# FEMALES


anno <- rio::import('/Users/gabba126/Documents/PhD_projects/Microbiome/Taxonomy/Taxonomy_CHAMP.tsv')
res <- rio::import('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/meta_species__sa_stratified_allantibiotics_female.tsv')
res <- merge(res, anno[, c("MGS", "species")], by.x = "outcome", by.y = "MGS")
setnames(res, "q.value", "qvalue")

rownames(anno) <- anno[, 1]
anno <- anno[, -1]
anno <- anno[, rev(colnames(anno))]


res <- res[which(res$model == "female"), ]
species <- unique(res$species)
res$beta[which(res$qvalue >= 0.05)] <- NA

res$exposure <- stringr::str_extract(res$exposure, "Class_.*$")
res$period <- stringr::str_extract(res$exposure, "\\d.*yr$")
res$time <- factor(res$period, levels = rev(c('1yr', '1_4yr', '4_8yr')), labels = rev(c("<1y", "1-4y", "4-8y")))

x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                       "cephalosporins","macrolides",
                       "lincosamides","TCLs","FQs","SMZTMP","NIT"))

lab = c("Penicillin V", "Flucloxacillin", "Penicillins ES",  "Amox-clav", 
        "Cephalosporins",
        "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")

res$abx <- factor(gsub("_\\d.*$", "", res$exposure), x, lab)
res$abx <- factor(res$abx, order_atb)
res$antibiotic <- factor(res$abx, order_atb)


tree <- phyloseq::tax_table(as.matrix(anno[which(anno$species %in% species), ]))

####
tax_df <- as.data.frame(anno, stringsAsFactors = FALSE)
tax_df <- tax_df[tax_df$species %in% species, , drop = FALSE]
species_ordered <- order_species_by_taxa(tax_df, species,
                                         tax_levels = c("kingdom","phylum","class","order","family","genus","species"),
                                         species_col = "species")
####

res$species <- factor(res$species, levels = species_ordered)
res <- res[order(res$species), ]


### users_female ####
users <- rio::import("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__work/count_antibiotics_nonzeros_bysex_allcohorts.tsv")
users_male <- users[, grep("abx|_male", colnames(users))]
users <- users[, grep("abx|_female", colnames(users))]
users$abx <- gsub("Pen\\.", "Penicillin", users$abx)
users$abx <- gsub("Marcrolides", "Macrolides", users$abx)
users$abx <- gsub("Penicillin ES", "Penicillins ES", users$abx)

users[users<=5] <- 0
users$users <- rowSums(users[, -1], na.rm = T)/7591
users$abx <- gsub("0_1yr", "1yr", users$abx)
users$period <- stringr::str_extract(users$abx, "\\d.*yr")
users$abx <- trimws(gsub("\\d.*yr", "", users$abx))


res <- merge(res, users[, c("abx", "period", "users")], by = c("abx", "period"))
all(order_atb %in% res$abx)
all(order_atb %in% users$abx)

res$abx <- factor(res$abx, order_atb)
users$abx <- factor(users$abx, order_atb)
res$antibiotic <- factor(res$antibiotic, order_atb)

summary(res$user)
plot <- ggplot(res, aes(users)) + geom_density(fill = "grey") + theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0.15, linetype = "dashed")
# plot

res$users[which(res$users > 0.15)] <- 0.15
user_breaks <- c(0, 0.05, 0.10, 0.15)
user_labels <- c("0.00", 0.05, 0.10, ">0.15")

max(res$beta[which(res$qvalue < 0.05 & res$beta < 0)])
min(res$beta[which(res$qvalue < 0.05 & res$beta < 0)])

min(res$beta[which(res$qvalue < 0.05 & res$beta > 0)])
max(res$beta[which(res$qvalue < 0.05 & res$beta > 0)])


res$beta[which(res$beta < -1)] <- -1 
res$beta[which(res$beta > 1)] <- 1

beta_breaks <- c(-1, -0.5, max(res$beta[which(res$qvalue < 0.05 & res$beta < 0)]), min(res$beta[which(res$qvalue < 0.05 & res$beta > 0)]), 0.5, 1)
beta_labels <- c("<-1", "-0.5", "-0.04      ", "      0.04", "0.5", ">1")

anno <- anno[match(species_ordered, anno$species), ]
anno$order[which(!anno$phylum == "Bacillota A")] <- NA
phylum <- as.data.frame(table(anno$phylum))
anno$phylum[which(anno$phylum %in% phylum$Var1[which(phylum$Freq < 50)])] <- "Other"
order <- as.data.frame(table(anno$order))
anno$order[which(anno$order %in% order$Var1[which(order$Freq < 50)])] <- "Other"

res$order <- anno$order[match(res$species, anno$species)]
res$phylum <- anno$phylum[match(res$species, anno$species)]

phylum_lines <- which(!duplicated(anno$phylum))
order_lines <- which(!duplicated(anno$order))
lines <- c(phylum_lines[-1], order_lines[-1])
phylum_lines <- diff(c(phylum_lines, nrow(anno))) / 2 + phylum_lines
order_lines <- diff(c(order_lines, nrow(anno))) / 2 + order_lines
res$phylum_label <- res$phylum
res$phylum_label[which(!res$species %in% anno$species[phylum_lines])] <- NA
res$phylum_label[which(duplicated(res$phylum_label))] <- NA
res$order_label <- res$order
res$order_label[which(!res$species %in% anno$species[order_lines])] <- NA
res$order_label[which(duplicated(res$order_label))] <- NA
res$order_label[which(res$order_label == "Christensenellales")] <- "Christens."
res$phylum_label[which(res$phylum_label == "Actinomycetota")] <- "Actino."
res$phylum_label[which(res$phylum_label == "Bacteroidota")] <- "Bacteroid."

plot_f1 <- ggplot(res, aes(x = species, y = time, fill = beta)) + theme_grey(base_size = 6) + facet_wrap(~abx, ncol = 1, strip.position = "left") + geom_raster() + labs(x = NULL, y = NULL, fill = "            Regression coefficient (q-value < 0.05)            ") + scale_fill_gradient2(breaks = beta_breaks, labels = beta_labels, high = "red4", low = "navyblue", na.value = "white") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(0.2, "lines"), panel.border = element_rect(fill = NA)) + theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5), legend.key.width = unit(10, "mm"), legend.key.height = unit(2.5, "mm")) + geom_vline(xintercept = lines, linetype = "dashed", linewidth = 0.2)
plot_f2 <- ggplot(res[which(!is.na(res$abx)), ], aes(x = 1, y = time, fill = users)) + theme_grey(base_size = 6) + facet_wrap(~abx, ncol = 1, strip.position = "left") + geom_raster() + labs(x = NULL, y = NULL, fill = "Proportion users") + scale_fill_gradient(breaks = user_breaks, labels = user_labels, low = "whitesmoke", high = "darkgreen", limits = c(0, max(user_breaks))) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6), axis.text.y = element_text(size = 6, color = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(0.2, "lines")) + theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5), legend.key.width = unit(4, "mm"), legend.key.height = unit(2.5, "mm"))

plot_f3 <- ggplot(res[which(!is.na(res$abx)), ], aes(x = species, y = 1, fill = order, label = order_label)) + theme_grey(base_size = 6) + geom_raster() + geom_text(size = 6/.pt) + labs(x = NULL, y = NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + guides(fill = guide_legend(nrow = 1, label.position = "bottom"), color = "none", linewidth = "none") + guides(fill = "none") + scale_fill_discrete(na.value = "white")
plot_f4 <- ggplot(res[which(!is.na(res$abx)), ], aes(x = species, y = 1, fill = phylum, label = phylum_label)) + theme_grey(base_size = 6) + geom_raster() + geom_text(size = 6/.pt) + labs(x = NULL, y = NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + guides(fill = guide_legend(nrow = 1, label.position = "bottom"), color = "none", linewidth = "none") + guides(fill = "none") + scale_fill_discrete(na.value = "white")

order_title <- ggplot(data.frame(x1="Order"), aes(y="Order", x = 1)) + geom_raster(fill = NA) +  theme_void() + theme(axis.text.y = element_text(size=6, margin = margin(r=-23), hjust = 1, face = "bold"))
phylum_title <- ggplot(data.frame(x1="Phylum"), aes(y="Phylum", x = 1)) + geom_raster(fill = NA) +  theme_void() + theme(axis.text.y = element_text(size=6, margin = margin(r=-23), hjust = 1, face = "bold"))

# MALES ####

anno_male <- rio::import('/Users/gabba126/Documents/PhD_projects/Microbiome/Taxonomy/Taxonomy_CHAMP.tsv')
res <- rio::import('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__results/meta_species__sa_stratified_allantibiotics_male.tsv')
res <- merge(res, anno_male[, c("MGS", "species")], by.x = "outcome", by.y = "MGS")
setnames(res, "q.value", "qvalue")

res <- res[which(res$model == "male"), ]

species <- unique(res$species)
res$beta[which(res$qvalue >= 0.05)] <- NA

res$exposure <- stringr::str_extract(res$exposure, "Class_.*$")
res$period <- stringr::str_extract(res$exposure, "\\d.*yr$")
res$time <- factor(res$period, levels = rev(c('1yr', '1_4yr', '4_8yr')), labels = rev(c("<1y", "1-4y", "4-8y")))

x <- paste0("Class_",c("Peni_BetaS","Peni_BetaR","Peni_Ext","Peni_Comb",
                       "cephalosporins","macrolides",
                       "lincosamides","TCLs","FQs","SMZTMP","NIT"))

lab = c("Penicillin V", "Flucloxacillin", "Penicillins ES",  "Amox-clav", 
        "Cephalosporins",
        "Macrolides", "Clindamycin", "Tetracyclines", "Fluoroquinolones","SMZ-TMP","Nitrofurantoin")

res$abx <- factor(gsub("_\\d.*$", "", res$exposure), x, lab)
res$abx <- factor(res$abx, order_atb)
res$antibiotic <- factor(res$abx, order_atb)

res$species <- factor(res$species, levels = species_ordered)
res <- res[order(res$species), ]


# Users_male ####
users <- users_male
users$abx <- gsub("Pen\\.", "Penicillin", users$abx)
users$abx <- gsub("Marcrolides", "Macrolides", users$abx)
users$abx <- gsub("Penicillin ES", "Penicillins ES", users$abx)

users[users<=5] <- 0
users$users <- rowSums(users[, -1], na.rm = T)/7383
users$abx <- gsub("0_1yr", "1yr", users$abx)
users$period <- stringr::str_extract(users$abx, "\\d.*yr")
users$abx <- trimws(gsub("\\d.*yr", "", users$abx))


res <- merge(res, users[, c("abx", "period", "users")], by = c("abx", "period"))
all(order_atb %in% res$abx)
all(order_atb %in% users$abx)

res$abx <- factor(res$abx, order_atb)
users$abx <- factor(users$abx, order_atb)
res$antibiotic <- factor(res$antibiotic, order_atb)


res$users[which(res$users > 0.15)] <- 0.15
user_breaks <- c(0, 0.05, 0.10, 0.15)
user_labels <- c("0.00", 0.05, 0.10, ">0.15")
summary(res$beta)
plot <- ggplot(res, aes(beta)) + geom_density(fill = "grey") + theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + geom_vline(xintercept = c(-1, 1), linetype = "dashed")


res$beta[which(res$beta < -1)] <- -1 
res$beta[which(res$beta > 1)] <- 1

res$order <- anno$order[match(res$species, anno$species)]
res$phylum <- anno$phylum[match(res$species, anno$species)]


res$phylum_label <- res$phylum
res$phylum_label[which(!res$species %in% anno$species[phylum_lines])] <- NA
res$phylum_label[which(duplicated(res$phylum_label))] <- NA
res$order_label <- res$order
res$order_label[which(!res$species %in% anno$species[order_lines])] <- NA
res$order_label[which(duplicated(res$order_label))] <- NA
res$order_label[which(res$order_label == "Christensenellales")] <- "Christens."
res$phylum_label[which(res$phylum_label == "Actinomycetota")] <- "Actino."
res$phylum_label[which(res$phylum_label == "Bacteroidota")] <- "Bacteroid."

plot_m1 <- ggplot(res, aes(x = species, y = time, fill = beta)) + theme_grey(base_size = 6) + facet_wrap(~abx, ncol = 1, strip.position = "left") + geom_raster() + labs(x = NULL, y = NULL, fill = "            Regression coefficient (q-value < 0.05)            ") + scale_fill_gradient2(breaks = beta_breaks, labels = beta_labels, high = "red4", low = "navyblue", na.value = "white") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(0.2, "lines"), panel.border = element_rect(fill = NA)) + theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5), legend.key.width = unit(10, "mm"), legend.key.height = unit(2.5, "mm")) + geom_vline(xintercept = lines, linetype = "dashed", linewidth = 0.2)
plot_m2 <- ggplot(res[which(!is.na(res$antibiotic)), ], aes(x = 1, y = time, fill = users)) + theme_grey(base_size = 6) + facet_wrap(~abx, ncol = 1, strip.position = "left") + geom_raster() + labs(x = NULL, y = NULL, fill = "Proportion users") + scale_fill_gradient(breaks = user_breaks, labels = user_labels, low = "whitesmoke", high = "darkgreen", limits = c(0, max(user_breaks))) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6), axis.text.y = element_text(size = 6, color = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(0.2, "lines")) + theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5), legend.key.width = unit(4, "mm"), legend.key.height = unit(2.5, "mm"))

plot_m3 <- ggplot(res[which(!is.na(res$antibiotic)), ], aes(x = species, y = "Order", fill = order, label = order_label)) + theme_grey(base_size = 6) + geom_raster() + geom_text(size = 6/.pt) + labs(x = NULL, y = NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank(),  axis.text.y = element_blank() , axis.ticks.y.left = element_blank(), axis.title.y = element_text(angle = 0, hjust = 1, size = 6, color = "black", margin = margin(r=-1)), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + guides(fill = guide_legend(nrow = 1, label.position = "bottom"), color = "none", linewidth = "none") + guides(fill = "none") + scale_fill_discrete(na.value = "white")
plot_m4 <- ggplot(res[which(!is.na(res$antibiotic)), ], aes(x = species, y = "Phylum", fill = phylum, label = phylum_label)) + theme_grey(base_size = 6) + geom_raster() + geom_text(size = 6/.pt) + labs(x = NULL, y = NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank(),  axis.text.y = element_blank() , axis.ticks.y.left = element_blank(), axis.title.y = element_text(angle = 0, hjust = 1, size = 6, color = "black", margin = margin(r=-1)), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + guides(fill = guide_legend(nrow = 1, label.position = "bottom"), color = "none", linewidth = "none") + guides(fill = "none") + scale_fill_discrete(na.value = "white")

order_title <- ggplot(data.frame(x1="Order"), aes(y="Order", x = 1)) + geom_raster(fill = NA) +  theme_void() + theme(axis.text.y = element_text(size=6, margin = margin(r=-23), hjust = 1))
phylum_title <- ggplot(data.frame(x1="Phylum"), aes(y="Phylum", x = 1)) + geom_raster(fill = NA) +  theme_void() + theme(axis.text.y = element_text(size=6, margin = margin(r=-23), hjust = 1))

p_f2_tagged <- plot_f2 + labs(tag = "a")
p_m2_tagged <- plot_m2 + labs(tag = "b")

p_f2_tagged <- p_f2_tagged + theme(plot.tag.position = c(0.15, .99),
                                   plot.tag = element_text(size = 7, face = "bold"))
p_m2_tagged <- p_m2_tagged + theme(plot.tag.position = c(0.15, .999),
                                   plot.tag = element_text(size = 7, face = "bold"))


plot <- p_f2_tagged + plot_f1 + plot_spacer() + plot_spacer () + p_m2_tagged  + plot_m1 + order_title +  plot_m3 + phylum_title + plot_m4 + plot_layout(ncol = 2, widths = c(1, 100), heights = c(50, 7, 50, 1.8, 1.8), guides = "collect") & theme(legend.position = "bottom", plot.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 5))
plot

ggsave('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_3/SupplFigures/fig_suppl_stratif_sex.pdf', plot, width = 180, height = 180, unit = "mm", dpi = 400)

