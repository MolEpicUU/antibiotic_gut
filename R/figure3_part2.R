#library(MicrobiotaProcess)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(patchwork)

source('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__atb_scripts/order_species_by_taxa.R')

anno <- rio::import('/Users/gabba126/Documents/PhD_projects/Microbiome/Taxonomy/Taxonomy_CHAMP.tsv')
rownames(anno) <- anno[, 1]
anno <- anno[, -1]
anno <- anno[, rev(colnames(anno))]

res <- rio::import('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__tables/species_results_for_fig2.tsv')
res <- res[which(res$model == "full model"), ]
res$species <- res$outcome
species <- unique(res$species)
res$beta[which(res$qvalue >= 0.05)] <- NA
res$time <- factor(gsub(" |ear|ears", "", res$period), levels = rev(c("<1y", "1-4y", "4-8y")))
res$antibiotic[which(res$antibiotic == "Penicillins extended spectrum")] <- "Penicillins ES"
res$abx <- res$antibiotic

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

users <- rio::import("/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_2/revision__work/antibiotics_number_of_users.tsv")
users$antibiotic[which(users$antibiotic == "Penicillin extended spectrum")] <- "Penicillins ES"
res$users <- users$proportion_users[match(paste(res$antibiotic, res$time), paste(users$antibiotic, users$period))]
res$abx_order <- users$antibiotic_order[match(res$antibiotic, users$antibiotic)]
res$abx <- factor(res$abx, levels = unique(res$abx[order(res$abx_order)]))
summary(res$user)

# Define proportion of users labels based on the distribution
plot <- ggplot(res, aes(users)) + geom_density(fill = "grey") + theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + geom_vline(xintercept = 0.15, linetype = "dashed")
print(plot)

res$users[which(res$users > 0.15)] <- 0.15
user_breaks <- c(0, 0.05, 0.10, 0.15)
user_labels <- c("0.00", 0.05, 0.10, ">0.15")
summary(res$beta)

# Define beta labels based on the distribution
plot <- ggplot(res, aes(beta)) + geom_density(fill = "grey") + theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + geom_vline(xintercept = c(-1, 1), linetype = "dashed")

print(plot)

res$beta[which(res$beta < -1)] <- -1 
res$beta[which(res$beta > 1)] <- 1
beta_breaks <- c(-1, -0.5, max(res$beta[which(res$qvalue < 0.05 & res$beta < 0)]), min(res$beta[which(res$qvalue < 0.05 & res$beta > 0)]), 0.5, 1)
beta_labels <- c("<-1.00", "-0.50", "-0.04       ", "       0.04", "0.50", ">1.00")

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

plot1 <- ggplot(res, aes(x = species, y = time, fill = beta)) + theme_grey(base_size = 7) + facet_wrap(~abx, ncol = 1, strip.position = "left") + geom_raster() + labs(x = NULL, y = NULL, fill = "            Regression coefficient (q-value < 0.05)            ") + scale_fill_gradient2(breaks = beta_breaks, labels = beta_labels, high = "red4", low = "navyblue", na.value = "white") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(0.5, "lines"), panel.border = element_rect(fill = NA), plot.margin = margin(t=10)) + theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5, size = 7), legend.key.width = unit(10, "mm"), legend.key.height = unit(2.5, "mm")) + geom_vline(xintercept = lines, linetype = "dashed", linewidth = 0.2)
plot2 <- ggplot(res[which(!is.na(res$antibiotic)), ], aes(x = 1, y = time, fill = users)) + theme_grey(base_size = 7) + facet_wrap(~abx, ncol = 1, strip.position = "left") + geom_raster() + labs(x = NULL, y = NULL, fill = "Proportion users") + scale_fill_gradient(breaks = user_breaks, labels = user_labels, low = "whitesmoke", high = "darkgreen", limits = c(0, max(user_breaks))) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7), axis.text.y = element_text(size = 7, color = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(0.5, "lines")) + theme(legend.position = "bottom", legend.title.position = "top", legend.title = element_text(hjust = 0.5), legend.key.width = unit(4, "mm"), legend.key.height = unit(2.5, "mm"))

plot3 <- ggplot(res[which(!is.na(res$antibiotic)), ], aes(x = species, y = "Order", fill = order, label = order_label)) + theme_grey(base_size = 7) + geom_raster() + geom_text(size = 7/.pt) + labs(x = NULL, y = NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank(),  axis.text.y = element_blank() , axis.ticks.y.left = element_blank(), axis.title.y = element_text(angle = 0, hjust = 1, size = 7, color = "black", margin = margin(r=-1)), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + guides(fill = guide_legend(nrow = 1, label.position = "bottom"), color = "none", linewidth = "none") + guides(fill = "none") + scale_fill_discrete(na.value = "white")
plot4 <- ggplot(res[which(!is.na(res$antibiotic)), ], aes(x = species, y = "Phylum", fill = phylum, label = phylum_label)) + theme_grey(base_size = 7) + geom_raster() + geom_text(size = 7/.pt) + labs(x = NULL, y = NULL) + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + theme(strip.background = element_blank(), strip.text.y.left = element_blank(),  axis.text.y = element_blank() , axis.ticks.y.left = element_blank(), axis.title.y = element_text(angle = 0, hjust = 1, size = 7, color = "black", margin = margin(r=-1)), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), strip.placement = "outside", panel.spacing = unit(1, "lines")) + guides(fill = guide_legend(nrow = 1, label.position = "bottom"), color = "none", linewidth = "none") + guides(fill = "none") + scale_fill_discrete(na.value = "white")

order_title <- ggplot(data.frame(x1="Order"), aes(y="Order", x = 1)) + geom_raster(fill = NA) +  theme_void() + theme(axis.text.y = element_text(size=7, margin = margin(r=-23), hjust = 1))
phylum_title <- ggplot(data.frame(x1="Phylum"), aes(y="Phylum", x = 1)) + geom_raster(fill = NA) +  theme_void() + theme(axis.text.y = element_text(size=7, margin = margin(r=-23), hjust = 1))

plot <- plot2 + plot1 + order_title + plot3 + phylum_title + plot4 + plot_layout(ncol = 2, widths = c(1, 100), heights = c(50, 1, 1), guides = "collect") & theme(legend.position = "bottom", plot.margin = margin(t=0.2, 0, 0, 0, unit = "mm"), legend.text = element_text(size = 7))
#plot

ggsave('/Users/gabba126/Documents/PhD_projects/3.Antibiotic/revision_NatMed/Revision_3/Figure3.pdf', plot, width = 180, height = 160, unit = "mm")
