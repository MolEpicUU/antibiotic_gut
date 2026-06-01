# This script will create the plot showing the FDA results. 

library(ggplot2)
library(patchwork)

abx <- rio::import("/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/atc_abx_classes.csv")
abx$abx_class[6] <- "Penicillins ES"

files <- list.files("/proj/nobackup/sens2019512/data_for_Koen/", pattern = "shannon", full.names = TRUE)
data <- lapply(files, function(x) {
  data <- rio::import(x)
  data$abx <- gsub("^.*_|.dta", "", x)
  data
})
data <- do.call(rbind, data)
data$lower <- data$mean - 1.96 * data$se
data$upper <- data$mean + 1.96 * data$se
data$abx <- abx$abx_class[match(data$abx, abx$ATC)]
data$abx <- factor(data$abx, levels = abx$abx_class)

plot1 <- ggplot(data, aes(x = month, y = mean, ymin = lower, ymax = upper, color = abx, fill = abx)) + geom_line() + geom_ribbon(alpha = 0.3, color = NA) + theme_bw() + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 96, length.out = 9), labels = 0:8) + geom_hline(yintercept = 0) + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) + labs(x = NULL, y = "Shannon index", fill = "Antibiotic", color = "Antibiotic") + ylim(-0.55, 0.2) + facet_wrap(~abx, nrow = 1) + theme(strip.background = element_rect(fill = NA)) + guides(color = "none", fill = "none")

files <- list.files("/proj/nobackup/sens2019512/data_for_Koen/", pattern = "richness", full.names = TRUE)
data <- lapply(files, function(x) {
  data <- rio::import(x)
  data$abx <- gsub("^.*_|.dta", "", x)
  data
})
data <- do.call(rbind, data)
data$lower <- data$mean - 1.96 * data$se
data$upper <- data$mean + 1.96 * data$se
data$abx <- abx$abx_class[match(data$abx, abx$ATC)]
data$abx <- factor(data$abx, levels = abx$abx_class)

plot2 <- ggplot(data, aes(x = month, y = mean, ymin = lower, ymax = upper, color = abx, fill = abx)) + geom_line() + geom_ribbon(alpha = 0.3, color = NA) + theme_bw() + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 96, length.out = 9), labels = 0:8) + geom_hline(yintercept = 0) + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) + labs(x = NULL, y = "Species richness", fill = "Antibiotic", color = "Antibiotic") + ylim(-110, 40) + facet_wrap(~abx, nrow = 1) + theme(strip.background = element_rect(fill = NA)) + guides(color = "none", fill = "none")

files <- list.files("/proj/nobackup/sens2019512/data_for_Koen/", pattern = "invsimpson", full.names = TRUE)
data <- lapply(files, function(x) {
  data <- rio::import(x)
  data$abx <- gsub("^.*_|.dta", "", x)
  data
})
data <- do.call(rbind, data)
data$lower <- data$mean - 1.96 * data$se
data$upper <- data$mean + 1.96 * data$se
data$abx <- abx$abx_class[match(data$abx, abx$ATC)]
data$abx <- factor(data$abx, levels = abx$abx_class)

plot3 <- ggplot(data, aes(x = month, y = mean, ymin = lower, ymax = upper, color = abx, fill = abx)) + geom_line() + geom_ribbon(alpha = 0.3, color = NA) + theme_bw() + scale_x_continuous(expand = c(0, 0), breaks = seq(0, 96, length.out = 9), labels = 0:8) + geom_hline(yintercept = 0) + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank()) + labs(x = "Years since antibiotic course", y = "Inv. Simpson", fill = "Antibiotic", color = "Antibiotic") + ylim(-11, 4) + facet_wrap(~abx, nrow = 1) + theme(strip.background = element_rect(fill = NA)) + guides(color = "none", fill = "none")

plot <- plot1 + plot2 + plot3 + plot_layout(guides = "collect", ncol = 1)

ggsave("/proj/sens2019512/nobackup/wharf/koede543/koede543-sens2019512/transfer/fda.pdf", plot, width = 8.3, height = 11.7)