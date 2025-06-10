rm(list=ls())
#################################
library(tidyverse)
library(cowplot) 
library(readxl)
library(icesTAF)
library(plyr)
library(dplyr)
library(writexl)
library(deSolve)
library(patchwork)
library(rstudioapi)
options(width = 10000)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
#################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

################################################################################
# LOAD DATA
keep_data     = readRDS('./SENSITIVITY_TH_RANGE.rds')
combined_data_aff = keep_data %>% pivot_longer(!c(range_tf_in, range_apop_in), names_to = "taxa", values_to = "affinity")

combined_data_aff_e = combined_data_aff %>% filter(taxa == 'Enterobacteriaceae')
combined_data_aff_b = combined_data_aff %>% filter(taxa == 'Bifidobacteriaceae')
combined_data_aff_bc = combined_data_aff %>% filter(taxa == 'Bacteroidaceae')
combined_data_aff_c = combined_data_aff %>% filter(taxa == 'Clostridiales')

# Define color vectors
c_vec_e = c("#ffb9c4", "#f20026")
c_vec_b = c("#7eb5fa", "#074fa8")
c_vec_bc = c("#aae8ef", "#2bc5d8")
c_vec_c = c("#ffd7d9", "#ff9aa1")

# Function to determine text color based on fill color brightness
get_text_color = function(fill_color) {
  rgb_values = col2rgb(fill_color)
  brightness = (rgb_values[1] * 299 + rgb_values[2] * 587 + rgb_values[3] * 114) / 1000
  ifelse(brightness > 128, "black", "white")
}

# Function to calculate the midpoint color
calculate_mid_color = function(low_color, high_color) {
  low_rgb = col2rgb(low_color)
  high_rgb = col2rgb(high_color)
  mid_rgb = (low_rgb + high_rgb) / 2
  rgb(mid_rgb[1], mid_rgb[2], mid_rgb[3], maxColorValue = 255)
}

# Apply color gradients and determine text color
apply_gradient = function(data, low_color, high_color) {
  mid_color = calculate_mid_color(low_color, high_color)
  data %>%
    mutate(fill_color = ifelse(is.na(affinity), "gray", colorRampPalette(c(low_color, mid_color, high_color))(100)[as.numeric(cut(affinity, breaks = 100))]),
           text_color = sapply(fill_color, get_text_color))
}

# Apply gradients to data
combined_data_aff_e = apply_gradient(combined_data_aff_e, c_vec_e[1], c_vec_e[2])
combined_data_aff_b = apply_gradient(combined_data_aff_b, c_vec_b[1], c_vec_b[2])
combined_data_aff_bc = apply_gradient(combined_data_aff_bc, c_vec_bc[1], c_vec_bc[2])
combined_data_aff_c = apply_gradient(combined_data_aff_c, c_vec_c[1], c_vec_c[2])

# Define plot functions
create_heatmap = function(data, title, low_color, high_color) {
  mid_color = calculate_mid_color(low_color, high_color)
  ggplot(data, aes(x = range_tf_in, y = range_apop_in, fill = fill_color)) +
    geom_tile(color = "white") +
    scale_fill_identity(na.value = "gray", name = "Affinity") +
    geom_text(aes(label = ifelse(is.na(affinity), "NA", round(affinity, 2)), color = text_color)) +
    scale_color_identity() +
    labs(title = title, x = expression(th[range]), y = expression(th[apop])) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      plot.title = element_text(hjust = 0.5, size = 14)
    ) +
    # Add bold borders for the specific box using annotate
    annotate("rect", xmin = 0.2, xmax = 0.3, ymin = 0.2, ymax = 0.3, 
             color = "black", size = 1.5, fill = NA)
}

# Create heatmaps
hm_e = create_heatmap(combined_data_aff_e, "Enterobacteriaceae", c_vec_e[1], c_vec_e[2])
hm_b = create_heatmap(combined_data_aff_b, "Bifidobacteriaceae", c_vec_b[1], c_vec_b[2])
hm_bc = create_heatmap(combined_data_aff_bc, "Bacteroidaceae", c_vec_bc[1], c_vec_bc[2])
hm_c = create_heatmap(combined_data_aff_c, "Clostridiales", c_vec_c[1], c_vec_c[2])

graphics.off()
row0 = plot_grid(hm_e,hm_b,hm_bc,hm_c,nrow = 2,rel_widths = c(1,1,1,1),labels=c('A','B','C','D'), label_fontface = 'bold')
png(file =paste0("./FIGURE_S9.png"),    # The directory you want to save the file in
    width     = 14,
    height    = 13,
    units     = "in",
    res       = 300)
print(row0)
dev.off()

