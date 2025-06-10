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
library(ggtext)
options(width = 10000)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
#################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')

keep_data = readRDS('./SENSITIVITY_TAUDELTA.rds')
combined_data_aff = keep_data[c(1,2,3,4,5)] 
combined_data_aff = combined_data_aff %>% pivot_longer(!taudelta_in, names_to = "taxa", values_to = "affinity")

p1 = ggplot(combined_data_aff, aes(x = taudelta_in, y = affinity, color = taxa)) +
  geom_line(size = 0.8) +  theme_minimal() +
  labs(
    title = "Average eSlgA Affinity",
    x = "τ<sup>δ</sup>",
    y = ""
  ) +
  scale_color_manual(values = my_colors) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 10),
    axis.title.x = ggtext::element_markdown(size = 10)
  ) +
  scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = 0.42, col = 'red', linetype = 'dashed')


keep_data = readRDS('./SENSITIVITY_SMALLCN.rds')
combined_data_aff = keep_data[c(1,2,3,4,5)] 
combined_data_aff = combined_data_aff %>% pivot_longer(!cn_in, names_to = "taxa", values_to = "affinity")

p2 = ggplot(combined_data_aff, aes(x = cn_in, y = affinity, color = taxa)) +
  geom_line(size = 0.8) +  theme_minimal() +
  labs(
    title = "Average eSlgA Affinity",
    x = "c<sub>n</sub>",
    y = ""
  ) +
  scale_color_manual(values = my_colors) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 10),
    axis.title.x = ggtext::element_markdown(size = 10)
  )  +
  scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = 6.34*1e-4, col = 'red', linetype = 'dashed')


keep_data = readRDS('./SENSITIVITY_BIGCN.rds')
combined_data_aff = keep_data[c(1,2,3,4,5)] 
combined_data_aff = combined_data_aff %>% pivot_longer(!C_n_in, names_to = "taxa", values_to = "affinity")

p3 = ggplot(combined_data_aff, aes(x = C_n_in, y = affinity, color = taxa)) +
  geom_line(size = 0.8) +  theme_minimal() +
  labs(
    title = "Average eSlgA Affinity",
    x = "C<sub>n</sub>",
    y = ""
  ) +
  scale_color_manual(values = my_colors) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 10),
    axis.title.x = ggtext::element_markdown(size = 10)
  ) +
  scale_x_log10() +
  scale_y_log10() +
  geom_vline(xintercept = 0.0344, col = 'red', linetype = 'dashed')


p1 = p1 + theme(legend.position = 'none') 
p3 = p3 + theme(legend.position = 'none') 

row = plot_grid(p1,p3,p2,align='h',ncol = 3, rel_widths = c(1,1,1.5),labels=c('A)','B)','C)'))

graphics.off()
png(file =paste0("./FIGURE_S11.png"),    # The directory you want to save the file in
    width     = 12,
    height    = 3,
    units     = "in",
    res       = 300)
row
dev.off()
