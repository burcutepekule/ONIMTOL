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

###############################################################################
keep_data = readRDS('./SENSITIVITY_MCELL_TIMING.rds')
my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
names(my_colors) = c('Bacteroidaceae', 'Bifidobacteriaceae', 'Clostridiales', 'Enterobacteriaceae')

keep_data$mcells_open_days_relative=keep_data$mcells_open_days-129
# keep_data = keep_data %>% filter(mcells_open_days_relative>-40)

combined_data_aff = keep_data[c(2,3,4,5,11)] # mcells_open_days_relative column is the last
combined_data_aff = combined_data_aff %>% pivot_longer(!mcells_open_days_relative, names_to = "taxa", values_to = "affinity")

p1 = ggplot(combined_data_aff, aes(x = mcells_open_days_relative, y = affinity, color = taxa)) +
  geom_line(size = 0.8) +  theme_minimal() +
  labs(title = "Average eSlgA Affinity", x = "M cell opening day relative\n to Baseline (DOL 129)", y = "") +
  scale_color_manual(values =my_colors) +
  theme(legend.title = element_blank(), plot.title = element_text(size = 10)) +
  scale_y_log10()+
  scale_x_continuous(breaks = seq(-90,20,10))

p2 = ggplot(keep_data, aes(x = mcells_open_days_relative, y = cumul_entero_load)) +
  geom_line(size = 0.8) +  theme_minimal() +
  labs(title = "Cumulative Enterobacteriaceae load\nin the gut lumen",
       x = "M cell opening day relative\n to Baseline (DOL 129)", 
       y = "(Cells/gLC) x 1e-11") + theme(legend.title = element_blank(), plot.title = element_text(size = 10)) 

p3 = ggplot(keep_data, aes(x = mcells_open_days_relative, y = cumul_entero_sampledPP)) +
  geom_line(size = 0.8) +  theme_minimal() +
  labs(title = "Cumulative Enterobacteriaceae load\nin the GALT inductive sites",
       x = "M cell opening day relative\n to Baseline (DOL 129)", 
       y = "(Cells/gLC) x 1e-11") + theme(legend.title = element_blank(), plot.title = element_text(size = 10)) 

graphics.off()
row0 = plot_grid(p1,p2,p3,nrow = 1,rel_widths = c(1.6,1,1),labels=c('A)','B)','C)'), label_fontface = 'bold')
png(file =paste0("./FIGURE_S10.png"),    # The directory you want to save the file in
    width     = 12,
    height    = 3,
    units     = "in",
    res       = 300)
print(row0)
dev.off()

