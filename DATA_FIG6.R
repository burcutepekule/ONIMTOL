#################################
rm(list=ls())
library(tidyverse)
library(cowplot)  # this is a plotting library but I'm unsure if you want to retain it
library(readxl)
library(plyr)
library(dplyr)
library(writexl)
library(deSolve)
library(ggplot2)
library(reshape2)
options(width = 10000)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

############ LOAD THE DATA WITH MCELL EXPERIMENT
load('./MCELL_EXPERIMENT.RDATA')

keep_res_level_1 = keep_res %>% dplyr::filter(alpha_b==combinations_df_2[1,]$alpha_b & alpha_bc==combinations_df_2[1,]$alpha_bc & alpha_c==combinations_df_2[1,]$alpha_c)
keep_res_level_2 = keep_res %>% dplyr::filter(alpha_b==combinations_df_2[2,]$alpha_b & alpha_bc==combinations_df_2[2,]$alpha_bc & alpha_c==combinations_df_2[2,]$alpha_c)
keep_res_level_3 = keep_res %>% dplyr::filter(alpha_b==combinations_df_2[3,]$alpha_b & alpha_bc==combinations_df_2[3,]$alpha_bc & alpha_c==combinations_df_2[3,]$alpha_c)

# Define your min, mid, and max values for the color scale
min_val = min(c(keep_res$aff_b_final, keep_res$aff_bc_final, keep_res$aff_c_final),na.rm = TRUE)
mid_val = 1  # Assuming 0 is your midpoint for gray
max_val = max(keep_res$aff_e_final)  # Setting the maximum value to 26 as requested

keep_res_in = keep_res_level_1
source('./misc/PLOT_MIGA_MCELL_ALL_HEATMAPS.R')
plot_e_level_1  = plot_e
plot_b_level_1  = plot_b
plot_bc_level_1 = plot_bc
plot_c_level_1  = plot_c
plot_ne_level_1 = plot_ne

keep_res_in = keep_res_level_2
source('./misc/PLOT_MIGA_MCELL_ALL_HEATMAPS.R')
plot_e_level_2  = plot_e
plot_b_level_2  = plot_b
plot_bc_level_2 = plot_bc
plot_c_level_2  = plot_c
plot_ne_level_2 = plot_ne

keep_res_in = keep_res_level_3
source('./misc/PLOT_MIGA_MCELL_ALL_HEATMAPS.R')
plot_e_level_3  = plot_e
plot_b_level_3  = plot_b
plot_bc_level_3 = plot_bc
plot_c_level_3  = plot_c
plot_ne_level_3 = plot_ne



plot_e_level_1  = plot_e_level_1 + theme(plot.title = element_text(face = "italic"))
plot_b_level_1  = plot_b_level_1 + theme(plot.title = element_text(face = "italic"))
plot_bc_level_1 = plot_bc_level_1 + theme(plot.title = element_text(face = "italic"))
plot_c_level_1  = plot_c_level_1 + theme(plot.title = element_text(face = "italic"))

row_level_1 = cowplot::plot_grid(plot_e_level_1,plot_b_level_1,plot_bc_level_1,plot_c_level_1,ncol = 2, nrow =2,align='h', rel_widths = c(1,1,1,1),labels=c('A)','B)','C)','D)'),label_fontface = "bold")
row_level_2 = cowplot::plot_grid(plot_e_level_2,plot_b_level_2,plot_bc_level_2,plot_c_level_2,ncol = 2, nrow =2,align='h', rel_widths = c(1,1,1,1),labels=c('A)','B)','C)','D)'),label_fontface = "bold")
row_level_3 = cowplot::plot_grid(plot_e_level_3,plot_b_level_3,plot_bc_level_3,plot_c_level_3,ncol = 2, nrow =2,align='h', rel_widths = c(1,1,1,1),labels=c('A)','B)','C)','D)'),label_fontface = "bold")


graphics.off()
png(file =paste0("./FIGURE_6.png"),
    width     = 8,
    height    = 7,
    units     = "in",
    res       = 300)
row_level_1
dev.off()

