#################################
rm(list=ls())
library(tidyverse)
library(cowplot)  # this is a plotting library but I'm unsure if you want to retain it
library(readxl)
library(plyr)
library(dplyr)
library(writexl)
library(deSolve)
options(width = 10000)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
my_colors   = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
spacer_grob = grid::nullGrob()

########################### BASELINE
load('./SIM_BASELINE.RDATA')
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
source('./misc/PLOT_ABUNDANCE_FIG5.R') #p_abs, p_rel, p_coating, p_micro
p_abs     = p_abs +scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_pcr = p_pcr + geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_pcr = p_pcr + geom_hline(yintercept = 0.50, color = "red", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_aff_log = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_coating = p_coating + scale_x_discrete(labels=c("C","BC","B","E"))
p_coating = p_coating + geom_hline(yintercept = 0.6, linetype = "dashed", color = "black", size = 1)
p_rel     = p_rel + theme(legend.position = "none")
row_0     = cowplot::plot_grid(p_rel,p_abs,p_aff_log,nrow = 1,labels=c('A)','B)','C)'), rel_widths = c(1.17,1.43,1.17))
p_pcr_0   = p_pcr
p_rel_0   = p_rel
p_abs_0   = p_abs
y_lumen_0 = y_lumen
aff_end_0 = aff_end

############################# ALLERGIC MUM
load('./SIM_ALLERGIC.RDATA')
source('./misc/PLOT_IGA_FUNCTION.R') #p_iga_function_merged
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
source('./misc/PLOT_ABUNDANCE_FIG5.R') #p_abs, p_rel, p_coating, p_micro
spacer_grob = grid::nullGrob()

p_abs = p_abs +scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_pcr = p_pcr + geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_pcr = p_pcr + geom_hline(yintercept = 0.50, color = "red", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_aff_log = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_coating = p_coating + scale_x_discrete(labels=c("C","BC","B","E"))
p_coating = p_coating + geom_hline(yintercept = 0.6, linetype = "dashed", color = "black", size = 1)
p_rel     = p_rel + theme(legend.position = "none")
row_1     = cowplot::plot_grid(p_rel,p_abs,p_aff_log,nrow = 1,labels=c('D)','E)','F)'), rel_widths = c(1.17,1.43,1.17))
p_pcr_1   = p_pcr
p_rel_1   = p_rel
p_abs_1   = p_abs
y_lumen_1 = y_lumen
aff_end_1 = aff_end

############################# IGA DEFICIENT MUM
load('./SIM_IGA_DEF.RDATA')
source('./misc/PLOT_IGA_FUNCTION.R') #p_iga_function_merged
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
source('./misc/PLOT_ABUNDANCE_FIG5.R') #p_abs, p_rel, p_coating, p_micro
spacer_grob = grid::nullGrob()

my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
p_abs     = p_abs +scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_pcr = p_pcr + geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_pcr = p_pcr + geom_hline(yintercept = 0.50, color = "red", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_aff_log = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_coating = p_coating + scale_x_discrete(labels=c("C","BC","B","E"))
p_coating = p_coating + geom_hline(yintercept = 0.6, linetype = "dashed", color = "black", size = 1)
p_rel     = p_rel + theme(legend.position = "none")
row_2     = cowplot::plot_grid(p_rel,p_abs,p_aff_log,nrow = 1,labels=c('G)','H)','I)'), rel_widths = c(1.17,1.43,1.17))
p_pcr_2   = p_pcr
p_rel_2   = p_rel
p_abs_2   = p_abs
y_lumen_2 = y_lumen
aff_end_2 = aff_end

############################# ONLY ECF
load('./SIM_IGA_ECF.RDATA')
source('./misc/PLOT_IGA_FUNCTION.R') #p_iga_function_merged
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
source('./misc/PLOT_ABUNDANCE_FIG5.R') #p_abs, p_rel, p_coating, p_micro
spacer_grob = grid::nullGrob()

my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
p_abs     = p_abs +scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_pcr = p_pcr + geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_pcr = p_pcr + geom_hline(yintercept = 0.50, color = "red", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_aff_log = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_coating = p_coating + scale_x_discrete(labels=c("C","BC","B","E"))
p_coating = p_coating + geom_hline(yintercept = 0.6, linetype = "dashed", color = "black", size = 1)
p_rel     = p_rel + theme(legend.position = "none")
row_3    = cowplot::plot_grid(p_rel,p_abs,p_aff_log,nrow = 1,labels=c('J)','K)','L)'), rel_widths = c(1.17,1.43,1.17))
p_pcr_3  = p_pcr
p_rel_3  = p_rel
p_abs_3  = p_abs
y_lumen_3 = y_lumen
aff_end_3 = aff_end

############################# ONLY ECF + PROBIOTICS
load('./SIM_IGA_ECF_PROBIOTICS.RDATA')
source('./misc/PLOT_IGA_FUNCTION.R') #p_iga_function_merged
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
source('./misc/PLOT_ABUNDANCE_FIG5.R') #p_abs, p_rel, p_coating, p_micro
spacer_grob = grid::nullGrob()

my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
p_abs     = p_abs +scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_pcr = p_pcr + geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_pcr = p_pcr + geom_hline(yintercept = 0.50, color = "red", linetype = "dashed", size = 1.5)  # Add a blue dashed horizontal line at y=0
p_aff_log = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_coating = p_coating + scale_x_discrete(labels=c("C","BC","B","E"))
p_coating = p_coating + geom_hline(yintercept = 0.6, linetype = "dashed", color = "black", size = 1)
p_rel     = p_rel + theme(legend.position = "none")
row_4    = cowplot::plot_grid(p_rel,p_abs,p_aff_log,nrow = 1,labels=c('M)','N)','O)'), rel_widths = c(1.17,1.43,1.17))
p_pcr_4  = p_pcr
p_rel_4  = p_rel
p_abs_4  = p_abs
y_lumen_4 = y_lumen
aff_end_4 = aff_end

label_plot_0 = ggdraw() + draw_label("Control", angle = 90, size = 14, hjust = 0.5)
label_plot_1 = ggdraw() + draw_label("Hyperreactive mSIgA in BM", angle = 90, size = 14, hjust = 0.5)
label_plot_2 = ggdraw() + draw_label("SIgA Deficient BM", angle = 90, size = 14, hjust = 0.5)
label_plot_3 = ggdraw() + draw_label("Only ECF", angle = 90, size = 14, hjust = 0.5)
label_plot_4 = ggdraw() + draw_label("Only ECF + Probiotics", angle = 90, size = 14, hjust = 0.5)

row_0 = cowplot::plot_grid(label_plot_0,row_0, nrow = 1, rel_widths = c(.15,10))
row_1 = cowplot::plot_grid(label_plot_1,row_1, nrow = 1, rel_widths = c(.15,10))
row_2 = cowplot::plot_grid(label_plot_2,row_2, nrow = 1, rel_widths = c(.15,10))
row_3 = cowplot::plot_grid(label_plot_3,row_3, nrow = 1, rel_widths = c(.15,10))
row_4 = cowplot::plot_grid(label_plot_4,row_4, nrow = 1, rel_widths = c(.15,10))

nested = cowplot::plot_grid(
  row_0, row_1, row_2, row_3,row_4,
  nrow = 5, 
  rel_heights = c(1, 1, 1, 1, 1), 
  vjust = 1.5# Adjust horizontal justification to create padding
)

graphics.off()
png(file =paste0("./FIGURE_S5.png"),     # The directory you want to save the file in
    width     = 14,
    height    = 12,
    units     = "in",
    res       = 300)
nested
dev.off()