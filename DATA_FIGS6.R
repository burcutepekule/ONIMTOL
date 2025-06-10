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
my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')

########### Control
load('IMAGE_DATA_S5A.RDATA')
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
p_aff_log_0   = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))

########### SIgA Deficient
load('IMAGE_DATA_S5B.RDATA')
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
p_aff_log_1   = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))

########### SIgA Deficient + GR- adjustment
load('IMAGE_DATA_S5C.RDATA')
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
p_aff_log_2 = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))

p_aff_log_0 = p_aff_log_0 + theme(legend.title = element_blank(),legend.position = c(1, 0.134))
p_aff_log_1 = p_aff_log_1 + theme(legend.title = element_blank(),legend.position = c(1, 0.134))
p_aff_log_2 = p_aff_log_2 + theme(legend.title = element_blank(),legend.position = c(1, 0.134))

row_final = cowplot::plot_grid(p_aff_log_0,p_aff_log_1,p_aff_log_2,nrow = 1,labels=c('A)','B)','C)'), rel_widths = c(1,1,1), label_fontface = "bold")

graphics.off()
png(file =paste0("./FIGURE_S6.png"),    # The directory you want to save the file in
    width     = 13,
    height    = 4,
    units     = "in",
    res       = 300)
print(row_final)
dev.off()

