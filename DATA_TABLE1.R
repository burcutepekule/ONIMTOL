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
load('./SIM_BASELINE.RDATA') # y_out_df is the dataframe of the simulation output!
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
mIgA_in_0 = mIgA_in # mother's IgA affinities

############################# ALLERGIC MUM
load('./SIM_ALLERGIC.RDATA') # y_out_df is the dataframe of the simulation output!
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
mIgA_in_1 = mIgA_in # mother's IgA affinities

############################# IGA DEFICIENT MUM
load('./SIM_IGA_DEF.RDATA') # y_out_df is the dataframe of the simulation output!
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
mIgA_in_2 = rep(NA,4) # IgA deficient mum this case

############################# ONLY ECF
load('./SIM_IGA_ECF.RDATA') # y_out_df is the dataframe of the simulation output!
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
mIgA_in_3 = rep(NA,4) # no breastfeeding for this case

############################# ONLY ECF + PROBIOTICS
load('./SIM_IGA_ECF_PROBIOTICS.RDATA') # y_out_df is the dataframe of the simulation output!
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
mIgA_in_4 = rep(NA,4) # no breastfeeding for this case

################################################################################
aff_end_0$Feeding = 'Control'
aff_end_1$Feeding = 'Hyperreactive SIgA BM'
aff_end_2$Feeding = 'SIgA Deficient BM'
aff_end_3$Feeding = 'Only ECF'
aff_end_4$Feeding = 'Only ECF + Probiotics'

aff_end_all_offspring = as.data.frame(rbind(aff_end_0,aff_end_1,aff_end_2,aff_end_3,aff_end_4))
aff_end_all_offspring = aff_end_all_offspring %>% dplyr::mutate(neutralizing_rate = round(1-1/(1+value),2)) 
aff_end_all_offspring$Host  = 'Offspring' #This is the table for the offspring values

aff_end_all_mother = aff_end_all_offspring # repeat first, and then fill in 
aff_end_all_mother$Host  = 'Mother' #This is the table for the offspring values

# here, pay attention to the order. mIgA values are ordered as E, B, BC, C. 
# In the table, it is E, BC, B, C. So re-order before filling in.

mIgA_in_0 = mIgA_in_0[c(1,3,2,4)]
mIgA_in_1 = mIgA_in_1[c(1,3,2,4)]
mIgA_in_2 = mIgA_in_2[c(1,3,2,4)]
mIgA_in_3 = mIgA_in_3[c(1,3,2,4)]
mIgA_in_4 = mIgA_in_4[c(1,3,2,4)]

aff_end_all_mother$value = c(mIgA_in_0, mIgA_in_1, mIgA_in_2, mIgA_in_3, mIgA_in_4) #insert
aff_end_all_mother = aff_end_all_mother %>% dplyr::mutate(neutralizing_rate = round(1-1/(1+value),2)) #convert from affinity to neutralizing rate

TABLE_1 = rbind(aff_end_all_mother, aff_end_all_offspring)
TABLE_1 = TABLE_1[c('taxa','Feeding','neutralizing_rate', 'Host')]
print(TABLE_1)
