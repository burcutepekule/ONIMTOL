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
source('./misc/MODEL_BASE_FUNCTIONS.R')
source('./misc/LOAD_DATA.R')

add_pred          = 180 # Time to add for prediction (Total simulation time = 720 + add_pred days)
e_on_in           = 1 # Endogenous immune system on/off button (=1 endogenous system on, =0, endogenous system off) 
level_hmo_in      = 1 # Scales the HMOs levels, 0 to 1 (=1, HMOs level in the control case)
level_miga_in     = 1 # Scales the mSIgA levels, 0 to 1 (=0 mimics SIgA deficiency in the breastmilk)
mcell_th          = 0.5 # M cells open when mSIgA drops below (100*mcell_th)% of its maximum value. 
EBF_duration      = 154 # Exclusive Breastfeeding duration (days)
MF_duration       = 308 # Mixed Feeding duration (days)
mIgA_in           = NA # mSIgA affinities. NA is the baseline scenario, where mIgA_in vector is c(4.18,0.18,0.13,0.12) for E, B, BC, and C
# mIgA_in           = rep(4.18,4) # allergic/IBD case, hyperreactive against symbiotic commensals (using the affinity value for E for all)

source('./misc/LOAD_DATA_LIST.R')
source('./misc/DATA_LIST_2_PARMS.R')
source('./misc/INTEGRATE_DDE.R')
gc()

############# PLOTTING
source('./misc/PLOT_CALORIES.R') # p_cals
source('./misc/PLOT_IGA_FUNCTION.R') #p_iga_function_merged
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
source('./misc/PLOT_ABUNDANCE.R') #p_abs, p_rel, p_coating, p_micro

my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
p_abs     = p_abs +scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_aff_log = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_eiga    = p_eiga + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_rel     = p_rel + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_coating = p_coating + scale_x_discrete(labels=c("C","BC","B","E"))
p_iga     = p_iga + scale_x_discrete(labels=c("C","BC","B","E"))

spacer_grob = grid::nullGrob()
row1   = plot_grid(p_cals,p_iga_function_merged,align='h',rel_widths = c(1,1),labels=c('A)','B)'))
row2   = plot_grid(p_micro,p_aff_log,p_eiga,p_coating,p_iga, align='h',ncol = 5, rel_widths = c(0.4,0.4,0.4,0.35,0.3),labels=c('C)','D)','E)','F)','G)'))
row3   = plot_grid(p_rel,p_abs,align='h', rel_widths = c(1.1,1),labels=c('H)','I)'))
nested = plot_grid(row1, spacer_grob, row2, spacer_grob, row3, ncol = 1, rel_heights = c(1, 0.05, 1, 0.05, 1))

z_c  = z %>% filter(Measurement=='IgA_coated_fraction_c')
z_k  = z %>% filter(Measurement=='IgA_coated_fraction_k')
z_i  = z %>% filter(Measurement=='iga_index')
z_ck = z %>% filter(Measurement %in% c('IgA_coated_fraction_c','IgA_coated_fraction_k') & Value>1e-10)
print(paste0('Coating ratio : ',round(100*mean(z_ck$Value),2),' %'))

# ### IF YOU WANNA SAVE THE PLOT
graphics.off()
png(file =paste0("FIGURE_1_",EBF_duration,'_',MF_duration,".png"),    # The directory you want to save the file in
    width     = 13,
    height    = 10,
    units     = "in",
    res       = 300)
nested
dev.off()

p_aff_log
