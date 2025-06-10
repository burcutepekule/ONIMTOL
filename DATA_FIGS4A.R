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

add_pred          = 0 # Time to add for prediction (Total simulation time = 720 + add_pred days)
e_on_in           = 1 # Endogenous immune system on/off button (=1 endogenous system on, =0, endogenous system off) 
level_hmo_in      = 1 # Scales the HMOs levels, 0 to 1 (=1, HMOs level in the control case)
level_miga_in     = 1 # Scales the mSIgA levels, 0 to 1 (=0 mimics SIgA deficiency in the breastmilk)
mcell_th          = 0.5 # M cells open when mSIgA drops below (100*mcell_th)% of its maximum value. 
EBF_duration      = 154 # Exclusive Breastfeeding duration (days)
MF_duration       = 0 # Mixed Feeding duration (days) -> THIS IS ZERO FOR THIS PARTICULAR FIGURE
mIgA_in           = NA # mSIgA affinities. NA is the baseline scenario, where mIgA_in vector is c(4.18,0.18,0.13,0.12) for E, B, BC, and C

source('./misc/LOAD_DATA_LIST.R')
source('./misc/DATA_LIST_2_PARMS.R')
source('./misc/INTEGRATE_DDE.R')
gc()

plot1 = 0 # to plot the unseen data gray rectangle
############# PLOTTING
y_out_df = y_out_df[1:240,] # filter until day 240
source('./misc/PLOT_CALORIES.R') # p_cals
p_cals = p_cals+theme(
  plot.title = element_text(size = 12),
  legend.position = "bottom",
  legend.title = element_blank(),
  legend.box = "horizontal",
  legend.justification = 'center',
  legend.background = element_rect(colour = "gray", fill = "white")
) +
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE),
    linetype = guide_legend(nrow = 3, byrow = TRUE)
  )

### IF YOU WANNA SAVE THE PLOT
graphics.off()
png(file =paste0("./FIGURE_S4A.png"),    # The directory you want to save the file in
    width     = 5,
    height    = 5,
    units     = "in",
    res       = 300)
p_cals
dev.off()
