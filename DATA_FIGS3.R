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

######## DYSBIOSIS VIA LACK OF HMOS ################################ 
level_hmo_in      = 0.25 # Scales the HMOs levels, 0 to 1 (=1, HMOs level in the control case)
######## DYSBIOSIS VIA LACK OF HMOS ################################ 

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

plot1 = 0 # to plot the unseen data gray rectangle
############# PLOTTING
source('./misc/PLOT_CALORIES.R') # p_cals
source('./misc/PLOT_IGA_FUNCTION.R') #p_iga_function_merged
source('./misc/PLOT_AFF_EIGA.R') # p_aff_log, p_eiga
source('./misc/PLOT_ABUNDANCE.R') #p_abs, p_rel, p_coating, p_micro

my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
p_abs     = p_abs +scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_aff_log = p_aff_log + scale_color_manual(values = my_colors,labels = c("BC", "B", "C", "E"))
p_aff_log = p_aff_log + geom_vline(xintercept = 129, linetype = "dotted", color = "black")
p_aff_log = p_aff_log + geom_vline(xintercept = 720, linetype = "dotted", color = "black")
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

### MORE MINIMAL

# Generate ellipse coordinates
theta <- seq(0, 2*pi, length.out = 100)
x_center <- 80      # X-coordinate of ellipse center
y_center <- 0.3     # Y-coordinate of ellipse center
a <- 200             # Semi-major axis length (X direction)
b <- .3             # Semi-minor axis length (Y direction)

# Create data frame with ellipse coordinates
ellipse_data <- data.frame(
  x = x_center + a * cos(theta),
  y = y_center + b * sin(theta)
)

p_rel = p_rel + geom_polygon(
  inherit.aes = FALSE,  # This is key - don't inherit aesthetics from the main plot
  data = ellipse_data,
  mapping = aes(x = x, y = y),
  fill = "magenta",
  alpha = 0.0,
  color = "magenta",
  size = 1  # Increase this value to make the contour thicker (default is typically 0.5)
)

# Define the rectangle parameters
x_min <- 100         # Left edge of rectangle
x_max <- 750         # Right edge of rectangle
y_min <- 0.3        # Bottom edge of rectangle
y_max <- 4        # Top edge of rectangle

# Create rectangle points (vertices in clockwise order)
rectangle_points <- data.frame(
  x = c(x_min, x_max, x_max, x_min, x_min),  # Last point connects back to first
  y = c(y_min, y_min, y_max, y_max, y_min)   # Last point connects back to first
)

p_aff_log = p_aff_log + geom_polygon(
  data = rectangle_points,
  mapping = aes(x = x, y = y),
  inherit.aes = FALSE,
  fill = "magenta",
  alpha = 0.0,       # Transparent filling
  color = "magenta", # Magenta contour
  size = 1           # Thick contour
)

nested = plot_grid(p_rel, p_micro, p_aff_log, ncol = 3, rel_widths = c(1.6,1,1),labels=c('A)','B)','C)'))

# ### IF YOU WANNA SAVE THE PLOT
graphics.off()
png(file = "FIGURE_S3.png",    # The directory you want to save the file in
    width     = 12,
    height    = 3,
    units     = "in",
    res       = 300)
nested
dev.off()

