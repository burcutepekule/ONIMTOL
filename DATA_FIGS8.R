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

### LOAD DATA
keep_aff_sims_1 = readRDS('./DATA_BASELINE_AFFINITY.rds')
keep_aff_sims_2 = readRDS('./DATA_2FOLD_AFFINITY.rds')
keep_aff_sims_3 = readRDS('./DATA_5FOLD_AFFINITY.rds')

keep_aff_sims_1$alpha_idx = 'Baseline'
keep_aff_sims_2$alpha_idx = '2-fold increase'
keep_aff_sims_3$alpha_idx = '5-fold increase'

keep_aff_sims = rbind(keep_aff_sims_1,keep_aff_sims_2,keep_aff_sims_3)
colnames(keep_aff_sims)[1] = 'taxon' 

taxon_colors <- c(
  "Enterobacteriaceae" = "#f20026",  # Red
  "Bifidobacteriaceae" = "#074fa8",  # Blue
  "Bacteroidaceae" = "#2bc5d8",      # Green
  "Clostridiales" = "#ff9aa1"        # Purple
)

keep_aff_sims$alpha_idx_factor <- factor(keep_aff_sims$alpha_idx, 
                                         levels = c('Baseline', '2-fold increase', '5-fold increase'))  # Change order as needed


keep_aff_sims_ne = keep_aff_sims %>% filter(taxon %in% c("Bifidobacteriaceae","Bacteroidaceae","Clostridiales"))

p=ggplot(keep_aff_sims_ne, aes(x = (mIgA), y = (eIgA), color = taxon, shape = factor(alpha_idx_factor))) +
  geom_point(size = 3, alpha = 1) +
  scale_shape_manual(values = c(16, 5, 17), name = "Invasiveness of Symbiotic\nCommensals (B, BC, C)") +
  scale_color_manual(values = taxon_colors) +
  labs(
    title = "mSIgA vs eSIgA affinities by taxonomic group and invasiveness",
    x = "mSIgA Affinity",
    y = "eSIgA Affinity",
    color = "Taxonomic Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

graphics.off()
png(file =paste0("FIGURE_S8.png"),    # The directory you want to save the file in
    width     = 8,
    height    = 5,
    units     = "in",
    res       = 300)
print(p)
dev.off()




