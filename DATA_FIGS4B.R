#################################
rm(list=ls())
library(tidyverse)
library(cowplot)  # this is a plotting library but I'm unsure if you want to retain it
library(readxl)
library(icesTAF)
library(plyr)
library(dplyr)
library(writexl)
library(loo)
library(mvtnorm)
library(MASS)
library(deSolve)
library(relaimpo)
library(rwa)
library(gridExtra)
library(ggtext)
library(WRS2)
library(boot)
library(effsize)
library(rstatix)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

pdc_data_use = readRDS('pdc_data_use.rds') # simulation data

## TURN TO NUMERIC
pdc_data_use$mf_length = as.numeric(as.character(pdc_data_use$mf_length))
pdc_data_use$bm_cutoff = as.numeric(as.character(pdc_data_use$bm_cutoff))

pdc_data_use_long = pdc_data_use[c('bm_cutoff','mf_length','aff_e','aff_b','aff_bc','aff_c')]
pdc_data_use_long = pdc_data_use_long %>% filter(bm_cutoff %in% seq(0,180,15))
unique_bm_cutoff  = unique(pdc_data_use_long$bm_cutoff)
scaledSize        = seq(1,length(unique_bm_cutoff),1)
colnames(pdc_data_use_long)[3:6] = c('Enterobacteriaceae','Bifidobacteriaceae','Bacteroidaceae','Clostridiales') 
pdc_data_use_long = pdc_data_use_long %>% pivot_longer(!c(bm_cutoff,mf_length),values_to = 'aff', names_to='taxa')

pdc_data_use_long$mf_length = as.numeric(pdc_data_use_long$mf_length)
pdc_data_use_long$bm_cutoff = as.numeric(pdc_data_use_long$bm_cutoff)

pdc_data_use_long = pdc_data_use_long %>% dplyr::rowwise() %>% dplyr::mutate(total_duration = bm_cutoff+mf_length)

# Define colors for consistency
colors <- c(
  'Enterobacteriaceae' = "#f20026",
  'Bifidobacteriaceae' = "#074fa8",
  'Bacteroidaceae' = "#2bc5d8",
  'Clostridiales' = "#ff9aa1"
)

graphics.off()
p_aff_idx = ggplot(pdc_data_use_long) +
  geom_point(aes(x = total_duration, y = aff, fill=taxa, size = bm_cutoff), color = "black", shape = 21, alpha = 0.7) +
  scale_fill_manual(values = colors,
                    labels = c("Bifidobacteriaceae", "Bacteroidaceae", "Clostridiales","Enterobacteriaceae"),
                    name = "Taxa") + # Update this line) +
  scale_size_continuous(range = c(1, max(scaledSize)), breaks = unique_bm_cutoff, name = "EBF Duration (days)") +
  labs(x = "Total feeding duration (EBF+MF, days)", y = "log (eSIgA Affinity)", title = "") +
  theme_minimal()


p_aff_idx = p_aff_idx +  geom_hline(yintercept =  1, colour="gray10", linetype = "longdash")

pdc_data_keep_allgood     = pdc_data_use %>% filter(aff_b <1 & aff_bc <1 & aff_c <1) #50% coating
pdc_data_keep_allgood_2   = pdc_data_use %>% filter(aff_b <0.33 & aff_bc <0.33 & aff_c <0.33) #75% coating
pdc_data_keep_allgood_3   = pdc_data_use %>% filter(aff_b <0.66 & aff_bc <0.66 & aff_c <0.66) #60% coating
pdc_data_use_mf_0         = pdc_data_use %>% filter(mf_length==0)

pdc_data_keep_allbad    = pdc_data_use %>% filter(aff_b >1 & aff_bc >1 & aff_c >1)
pdc_data_keep_bifigood  = pdc_data_use %>% filter(aff_b <1 & aff_bc >1 & aff_c >1)

pdc_data_keep_mf0       = pdc_data_use %>% filter(mf_length==0)
pdc_data_keep_bm0       = pdc_data_use %>% filter(bm_cutoff==0)

vline_1 = min(pdc_data_keep_allgood$total_duration)
vline_2 = min(pdc_data_keep_allgood_2$total_duration)
vline_3 = min(pdc_data_keep_allgood_3$total_duration)
vline_4 = max(pdc_data_use_mf_0$total_duration)


pdc_data_keep_allgood_1 = pdc_data_keep_allgood %>% filter(total_duration==vline_1)

aff_eff =1/(1+c(mean(pdc_data_keep_allgood_1$aff_b), mean(pdc_data_keep_allgood_1$aff_bc), mean(pdc_data_keep_allgood_1$aff_c)))
aff_eff = round(aff_eff,4)

p_aff_idx = p_aff_idx + scale_y_continuous(trans = 'log',  # Define log-spaced breaks
                                           labels = function(x) round(x, 1),
                                           breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 20))  # Round the labels to show on the axis



linename = "No Mixed Feeding \nMF=0 days"
p_aff_idx = p_aff_idx + 
  geom_line(data = subset(pdc_data_use, mf_length == 0), aes(x = total_duration, y = aff_e, color = "No Mixed Feeding \nMF=0 days"), alpha = 1, size=1.25) +
  geom_line(data = subset(pdc_data_use, mf_length == 0), aes(x = total_duration, y = aff_b, color = "No Mixed Feeding \nMF=0 days"), alpha = 1, size=1.25) +
  geom_line(data = subset(pdc_data_use, mf_length == 0), aes(x = total_duration, y = aff_bc, color = "No Mixed Feeding \nMF=0 days"), alpha = 1, size=1.25) +
  geom_line(data = subset(pdc_data_use, mf_length == 0), aes(x = total_duration, y = aff_c, color = "No Mixed Feeding \nMF=0 days"), alpha = 1, size=1.25) +
  scale_color_manual(name = "",values = c("No Mixed Feeding \nMF=0 days" = "gray10"))

# SAVE PNG
graphics.off()
png(file =paste0("./FIGURE_S4B.png"),   # The directory you want to save the file in
    width     = 8,
    height    = 8,
    units     = "in",
    res       = 300)
p_aff_idx
dev.off()