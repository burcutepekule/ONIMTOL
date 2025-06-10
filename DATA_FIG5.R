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

################################################################################

y_lumen_0$Feeding = 'Control'
y_lumen_1$Feeding = 'Hyperreactive SIgA BM'
y_lumen_2$Feeding = 'SIgA Deficient BM'
y_lumen_3$Feeding = 'Only ECF'
y_lumen_4$Feeding = 'Only ECF + Probiotics'
y_lumen_all = as.data.frame(rbind(y_lumen_0,y_lumen_1,y_lumen_2,y_lumen_3,y_lumen_4))

color_vec = c(
  "#1E7A1E",  
  "magenta", 
  "#FFB400",
  "#66CDAA", 
  "#800080"  
)

# Step 1: Separate data for Enterobacteriaceae and others
entero_data <- y_lumen_all %>%  
  filter(taxa == "Enterobacteriaceae")

colnames(entero_data) = c('days','taxa','killed_lumen_e','uncoated_lumen_e',
                          'coated_lumen_e','total_feces_e','igaplus_feces_e',
                          'coated_feces_e','killed_feces_e','uncoated_feces_e','Feeding')

# Group and summarize for all other taxa
other_taxa_sum <- y_lumen_all %>%
  dplyr::filter(taxa != "Enterobacteriaceae") %>%
  dplyr::group_by(days, Feeding) %>%
  dplyr::summarize(killed_lumen_ne   = sum(median_k),
                   uncoated_lumen_ne = sum(median_uc),
                   coated_lumen_ne   = sum(median_c),
                   total_feces_ne    = sum(median_tot),
                   igaplus_feces_ne  = sum(median_ck_tot),
                   coated_feces_ne   = sum(median_c_tot),
                   killed_feces_ne   = sum(median_k_tot),
                   uncoated_feces_ne = sum(median_uc_tot),.groups = "drop")

merged_data = merge(entero_data,other_taxa_sum, by=c('days','Feeding'))
merged_data = merged_data %>% dplyr::filter(days<=30 & days>=1)

total_abundance_scale    = readRDS('total_abundance_scale.rds')
total_abundance_scale_df = data.frame(days=seq(0,length(total_abundance_scale)-1,1), mult=total_abundance_scale)
total_abundance_scale_df = total_abundance_scale_df %>% dplyr::filter(days<=30 & days>=1)

merged_data = merge(merged_data,total_abundance_scale_df,by=c('days'))
merged_data = merged_data %>% dplyr::rowwise() %>% dplyr::mutate(ratio_1 = mult*uncoated_feces_e/(total_feces_e+total_feces_ne))
merged_data = merged_data %>% dplyr::rowwise() %>% dplyr::mutate(ratio_2 = uncoated_lumen_e/(uncoated_lumen_ne+coated_lumen_ne))
merged_data = merged_data %>% dplyr::rowwise() %>% dplyr::mutate(ratio_c = mult*igaplus_feces_e/(total_feces_e+total_feces_ne))
merged_data = merged_data %>% dplyr::rowwise() %>% dplyr::mutate(ratio_t = mult*(total_feces_e)/(total_feces_e+total_feces_ne))

p_rat = ggplot(merged_data, aes(x = days, y = ratio_1, color = Feeding, linetype = Feeding)) +
  geom_line(size=1.0) +
  scale_color_manual(values = color_vec,
                     breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics"),
                     guide = guide_legend(override.aes = list(color = "white"))) +
  scale_linetype_manual(values = c("solid", "11", "dashed", "dotdash", "B5"), # Adjust these linetypes as needed
                        breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics")) +
  theme_minimal() + theme(axis.title.x = element_text(size =12), axis.title.y = element_text(size =12), 
                          plot.title = element_text(size = 12), 
                          legend.text = element_text(size = 12, color = "transparent"), 
                          legend.title = element_text(color = "transparent"),
                          legend.background = element_rect(fill = NA, colour = NA), # Make legend background transparent
                          legend.key = element_rect(fill = NA, colour = NA))+
  labs(x = "DOL", y = "Relative abundance (%)", title = "SIgA- Enterobacteriaceae, fecal samples")



merged_data_10 = merged_data %>% dplyr::filter(days==10)
labels = paste0(round(merged_data_10$ratio_1,1), '%')
pts = merged_data_10$ratio_1


p_rat_t = ggplot(merged_data, aes(x = days, y = ratio_t, color = Feeding, linetype = Feeding)) +
  geom_line(size=1.0) +
  scale_color_manual(values = color_vec,
                     breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics"))+
  scale_linetype_manual(values = c("solid", "11", "dashed", "dotdash", "B5"), # Adjust these linetypes as needed
                        breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics")) +
  theme_minimal() + theme(axis.title.x = element_text(size =12), axis.title.y = element_text(size =12), plot.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  labs(x = "DOL", y = "Relative abundance (%)", title = "Enterobacteriaceae, fecal samples")

labels_t = paste0(round(merged_data_10$ratio_t,1), '%')
pts_t = merged_data_10$ratio_t


p_rat_c = ggplot(merged_data, aes(x = days, y = ratio_c, color = Feeding, linetype = Feeding)) +
  geom_line(size=1.0) +
  scale_color_manual(values = color_vec,
                     breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics"))+
  scale_linetype_manual(values = c("solid", "11", "dashed", "dotdash", "B5"), # Adjust these linetypes as needed
                        breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics")) +
  theme_minimal() + theme(axis.title.x = element_text(size =12), axis.title.y = element_text(size =12), plot.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  labs(x = "DOL", y = "Relative abundance (%)", title = "SIgA+ Enterobacteriaceae, fecal samples")

labels_c = paste0(round(merged_data_10$ratio_c,1), '%')
pts_c = merged_data_10$ratio_c

p_abs = ggplot(merged_data, aes(x = days, y = uncoated_lumen_e, color = Feeding, linetype=Feeding)) +
  geom_line(size=1.0) +
  scale_color_manual(values = color_vec,
                     breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics"))+
  scale_linetype_manual(values = c("solid", "11", "dashed", "dotdash", "B5"), # Adjust these linetypes as needed
                        breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics")) +
  theme_minimal() + theme(axis.title.x = element_text(size =12), axis.title.y = element_text(size =12), plot.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  labs(x = "DOL", y = "(Cells/g LC) x 1e-11", title = "SIgA- Enterobacteriaceae, gut lumen")

labels_2 = round(merged_data_10$uncoated_lumen_e,3)
pts_2 = merged_data_10$uncoated_lumen_e

p_lum = ggplot(merged_data, aes(x = days, y = ratio_2, color = Feeding, linetype=Feeding)) +
  geom_line(size=1.0) +
  scale_color_manual(values = color_vec,
                     breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics"))+
  scale_linetype_manual(values = c("solid", "11", "dashed", "dotdash", "B5"), # Adjust these linetypes as needed
                        breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics")) +
  theme_minimal() + theme(axis.title.x = element_text(size =12), axis.title.y = element_text(size =12), plot.title = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  labs(x = "DOL", y = "Ratio", title = "Pathogenic/Symbiotic commensals, lumen")+
  scale_y_log10() 
# annotate(geom = "rect", xmin = 0, xmax = 30, ymin = 0, ymax = 1, alpha = 0.1, fill = "#176B87") # Add this line
p_lum

labels_2 = round(merged_data_10$uncoated_lumen_e,3)
pts_2    = merged_data_10$uncoated_lumen_e
p_abs    = p_abs + theme(legend.position = "none")
p_rat_c  = p_rat_c + theme(legend.position = "none")
p_lum    = p_lum + theme(legend.position = "none")

##### REPEAT p_rat with legends #########

p_rat = ggplot(merged_data, aes(x = days, y = ratio_1, color = Feeding, linetype = Feeding)) +
  geom_line(size=1.0) +
  scale_color_manual(values = color_vec,
                     breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics")) +
  scale_linetype_manual(values = c("solid", "11", "dashed", "dotdash", "B5"), # Adjust these linetypes as needed
                        breaks = c("Control", "Hyperreactive SIgA BM", "SIgA Deficient BM", "Only ECF", "Only ECF + Probiotics")) +
  theme_minimal() + theme(axis.title.x = element_text(size =12), axis.title.y = element_text(size =12), 
                          plot.title = element_text(size = 12), legend.text = element_text(size=12))+
  labs(x = "DOL", y = "Relative abundance (%)", title = "SIgA- Enterobacteriaceae, fecal samples")


row_final = cowplot::plot_grid(p_lum,p_rat_t,p_rat_c,p_rat,nrow = 2,labels=c('A','B','C','D'), rel_widths = c(1,1.5,1,1.5), label_fontface = "bold")
row_final = cowplot::plot_grid(p_lum,p_rat_c,p_rat,nrow = 1,labels=c('A','B','C'), rel_widths = c(1.1,1,1.6), label_fontface = "bold")

graphics.off()
png(file =paste0("./FIGURE_5.png"),    # The directory you want to save the file in
    width     = 12.9,
    height    = 3.5,
    units     = "in",
    res       = 300)
print(row_final)
dev.off()

