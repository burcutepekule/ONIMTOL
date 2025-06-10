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
library(randomForest)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load('./IMAGE_DATA4.RDATA') # Combined files of simulations over varying EBF and MF durations, takes a few seconds to load

pdc_data_use = readRDS('pdc_data_use.rds') # simulation data
######## RANDOM FORESTS RELATIVE IMPORTANCE
##RF IS A STOCHASTIC METHOD, THEREFORE SET THE SEED FOR REPRODUCILBILITY
set.seed(42)

## TURN TO NUMERIC
pdc_data_use$mf_length = as.numeric(as.character(pdc_data_use$mf_length))
pdc_data_use$bm_cutoff = as.numeric(as.character(pdc_data_use$bm_cutoff))

pdc_data_use    = pdc_data_use %>% rowwise() %>% mutate(y0_b = y0_b_c+y0_b_uc)
pdc_data_use_2  = pdc_data_use %>% filter(bm_cutoff>0) # for this one you need EBF to exist

rf_model_1_e    = randomForest(aff_e ~ bm_cutoff + mf_length, data = pdc_data_use)
rf_model_1_b    = randomForest(aff_b ~ bm_cutoff + mf_length, data = pdc_data_use)
rf_model_1_bc   = randomForest(aff_bc ~ bm_cutoff + mf_length, data = pdc_data_use)
rf_model_1_c    = randomForest(aff_c ~ bm_cutoff + mf_length, data = pdc_data_use)

rf_model_2_bc   = randomForest(aff_bc ~ inf_t0 + y0_e_uc + y0_b, data = pdc_data_use_2)
rf_model_2_c    = randomForest(aff_c ~ inf_t0 + y0_e_uc + y0_b, data = pdc_data_use_2)

pdc_data_numeric_cor   = pdc_data_use_2[c('aff_bc','inf_t0','y0_e_uc','y0_b')]
cor_vec                = sign(cor(pdc_data_numeric_cor)[1,2:4])
importance_scores      = importance(rf_model_2_bc)[, "IncNodePurity"]
total_importance       = sum(importance_scores)
importance_percentages = (importance_scores / total_importance) * 100
importance_df_bc_2 = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = importance_percentages,
  taxa ='Bacteroidaceae',
  Direction = cor_vec[names(importance_percentages)]
)

pdc_data_numeric_cor   = pdc_data_use_2[c('aff_c','inf_t0','y0_e_uc','y0_b')]
cor_vec                = sign(cor(pdc_data_numeric_cor)[1,2:4])
importance_scores      = importance(rf_model_2_c)[, "IncNodePurity"]
total_importance       = sum(importance_scores)
importance_percentages = (importance_scores / total_importance) * 100
importance_df_c_2 = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = importance_percentages,
  taxa ='Clostridiales',
  Direction = cor_vec[names(importance_percentages)]
)

importance_df_cbc_2 = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = 0.5*(importance_df_c_2$AbsImportance+importance_df_bc_2$AbsImportance),
  taxa ='Bacteroidaceae \nand Clostridiales',
  Direction = cor_vec[names(importance_percentages)] # SAME ANYWAY
)

pdc_data_numeric_cor   = pdc_data_use[c('aff_e','bm_cutoff','mf_length')]
cor_vec                = sign(cor(pdc_data_numeric_cor)[1,2:3])
importance_scores      = importance(rf_model_1_e)[, "IncNodePurity"]
total_importance       = sum(importance_scores)
importance_percentages = (importance_scores / total_importance) * 100
importance_df_e = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = importance_percentages,
  taxa ='Enterobacteriaceae',
  Direction = cor_vec[names(importance_percentages)] 
)

pdc_data_numeric_cor   = pdc_data_use[c('aff_b','bm_cutoff','mf_length')]
cor_vec                = sign(cor(pdc_data_numeric_cor)[1,2:3])
importance_scores      = importance(rf_model_1_b)[, "IncNodePurity"]
total_importance       = sum(importance_scores)
importance_percentages = (importance_scores / total_importance) * 100
importance_df_b = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = importance_percentages,
  taxa ='Bifidobacteriaceae',
  Direction = cor_vec[names(importance_percentages)] 
)

pdc_data_numeric_cor   = pdc_data_use[c('aff_bc','bm_cutoff','mf_length')]
cor_vec                = sign(cor(pdc_data_numeric_cor)[1,2:3])
importance_scores      = importance(rf_model_1_bc)[, "IncNodePurity"]
total_importance       = sum(importance_scores)
importance_percentages = (importance_scores / total_importance) * 100
importance_df_bc = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = importance_percentages,
  taxa ='Bacteroidaceae',
  Direction = cor_vec[names(importance_percentages)] 
)

pdc_data_numeric_cor   = pdc_data_use[c('aff_c','bm_cutoff','mf_length')]
cor_vec                = sign(cor(pdc_data_numeric_cor)[1,2:3])
importance_scores      = importance(rf_model_1_c)[, "IncNodePurity"]
total_importance       = sum(importance_scores)
importance_percentages = (importance_scores / total_importance) * 100
importance_df_c = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = importance_percentages,
  taxa ='Clostridiales',
  Direction = cor_vec[names(importance_percentages)] 
)

importance_df_cbc = data.frame(
  Variables = names(importance_percentages),
  AbsImportance = 0.5*(importance_df_c$AbsImportance+importance_df_bc$AbsImportance),
  taxa ='Bacteroidaceae \nand Clostridiales',
  Direction = cor_vec[names(importance_percentages)] # SAME ANYWAY
)

importance_df = rbind(importance_df_e,importance_df_b,importance_df_bc, importance_df_c)
importance_df = importance_df %>% dplyr::rowwise() %>% dplyr::mutate(Variables = ifelse(Variables=='bm_cutoff','Exclusive Breastfeeding Duration','Mixed Feeding Duration'))
importance_df = importance_df %>% add_row(Variables='Invisible',AbsImportance=0,taxa='',Direction=-1)
importance_df$alpha = ifelse(importance_df$Variables == "Invisible", 0, 1)
importance_df$taxa  = factor(importance_df$taxa, levels = c("Clostridiales","Bacteroidaceae","Bifidobacteriaceae","Enterobacteriaceae" ,""))
importance_df = importance_df%>%dplyr::rowwise()%>%dplyr::mutate(Importance = Direction*AbsImportance)

importance_df_cbc_2 = rbind(importance_df_c_2, importance_df_bc_2) # EDIT HERE TO PUT BOTH C AND BC

importance_df_cbc_2 = importance_df_cbc_2 %>% dplyr::rowwise() %>% dplyr::mutate(Variables = ifelse(Variables=='inf_t0','Microenvironmental Stimulation',
                                                                                                    ifelse(Variables=='y0_e_uc','Enterobacteriaceae Abundance',
                                                                                                           'Bifidobacteriaceae Abundance')))

importance_df_cbc_2 = importance_df_cbc_2 %>% dplyr::rowwise() %>% dplyr::mutate(alpha=ifelse(Variables %in% c("Invisible1","Invisible2","Invisible3"), 0, 1))
importance_df_cbc_2 = importance_df_cbc_2%>%dplyr::rowwise()%>%dplyr::mutate(Importance = Direction*AbsImportance)


importance_df_cbc_2$Variables = factor(importance_df_cbc_2$Variables, levels = c("Microenvironmental Stimulation",
                                                                                 "Bifidobacteriaceae Abundance",
                                                                                 "Enterobacteriaceae Abundance",
                                                                                 "Invisible1"))


####### plot

p_imp_1_pos = ggplot(importance_df, aes(x = taxa, y = abs(Importance), fill = Variables)) +
  geom_col(position = position_dodge()) +
  coord_flip() +
  labs(y = "Relative Importance (%)", x = "", title = "") +
  theme_minimal() +
  scale_fill_manual(values = c('Exclusive Breastfeeding Duration' = '#FF006E',
                               'Mixed Feeding Duration' = '#12db99',
                               'Invisible' = '#FFFFFF'),  # Assuming 'Invisible' is still needed
                    breaks = c('Exclusive Breastfeeding Duration', 'Mixed Feeding Duration')) +  # Specify visible legend entries
  scale_alpha_identity() +  # Use the alpha values as is
  # geom_hline(yintercept = 50, linetype = "dashed")+
  theme(legend.background = element_rect(fill = "white", colour = "black"),legend.position = c(0.25, 0.92)) # Set legend background to white

p_imp_2_pos = ggplot(importance_df_cbc_2, aes(x = taxa, y = abs(Importance), fill = Variables)) +
  geom_col(position = position_dodge()) +
  coord_flip() +
  labs(y = "Relative Importance (%)", x = "", title = "", fill = "Variables (at start of weaning)", size=1) +  # Add your legend title here
  theme_minimal() +
  scale_fill_manual(values = c("Microenvironmental Stimulation" = '#3BF4FB',
                               "Enterobacteriaceae Abundance"='#f20026',
                               "Bifidobacteriaceae Abundance"="#074fa8",
                               'Invisible1' = '#FFFFFF'),                    
                    breaks = c('Microenvironmental Stimulation', 'Enterobacteriaceae Abundance', "Bifidobacteriaceae Abundance")) +  # Specify visible legend entries
  scale_alpha_identity() +  # Use the alpha values as is
  theme(legend.background = element_rect(fill = "white", colour = "black"),legend.position = c(0.65, 0.90)) # Set legend background to white

row_2_pos  = cowplot::plot_grid(p_imp_1_pos,p_imp_2_pos, ncol = 2, rel_widths  = c(1,1.1), label_fontface = "bold")

# SAVE PNG
graphics.off()
png(file = paste0("./FIGURE_4.png"),   # The directory you want to save the file in
    width     = 10,
    height    = 4,
    units     = "in",
    res       = 300)
row_2_pos
dev.off()

