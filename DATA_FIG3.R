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
source('./misc/PLOT_ROW_1_SOFT.R')#p_bp, p_bp_dis, p_hm, p_hm_cts
row_1   = cowplot::plot_grid(p_hm_cts,p_bp,ncol = 2,align='h', rel_widths = c(1,1), labels=c('A)','B)'),label_fontface = "bold")

### PREDICTIVE CAPACITY
results = readRDS('logreg_auc_results_Enterobacteriaceae.rds') #AUC results
# Other plotting parameters
span_plot  = 0.05
add_pred   = 0
min_day    = 5
max_day    = 720+add_pred
precision  = 4
taxon_pick = 'Enterobacteriaceae'

source('./misc/PLOT_PREDICTIVE_CAPACTIY_AUC.R') # computes the predictive capacity (seed is fixed for reproducibility)
# # for each time point, output plots: p_e_pred_1, p_e_pred_2, p_e_pred_3

balanced_sampled_data = readRDS('balanced_sampled_data.rds') # load the already subsampled dataset
balanced_sampled_data = balanced_sampled_data %>% dplyr::rowwise() %>% dplyr::mutate(State=ifelse(final_status==1,'Tolerant','Hyperreactive'))
combined_data         = balanced_sampled_data
combined_data_total   = combined_data[c('days','feces_total','State')]
combined_data_iga     = combined_data[c('days','feces_iga_bound','State')]

### p-value funtions 
get_p_value_iga = function(sub_data) {
  
  wilcox_result = wilcox.test(feces_iga_bound ~ State, data = sub_data, exact = FALSE, correct = FALSE)
  d_t = sub_data %>% filter(State == 'Tolerant')
  d_h = sub_data %>% filter(State == 'Hyperreactive')
  
  N1 = dim(d_t)[1]
  N2 = dim(d_h)[1]
  N = nrow(sub_data)
  U = wilcox_result$statistic
  r_rb = 1 - 2 * U / (N1 * N2)
  mean_p_value = wilcox_result$p.value
  mean_r_value = r_rb
  return(c(mean_p_value,mean_r_value))
}

get_p_value_total = function(sub_data) {
  
  wilcox_result = wilcox.test(feces_total ~ State, data = sub_data, exact = FALSE, correct = FALSE)
  d_t = sub_data %>% filter(State == 'Tolerant')
  d_h = sub_data %>% filter(State == 'Hyperreactive')
  
  N1 = dim(d_t)[1]
  N2 = dim(d_h)[1]
  N = nrow(sub_data)
  U = wilcox_result$statistic
  r_rb = 1 - 2 * U / (N1 * N2)
  mean_p_value = wilcox_result$p.value
  mean_r_value = r_rb
  return(c(mean_p_value,mean_r_value))
  
}

p_values_iga = combined_data_iga %>%
  dplyr::group_by(days) %>%
  dplyr::summarise(
    values = get_p_value_iga(cur_data()), # Assuming cur_data() is a placeholder for the actual data passed to the function
    p_value = values[1],
    r_value = values[2]
  ) %>%
  mutate(label = ifelse(p_value < 0.05, sprintf("p = %.2e", p_value), "ns")) # Formatting p-value as scientific notation

p_values_total = combined_data_total %>%
  dplyr::group_by(days) %>%
  dplyr::summarise(
    values = get_p_value_total(cur_data()), # Assuming cur_data() is a placeholder for the actual data passed to the function
    p_value = values[1],
    r_value = values[2]
  ) %>%
  mutate(label = ifelse(p_value < 0.05, sprintf("p = %.2e", p_value), "ns")) # Formatting p-value as scientific notation

# Print to check the correctness of p-values
get_significance_stars = function(p_value) {
  if (p_value < 0.001) {
    return("***")  # Highly significant
  } else if (p_value < 0.01) {
    return("**")   # Very significant
  } else if (p_value < 0.05) {
    return("*")    # Significant
  } else {
    return("ns")     # Not significant
  }
}

# Apply this function to the p_values dataframe
p_values_iga$stars   = sapply(p_values_iga$p_value, get_significance_stars)
combined_data_agg_iga = aggregate(feces_iga_bound~days, combined_data, FUN=max)
colnames(combined_data_agg_iga)[2]='ymax'
p_values_iga = merge(p_values_iga,combined_data_agg_iga,by='days')

p_values_total$stars = sapply(p_values_total$p_value, get_significance_stars)
combined_data_agg_total = aggregate(feces_total~days, combined_data, FUN=max)
colnames(combined_data_agg_total)[2]='ymax'
p_values_total = merge(p_values_total,combined_data_agg_total,by='days')

### Title and axes labels
tit_1 = paste0(taxon_pick, " markers in fecal samples")
tit_2 = paste0(taxon_pick, " markers in fecal samples")
y_1   = "Normalized value"

############################ IGA - BOUND

combined_data_1_iga = combined_data_iga %>% filter(days %in% c(5,10,15,30,45,60))
p_values_1_iga      = p_values_iga %>% filter(days %in% c(5,10,15,30,45,60))
# Plot with significance stars
plot_r2_1_iga = ggplot(combined_data_1_iga, aes(x = factor(days), y = 100*feces_iga_bound, fill = State)) +
  geom_boxplot() +
  # geom_violin(trim=TRUE) +
  geom_text(data = p_values_1_iga, aes(x = factor(days), y = 100*(ymax+0.005), label = stars),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  geom_text(data = p_values_1_iga, aes(x = factor(days), y = 100*(ymax+0.012), label = round(r_value,2)),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  labs(title = tit_2,
       subtitle = "",
       x = "",
       y = y_1) +
  theme_minimal() +
  scale_fill_manual(values = c("Hyperreactive" = "#9A0680", "Tolerant" = "#A3FFD6")) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10))

combined_data_2_iga = combined_data_iga %>% filter(days %in% c(90, 120, 150, 180, 270))
p_values_2_iga      = p_values_iga %>% filter(days %in% c(90, 120, 150, 180, 270))
# Plot with significance stars
plot_r2_2_iga = ggplot(combined_data_2_iga, aes(x = factor(days), y = 100*feces_iga_bound, fill = State)) +
  geom_boxplot() +
  # geom_violin(trim=TRUE) +
  geom_text(data = p_values_2_iga, aes(x = factor(days), y = 100*(ymax+0.0008), label = stars),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  geom_text(data = p_values_2_iga, aes(x = factor(days), y = 100*(ymax+0.0020), label = round(r_value,2)),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  labs(title = "",
       subtitle = "",
       x = "DOL",
       y = "") +
  theme_minimal() +
  scale_fill_manual(values = c("Hyperreactive" = "#9A0680", "Tolerant" = "#A3FFD6")) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10))

combined_data_3_iga = combined_data_iga %>% filter(days %in% c(360, 450, 540, 630, seq(720,max_day,180)))
p_values_3_iga      = p_values_iga %>% filter(days %in% c(360, 450, 540, 630, seq(720,max_day,180)))
# Plot with significance stars
plot_r2_3_iga = ggplot(combined_data_3_iga, aes(x = factor(days), y = 100*feces_iga_bound, fill = State)) +
  geom_boxplot() +
  # geom_violin(trim=TRUE) +
  geom_text(data = p_values_3_iga, aes(x = factor(days), y = 100*(ymax+0.0003), label = stars),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  geom_text(data = p_values_3_iga, aes(x = factor(days), y = 100*(ymax+0.0008), label = round(r_value,2)),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  labs(title = "",
       subtitle = "",
       x = "DOL",
       y = "") +
  theme_minimal() +
  scale_fill_manual(values = c("Hyperreactive" = "#9A0680", "Tolerant" = "#A3FFD6")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10))

############################ TOTAL 

combined_data_1_total = combined_data_total %>% filter(days %in% c(5,10,15,30,45,60))
p_values_1_total      = p_values_total %>% filter(days %in% c(5,10,15,30,45,60))
# Plot with significance stars
plot_r2_1_total = ggplot(combined_data_1_total, aes(x = factor(days), y = 100*feces_total, fill = State)) +
  geom_boxplot() +
  # geom_violin(trim=TRUE) +
  geom_text(data = p_values_1_total, aes(x = factor(days), y = 100*(ymax+0.020), label = stars),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  geom_text(data = p_values_1_total, aes(x = factor(days), y = 100*(ymax+0.035), label = round(r_value,2)),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  labs(title = tit_1,
       subtitle = "",
       x = "",
       y = y_1) +
  theme_minimal() +
  scale_fill_manual(values = c("Hyperreactive" = "#9A0680", "Tolerant" = "#A3FFD6")) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10))

combined_data_2_total = combined_data_total %>% filter(days %in% c(90, 120, 150, 180, 270))
p_values_2_total      = p_values_total %>% filter(days %in% c(90, 120, 150, 180, 270))
# Plot with significance stars
plot_r2_2_total = ggplot(combined_data_2_total, aes(x = factor(days), y = 100*feces_total, fill = State)) +
  geom_boxplot() +
  # geom_violin(trim=TRUE) +
  geom_text(data = p_values_2_total, aes(x = factor(days), y = 100*(ymax+0.002), label = stars),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  geom_text(data = p_values_2_total, aes(x = factor(days), y = 100*(ymax+0.005), label = round(r_value,2)),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  labs(title = "",
       subtitle = "",
       x = "DOL",
       y = "") +
  theme_minimal() +
  scale_fill_manual(values = c("Hyperreactive" = "#9A0680", "Tolerant" = "#A3FFD6")) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10))

combined_data_3_total = combined_data_total %>% filter(days %in% c(360, 450, 540, 630, seq(720,max_day,180)))
p_values_3_total      = p_values_total %>% filter(days %in% c(360, 450, 540, 630, seq(720,max_day,180)))
# Plot with significance stars
plot_r2_3_total = ggplot(combined_data_3_total, aes(x = factor(days), y = 100*feces_total, fill = State)) +
  geom_boxplot() +
  # geom_violin(trim=TRUE) +
  geom_text(data = p_values_3_total, aes(x = factor(days), y = 100*(ymax+0.0006), label = stars),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  geom_text(data = p_values_3_total, aes(x = factor(days), y = 100*(ymax+0.002), label = round(r_value,2)),
            inherit.aes = FALSE, vjust = 0, color = "black") +  # Adjust vjust as needed
  labs(title = "",
       subtitle = "",
       x = "DOL",
       y = "") +
  theme_minimal() +
  scale_fill_manual(values = c("Hyperreactive" = "#9A0680", "Tolerant" = "#A3FFD6")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 10))

######### PREDICTIVE CAPACITY ##################################################
specific_days <- c(5,10,15,30,45,60)
stats_merged_plot <- stats_merged %>% dplyr::filter(days >= min(specific_days))

stats_merged_plot <- stats_merged_plot %>% 
  dplyr::mutate(Abundance = factor(Abundance, levels = c('Total', 'SIgA-bound', 'Total+SIgA-bound')))

# Define your specific days and their new positions
new_positions <- seq_along(specific_days)-1

# Create a data frame of these mappings
day_mapping <- data.frame(actual_day = specific_days, plot_day = new_positions)

# Add interpolated plot days for all days in your dataset
stats_merged_plot <- stats_merged_plot %>%
  dplyr::mutate(plot_day = approx(specific_days, new_positions, xout = days)$y)

# Plot for all days
p_e_pred_1a <- ggplot(stats_merged_plot, aes(x = plot_day)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper, fill = Abundance), alpha = 0.15) +
  geom_line(aes(y = smooth_median, color = Abundance, linetype = Abundance)) +
  scale_color_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_fill_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed")) +
  scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_continuous(breaks = new_positions, 
                     labels = specific_days, 
                     limits=c(min(new_positions),max(new_positions)),
                     expand = expansion(mult = c(0.12, 0.12), add = c(0.05, 0.025))) +
  labs(legend.position = "none", 
       title = paste0("Predictive Power of ",taxon_pick," Abundance over Time"),
       x = "",
       y = "Predictive Power (AUC)") +
  theme_minimal() + theme(legend.position = "none")


specific_days <- c(90, 120, 150, 180, 270)
stats_merged_plot <- stats_merged %>% dplyr::filter(days >= min(specific_days))

stats_merged_plot <- stats_merged_plot %>% 
  dplyr::mutate(Abundance = factor(Abundance, levels = c('Total', 'SIgA-bound', 'Total+SIgA-bound')))

# Define your specific days and their new positions
new_positions <- seq_along(specific_days)-1

# Create a data frame of these mappings
day_mapping <- data.frame(actual_day = specific_days, plot_day = new_positions)

# Add interpolated plot days for all days in your dataset
stats_merged_plot <- stats_merged_plot %>%
  dplyr::mutate(plot_day = approx(specific_days, new_positions, xout = days)$y)

# Plot for all days
p_e_pred_1b <- ggplot(stats_merged_plot, aes(x = plot_day)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper, fill = Abundance), alpha = 0.15) +
  geom_line(aes(y = smooth_median, color = Abundance, linetype = Abundance)) +
  scale_color_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_fill_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed")) +
  scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_continuous(breaks = new_positions, 
                     labels = specific_days, 
                     limits=c(min(new_positions),max(new_positions)),
                     expand = expansion(mult = c(0.12, 0.12), add = c(0.1, 0.05))) +
  labs(legend.position = "none", 
       title = "",
       x = "DOL",
       y = "") +
  theme_minimal() + theme(legend.position = "none") + theme(axis.text.y = element_blank())


specific_days <- c(360, 450, 540, 630, seq(720,max_day,180))
stats_merged_plot <- stats_merged %>% dplyr::filter(days >= min(specific_days))

stats_merged_plot <- stats_merged_plot %>% 
  dplyr::mutate(Abundance = factor(Abundance, levels = c('Total', 'SIgA-bound', 'Total+SIgA-bound')))

# Define your specific days and their new positions
new_positions <- seq_along(specific_days)-1

# Create a data frame of these mappings
day_mapping <- data.frame(actual_day = specific_days, plot_day = new_positions)

# Add interpolated plot days for all days in your dataset
stats_merged_plot <- stats_merged_plot %>%
  dplyr::mutate(plot_day = approx(specific_days, new_positions, xout = days)$y)

# Plot for all days
p_e_pred_1c <- ggplot(stats_merged_plot, aes(x = plot_day)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper, fill = Abundance), alpha = 0.15) +
  geom_line(aes(y = smooth_median, color = Abundance, linetype = Abundance)) +
  scale_color_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_fill_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed")) +
  scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_continuous(breaks = new_positions, 
                     labels = specific_days, 
                     limits=c(min(new_positions),max(new_positions)),
                     expand = expansion(mult = c(0.2, 0.05), add = c(0, 0))) +
  labs(title = "",
       x = "DOL",
       y = "") +
  theme_minimal() + theme(axis.text.y = element_blank())  # This removes y-axis text

specific_days <- c(5, 10, 15, 30, 45, 60, 90, 120, 150, 180, 270, 360, 450, 540, 630, seq(720,max_day,180))
stats_merged_plot <- stats_merged %>% dplyr::filter(days >= min(specific_days))

stats_merged_plot <- stats_merged_plot %>% 
  dplyr::mutate(Abundance = factor(Abundance, levels = c('Total', 'SIgA-bound', 'Total+SIgA-bound')))

# Define your specific days and their new positions
new_positions <- seq_along(specific_days)-1

# modify to align
# new_positions = c(seq(from = 1, to = 4.35, length.out = 6),  seq(from = 6, to = 9.35, length.out = 5), seq(from = 11.25, to = 14.5, length.out = 5))
new_positions = c(seq(from = 1, to = 4.35, length.out = 6),  seq(from = 6, to = 9.35, length.out = 5), seq(from = 11.25, to = 14.5, length.out = length(specific_days)-11))
# Create a data frame of these mappings
day_mapping <- data.frame(actual_day = specific_days, plot_day = new_positions)

# Add interpolated plot days for all days in your dataset
stats_merged_plot <- stats_merged_plot %>%
  dplyr::mutate(plot_day = approx(specific_days, new_positions, xout = days)$y)

# Plot for all days
p_e_pred_1d <- ggplot(stats_merged_plot, aes(x = plot_day)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper, fill = Abundance), alpha = 0.15) +
  geom_line(aes(y = smooth_median, color = Abundance, linetype = Abundance)) +
  scale_color_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_fill_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed")) +
  scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_continuous(breaks = new_positions, 
                     labels = specific_days, 
                     limits=c(min(new_positions),max(new_positions)),
                     expand = expansion(mult = c(0, 0), add = c(0.5, 0.2))) +
  labs(title = paste0("Predictive Power of ",taxon_pick," Abundance over Time"),
       x = "",
       y = "Predictive Power (AUC)") +
  theme_minimal() 


row_2 = cowplot::plot_grid(plot_r2_1_total,plot_r2_2_total,plot_r2_3_total, ncol = 3,align='h', rel_widths = c(0.5,0.5,0.7))
row_3 = cowplot::plot_grid(plot_r2_1_iga,plot_r2_2_iga,plot_r2_3_iga, ncol = 3,align='h', rel_widths = c(0.5,0.5,0.7))
row_4 = cowplot::plot_grid(p_e_pred_1d, ncol = 1,align='h', rel_widths = c(1))
nested_plot  = cowplot::plot_grid(row_1,row_2,row_3,row_4,nrow = 4,align='v', rel_heights  = c(1.2,1,1,1),labels=c('','C)','D)','E)'),label_fontface = "bold")
graphics.off()
png(file = paste0("./FIGURE_3.png"), 
    width     = 12,
    height    = 14,
    units     = "in",
    res       = 300)
print(nested_plot)
dev.off()



