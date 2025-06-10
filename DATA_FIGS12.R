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
library(ggtext)
options(width = 10000)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
#################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# LOAD AUC VALUES FROM REGRESSION
auc_values_keep = readRDS('AUC_MODEL.rds')
auc_values_data = readRDS('AUC_DATA.rds')

# REMOVE AUC=1 
sims_with_auc_1 <- auc_values_keep %>%
  group_by(sim) %>%
  filter(any(AUC == 1)) %>%
  distinct(sim) %>%
  pull(sim)
filtered_auc_values <- auc_values_keep %>%
  filter(!(sim %in% sims_with_auc_1))

sim_stats_filtered = filtered_auc_values %>%
  dplyr::group_by(Model) %>%
  dplyr::summarize(
    mean_AUC = mean(AUC),
    lower_CI_mean = mean(LowerCI),
    upper_CI_mean = mean(UpperCI)
  )
sim_stats_filtered = sim_stats_filtered[order(sim_stats_filtered$mean_AUC, decreasing=TRUE),]
month_pick         = 6

auc_values_model = sim_stats_filtered
colnames(auc_values_model)[2:4]=c('AUC','LowerCI','UpperCI')
auc_values_model$source = 'Model'
auc_values_model$rank   = rank(-auc_values_model$AUC)  # negative sign to rank in descending order
auc_values_model$month  = month_pick

auc_values_merged = as.data.frame(rbind(auc_values_model,auc_values_data))

auc_values_pick = auc_values_merged
auc_values_pick = auc_values_pick[c('Model','AUC','source','rank','month','LowerCI','UpperCI')]
all_values_use  = auc_values_pick
colnames(all_values_use)[1] =c('taxon')

# Add the mapping to your dataframe
all_values_compare = all_values_use

rank_comparison = all_values_compare %>%
  select(taxon, month, source, rank) %>%
  pivot_wider(names_from = source, values_from = rank)

# Add a column for rank difference
rank_comparison$rank_diff = rank_comparison$Data - rank_comparison$Model

rank_dot_data = all_values_compare %>% select(taxon, month, source, rank)

p1 = ggplot(rank_dot_data, aes(x = rank, y = reorder(taxon, -rank), color = source)) +
  # Add position_dodge to separate overlapping points
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  # Also apply the same dodge to the lines
  geom_line(aes(group = interaction(taxon, month)), 
            position = position_dodge(width = 0.3),
            color = "gray") +
  scale_color_manual(values = c("Data" = "#3498db", "Model" = "#e74c3c")) +
  scale_x_continuous(breaks = 1:5) +
  labs(title = "Comparison of Rankings: Data vs Model",
       # subtitle = "Connected points show rank difference for the same taxonomic group",
       x = "Rank",
       y = "Taxonomic Group",
       color = "Source") +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# Create a faceted plot that shows direct comparisons for each taxonomic group by month
p2 = ggplot(all_values_compare, aes(x = source, y = AUC, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  # Add error bars for confidence intervals
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), 
                position = position_dodge(width = 0.9),
                width = 0.25,  # Width of the error bars
                color = "black") +
  # Position text above the UpperCI
  geom_text(aes(y = UpperCI, label = sprintf("%.2f", AUC)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,  # Negative value moves text upward from the anchor point
            size = 3) +
  facet_grid(~ taxon, labeller = labeller(
    taxon = function(x) gsub("\\s", "\n", x)  # Add line breaks for long names
  )) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +  # Reference line
  scale_fill_manual(values = c("Data" = "#3498db", "Model" = "#e74c3c")) +
  labs(title = "Comparison of AUC Values: Model vs Data",
       x = "",
       y = "AUC (Area Under ROC Curve)",
       fill = "Source") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray90")
  ) +
  ylim(0, 1)  # Set y-axis limits from 0 to 1 for AUC values

row_final2 = cowplot::plot_grid(p1, p2,labels=c('A)','B)'), ncol=2, rel_widths = c(0.70,1), label_fontface = "bold")

graphics.off()
png(file ="./FIGURE_S12.png",    # The directory you want to save the file in
    width     = 12,
    height    = 5,
    units     = "in",
    res       = 300)
print(row_final2)
dev.off()

