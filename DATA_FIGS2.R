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

options(width = 10000)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

###############################################################################
library(randomForest)
library(pROC)
library(party)
library(scales)  # For scientific notation formatting
library(patchwork)  # For combining plots

results = readRDS('./DATA_VARIABLE_IMPORTANCE.rds') # this was with folds
results = results[c("time_point","i_e_aff","i_e_iga","i_b_tot","i_bc_tot","i_c_tot","i_rnd")]

results = results %>%
  rowwise() %>%
  mutate(
    min_val = min(c(i_e_aff, i_e_iga, i_b_tot, i_bc_tot, i_c_tot, i_rnd)),
    max_val = max(c(i_e_aff, i_e_iga, i_b_tot, i_bc_tot, i_c_tot, i_rnd)),
    i_e_aff  = (i_e_aff  - min_val) / (max_val - min_val),
    i_e_iga  = (i_e_iga  - min_val) / (max_val - min_val),
    i_b_tot  = (i_b_tot  - min_val) / (max_val - min_val),
    i_bc_tot = (i_bc_tot - min_val) / (max_val - min_val),
    i_c_tot  = (i_c_tot  - min_val) / (max_val - min_val),
    i_rnd    = (i_rnd    - min_val) / (max_val - min_val)
  ) %>%
  select(-min_val, -max_val)  # Clean up

# Reshape the data from wide to long format for plotting
long_data = results %>%
  select(time_point, i_e_aff, i_e_iga, i_b_tot, i_bc_tot, i_c_tot, i_rnd) %>%
  pivot_longer(
    cols = c(i_e_aff, i_b_tot, i_e_iga, i_bc_tot, i_c_tot, i_rnd),
    names_to = "predictor",
    values_to = "value"
  ) %>%
  # Create a factor for predictors to control order
  mutate(predictor = factor(predictor, 
                            levels = c("i_e_aff", "i_e_iga", "i_b_tot", "i_bc_tot", "i_c_tot", "i_rnd")))

min_value = min(long_data$value)
max_value = 1  # Maximum is 1 for normalized values
# Create prettier labels for the legend
predictor_labels = c(
  "i_e_aff" = "eSIgA affinity against<br>*Enterobacteriaceae*",
  "i_e_iga" = "eSIgA concentration against<br>*Enterobacteriaceae*",
  "i_b_tot" = "Total *Bifidobacteriaceae*<br>abundance in gut lumen",
  "i_bc_tot" = "Total *Bacteroidaceae*<br>abundance in gut lumen",
  "i_c_tot" = "Total *Clostridiales*<br>abundance in gut lumen",
  "i_rnd" = "Random predictor"
)

library(ggtext)

# 4. Heatmap-style visualization
heatmap_plot = ggplot(long_data, aes(x = as.factor(time_point), y = predictor, fill = value)) +
  geom_tile(color = "gray90", size = 0.5) +
  # Grayscale palette with 1=black and min_value=white
  scale_fill_gradient(low = "white", high = "black", 
                      name = "Value", 
                      limits = c(min_value, max_value)) +
  # Add text labels showing the actual values
  geom_text(aes(label = sprintf("%.2f", value)), 
            color = ifelse(long_data$value > (min_value + max_value)/2, "white", "black"),
            size = 3.5, fontface = "bold") +
  labs(
    title = "Heatmap of normalized predictor importance values<br>in predicting total *Enterobacteriaceae* abundance in gut lumen",
    subtitle = paste0("Scale: ", sprintf("%.2f", min_value), " (white) to 1.00 (black)"),
    x = "Day of Life (DOL)",
    y = "Predictor"
  )+
  theme_minimal() +
  theme(
    axis.text.y = element_markdown(size = 10, lineheight = 1.1),
    panel.grid = element_blank(),
    plot.title = element_markdown(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_y_discrete(labels = predictor_labels)


graphics.off()
png(file = paste0("./FIGURE_S2.png"), 
    width     = 9,
    height    = 5,
    units     = "in",
    res       = 300)
print(heatmap_plot)
dev.off()


