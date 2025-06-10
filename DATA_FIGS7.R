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
library(latex2exp)

options(width = 10000)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
#################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

melted_data = readRDS('./DATA_MORRIS.rds') # read Morris table

p=ggplot(melted_data, aes(x = parameter, y = mu.star)) +
  geom_point(size = 3) +
  labs(
    x = "Parameter",
    y = TeX("$\\mu^{*}$"),  # LaTeX syntax for mu*
    title = "Morris Sensitivity Analysis (Mu Star Values)",
    color = "Method"
  ) +
  scale_x_discrete(labels = c(
    "kappa_e" = TeX("$\\kappa_{E}$"),
    "alpha_b" = TeX("$\\alpha_{B}$"),
    "alpha_bc" = TeX("$\\alpha_{BC}$"),
    "alpha_c" = TeX("$\\alpha_{C}$"),
    "mcell_th" = TeX("$t^{m}$"),
    "range_apop" = TeX("$th_{apop}$"),
    "range_tf" = TeX("$th_{range}$"),
    "Ag_dependent_activation_m" = TeX("$\\psi^{m}$"),
    "tau_new" = TeX("$\\tau^{new}$"),
    "tau_c" = TeX("$\\tau^{c}$"),
    "tau_delta" = TeX("$\\tau^{\\delta}$"),
    "C_n" = TeX("$C_{n}$"), 
    "C_I" = TeX("$C_{I}$"),
    "c_n" = TeX("$c_{n}$"), 
    "coeff_binding" = TeX("$\\epsilon^{m}/\\epsilon^{uc}$")
  )) +
  theme_minimal() +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 14, family = "serif"),
        axis.title.y  = element_text(size = 18, family = "serif"))

# Get summary statistics
summary_stats <- melted_data %>%
  summarize(
    min = min(mu.star),
    median = median(mu.star),
    mean = mean(mu.star),
    q25 = quantile(mu.star, 0.25),
    q50 = quantile(mu.star, 0.50),
    q75 = quantile(mu.star, 0.75),
    q95 = quantile(mu.star, 0.95),
    max = max(mu.star)
  )

q25_threshold <- summary_stats$q25
q50_threshold <- summary_stats$q50
q75_threshold <- summary_stats$q75
q95_threshold <- summary_stats$q95

# p = p + geom_hline(yintercept = q25_threshold, linetype = "dashed", color = "red", size = 1) +
#   annotate("text", x = 12, y = q25_threshold, label = paste("Q25 =", round(q25_threshold, 3)),
#            vjust = -0.5, hjust = 0, color = "red")

p = p + geom_hline(yintercept = q50_threshold, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 12, y = q50_threshold, label = paste("Q50 =", round(q50_threshold, 3)),
           vjust = -0.5, hjust = 0, color = "red")

# p = p + geom_hline(yintercept = q75_threshold, linetype = "dashed", color = "red", size = 1) +
#   annotate("text", x = 12, y = q75_threshold, label = paste("Q75 =", round(q75_threshold, 3)),
#            vjust = -0.5, hjust = 0, color = "red")

p = p + geom_hline(yintercept = q95_threshold, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 12, y = q95_threshold, label = paste("Q95 =", round(q95_threshold, 3)),
           vjust = -0.5, hjust = 0, color = "red")

graphics.off()
png(file =paste0("./FIGURE_S7.png"),    # The directory you want to save the file in
    width     = 7,
    height    = 4,
    units     = "in",
    res       = 300)
print(p)
dev.off()

melted_data$mu      = round(melted_data$mu,4)
melted_data$mu.star = round(melted_data$mu.star,4)
melted_data$sigma   = round(melted_data$sigma,4)




