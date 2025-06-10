library(smotefamily)
library(pROC)

################################################################################
# Visualization section
################################################################################

# Create long format data for plotting
results_long_auc = tidyr::pivot_longer(
  results,
  cols = c(auc_total, auc_iga, auc_both),
  names_to = "predictor",
  values_to = "auc"
)

# Add standard error columns to long format
results_long_auc$se = NA
results_long_auc$se[results_long_auc$predictor == "auc_total"] = results$se_auc_total
results_long_auc$se[results_long_auc$predictor == "auc_iga"] = results$se_auc_iga
results_long_auc$se[results_long_auc$predictor == "auc_both"] = results$se_auc_both

# Rename predictors for better labels
results_long_auc$predictor = factor(results_long_auc$predictor,
                                    levels = c("auc_total", "auc_iga", "auc_both"),
                                    labels = c("Feces Total", "Feces IgA Bound", "Both"))

# Plot AUC vs time
p1 = ggplot(results_long_auc, aes(x = time_point, y = auc, color = predictor, group = predictor)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = auc - se, ymax = auc + se, fill = predictor), 
              alpha = 0.2, color = NA) +
  labs(title = "Predictive Power (AUC) at Each Time Point",
       x = "Days",
       y = "AUC (Area Under ROC Curve)",
       color = "Predictor",
       fill = "Predictor") +
  theme_minimal() +
  theme(legend.position = "bottom")

################################################################################
# Smoothed visualization for different time periods as in original code
################################################################################

# First convert column names to match the original analysis
colnames(results_long_auc)[which(colnames(results_long_auc) == "time_point")] = "days"

# Filter and process data
stats_total = results_long_auc %>% dplyr::filter(predictor == 'Feces Total')
stats_iga_alone = results_long_auc %>% dplyr::filter(predictor == 'Feces IgA Bound')
stats_iga = results_long_auc %>% dplyr::filter(predictor == 'Both')

# Rename columns to match original code
colnames(stats_total)[which(colnames(stats_total) == "auc")] = "median"
colnames(stats_iga_alone)[which(colnames(stats_iga_alone) == "auc")] = "median"
colnames(stats_iga)[which(colnames(stats_iga) == "auc")] = "median"

# Calculate lower and upper bounds
stats_total = stats_total %>% dplyr::rowwise() %>% dplyr::mutate(lower = median - se)
stats_total = stats_total %>% dplyr::rowwise() %>% dplyr::mutate(upper = median + se)

stats_iga_alone = stats_iga_alone %>% dplyr::rowwise() %>% dplyr::mutate(lower = median - se)
stats_iga_alone = stats_iga_alone %>% dplyr::rowwise() %>% dplyr::mutate(upper = median + se)

stats_iga = stats_iga %>% dplyr::rowwise() %>% dplyr::mutate(lower = median - se)
stats_iga = stats_iga %>% dplyr::rowwise() %>% dplyr::mutate(upper = median + se)

# Remove unneeded columns
stats_total = stats_total %>% dplyr::select(-se)
stats_iga_alone = stats_iga_alone %>% dplyr::select(-se)
stats_iga = stats_iga %>% dplyr::select(-se)

# Function to apply loess smoothing and extract smoothed predictions
smooth_values = function(data, x, y) {
  # Skip smoothing if insufficient data points
  if (length(x) < 4) {
    return(y)  # Return original values if not enough for smoothing
  }
  
  tryCatch({
    loess_fit = loess(y ~ x, data = data.frame(x = x, y = y), span = span_plot)
    predict(loess_fit, data.frame(x = x))
  }, error = function(e) {
    # Return original values if smoothing fails
    warning(paste("Smoothing failed:", e$message))
    return(y)
  })
}

smooth_values_iga = function(data, x, y) {
  # Skip smoothing if insufficient data points
  if (length(x) < 4) {
    return(y)  # Return original values if not enough for smoothing
  }
  
  tryCatch({
    loess_fit = loess(y ~ x, data = data.frame(x = x, y = y), span = span_plot)
    predict(loess_fit, data.frame(x = x))
  }, error = function(e) {
    # Return original values if smoothing fails
    warning(paste("Smoothing failed:", e$message))
    return(y)
  })
}

# Apply smoothing with error handling
stats_total$smooth_median = smooth_values_iga(stats_total, stats_total$days, stats_total$median)
stats_total$smooth_lower = smooth_values_iga(stats_total, stats_total$days, stats_total$lower)
stats_total$smooth_upper = smooth_values_iga(stats_total, stats_total$days, stats_total$upper)

stats_iga$smooth_median = smooth_values_iga(stats_iga, stats_iga$days, stats_iga$median)
stats_iga$smooth_lower = smooth_values_iga(stats_iga, stats_iga$days, stats_iga$lower)
stats_iga$smooth_upper = smooth_values_iga(stats_iga, stats_iga$days, stats_iga$upper)

stats_iga_alone$smooth_median = smooth_values_iga(stats_iga_alone, stats_iga_alone$days, stats_iga_alone$median)
stats_iga_alone$smooth_lower = smooth_values_iga(stats_iga_alone, stats_iga_alone$days, stats_iga_alone$lower)
stats_iga_alone$smooth_upper = smooth_values_iga(stats_iga_alone, stats_iga_alone$days, stats_iga_alone$upper)

# Add group identifiers
stats_total$Abundance = 'Total'
stats_iga$Abundance = 'Total+SIgA-bound'
stats_iga_alone$Abundance = 'SIgA-bound'

# Combine datasets
stats_merged = rbind(stats_total, stats_iga_alone, stats_iga)
stats_merged = as.data.frame(stats_merged)

# Define minimum day constant if not defined (using first day as fallback)
if (!exists("min_day")) {
  min_day = min(stats_merged$days)
  cat("Setting min_day to", min_day, "\n")
}

# Plot early days
specific_days = min_day:15
stats_merged_plot = stats_merged %>% dplyr::filter(days <= max(specific_days))
stats_merged_plot = stats_merged_plot %>% 
  dplyr::mutate(Abundance = factor(Abundance, levels = c('Total', 'SIgA-bound', 'Total+SIgA-bound')))

# Define new positions for the x-axis
new_positions = seq_along(specific_days)

# Create a data frame of these mappings
day_mapping = data.frame(actual_day = specific_days, plot_day = new_positions)

# Add interpolated plot days for all days in your dataset
stats_merged_plot = stats_merged_plot %>%
  dplyr::mutate(plot_day = approx(specific_days, new_positions, xout = days)$y)

# Plot for early days
p_e_pred_1 = ggplot(stats_merged_plot, aes(x = plot_day)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper, fill = Abundance), alpha = 0.15) +
  geom_line(aes(y = smooth_median, color = Abundance, linetype = Abundance)) +
  scale_color_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_fill_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed")) +
  scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_continuous(breaks = new_positions, labels = specific_days) +
  labs(title = paste0("Predictive Power of ",taxon_pick," Abundance over Time"),
       x = "Days of Life",
       y = "Predictive Power (AUC)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for later days
specific_days <- c(30, 45, 60, 90, 120, 150, 180, 270, 360, 450, 540, 630, 720, 900, 1080, 1260, 1440)
specific_days = specific_days[which(specific_days<max_day)]
stats_merged_plot = stats_merged %>% dplyr::filter(days >= min(specific_days))
stats_merged_plot = stats_merged_plot %>% 
  dplyr::mutate(Abundance = factor(Abundance, levels = c('Total', 'SIgA-bound', 'Total+SIgA-bound')))

# Define your specific days and their new positions
new_positions = seq_along(specific_days)

# Create a data frame of these mappings
day_mapping = data.frame(actual_day = specific_days, plot_day = new_positions)

# Add interpolated plot days for all days in your dataset
stats_merged_plot = stats_merged_plot %>%
  dplyr::mutate(plot_day = approx(specific_days, new_positions, xout = days)$y)

# Plot for later days
p_e_pred_2 = ggplot(stats_merged_plot, aes(x = plot_day)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper, fill = Abundance), alpha = 0.15) +
  geom_line(aes(y = smooth_median, color = Abundance, linetype = Abundance)) +
  scale_color_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_fill_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed")) +
  scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_continuous(breaks = new_positions, labels = specific_days) +
  labs(title = paste0("Predictive Power of ",taxon_pick," Abundance over Time"),
       x = "DOL",
       y = "Predictive Power (AUC)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for all days
specific_days <- unique(c(min_day, 5, 10, 15, 30, 45, 60, 90, 120, 150, 180, 270, 360, 450, 540, 630, 720, 900, 1080, 1260, 1440))
specific_days = specific_days[which(specific_days<=max_day)]
stats_merged_plot = stats_merged %>% dplyr::filter(days >= min(specific_days))
stats_merged_plot = stats_merged_plot %>% 
  dplyr::mutate(Abundance = factor(Abundance, levels = c('Total', 'SIgA-bound', 'Total+SIgA-bound')))

# Define your specific days and their new positions
new_positions = seq_along(specific_days)-1

# Create a data frame of these mappings
day_mapping = data.frame(actual_day = specific_days, plot_day = new_positions)

# Add interpolated plot days for all days in your dataset
stats_merged_plot = stats_merged_plot %>%
  dplyr::mutate(plot_day = approx(specific_days, new_positions, xout = days)$y)

# Plot for all days
p_e_pred_3 = ggplot(stats_merged_plot, aes(x = plot_day)) +
  geom_ribbon(aes(ymin = smooth_lower, ymax = smooth_upper, fill = Abundance), alpha = 0.15) +
  geom_line(aes(y = smooth_median, color = Abundance, linetype = Abundance)) +
  scale_color_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_fill_manual(values = c("#0C2D57", 'blue', "#f20026")) +
  scale_linetype_manual(values = c("solid", "dotdash", "dashed")) +
  scale_y_continuous(limits = c(0.5, 1.01)) +
  scale_x_continuous(breaks = new_positions, 
                     labels = specific_days, 
                     limits=c(min(new_positions),max(new_positions)),
                     expand = c(0, 0)) +  # This removes the padding
  labs(title = paste0("Predictive Power of ",taxon_pick," Abundance over Time"),
       x = "DOL",
       y = "Predictive Power (AUC)") +
  theme_minimal() 
