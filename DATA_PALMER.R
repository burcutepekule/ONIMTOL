# Load necessary libraries
library(tidyverse)
library(gridExtra)


# Averaged for subjects 7 and 9 in Fig. 2, Palmer et al. 2007
data_abs = tibble(
  Age_Num = c(2, 7, 30, 90, 180, 1095),  # Age in days
  Enterobacteriaceae = 10^c(7.5,9.8,10.1,10.4,10.4,10.7),
  Bifidobacteriaceae = 10^c(7.5,9.8,10.1,10.4,10.4,10.7),
  Bacteroidaceae = 10^c(7.5,9.8,10.1,10.4,10.4,10.7),
  Clostridales = 10^c(7.5,9.8,10.1,10.4,10.4,10.7)
)

# Function to interpolate each bacteria group
interpolate_bacteria = function(bacteria) {
  df = data_abs %>% select(Age_Num, bacteria)
  names(df)[2] = "Count"
  interpolated = approx(df$Age_Num, df$Count, n = 1095)$y
  return(data.frame(Age_Num = 1:1095, Bacteria = bacteria, Count = interpolated))
}

# Apply the function to each bacteria group and combine results
interpolated_data = bind_rows(
  interpolate_bacteria("Enterobacteriaceae"),
  interpolate_bacteria("Bifidobacteriaceae"),
  interpolate_bacteria("Bacteroidaceae"),
  interpolate_bacteria("Clostridales")
)
interpolated_data_all     = interpolated_data
interpolated_data_1       = interpolated_data

# Find the rows corresponding to day 720 for each bacteria
interpolated_data_1 = interpolated_data_1 %>% filter(Age_Num<=720)
day_720_data        = interpolated_data_1 %>% filter(Age_Num == 720)

# Create a new dataframe for days 721 to 735
additional_days = tibble(
  Age_Num = rep(721:735, times = nrow(day_720_data)),
  Bacteria = rep(day_720_data$Bacteria, each = 15),
  Count = rep(day_720_data$Count, each = 15)
)

# Combine the new data with the original dataframe
interpolated_data_extended = bind_rows(interpolated_data_1, additional_days)
saveRDS(interpolated_data_extended,'interpolated_data_extended.rds')

##### FOR PLOTTING
# Create a new dataframe for days 721 to 735
additional_days = tibble(
  Age_Num = rep(721:1095, times = nrow(day_720_data)),
  Bacteria = rep(day_720_data$Bacteria, each = 375),
  Count = rep(day_720_data$Count, each = 375)
)

# Combine the new data with the original dataframe
interpolated_data_extended_plot = bind_rows(interpolated_data_1, additional_days)

interpolated_data_all$type           = 'data'
interpolated_data_extended_plot$type = 'modified convergence'
interpolated_data_plot = rbind(interpolated_data_all,interpolated_data_extended_plot)
interpolated_data_plot = unique(interpolated_data_plot[c(1,3,4)])

p=ggplot() +
  geom_line(data = interpolated_data_plot, aes(x = Age_Num, y = Count, linetype = type), size = 0.5) +
  scale_linetype_manual(values = c('data' = 'solid', 'modified convergence' = 'dashed')) +
  labs(title = "",
       x = "Days",
       y = "log10(Cells/g Feces)") +
  theme_minimal() + 
  theme(legend.title = element_blank()) + # This line removes the legend title
  scale_y_log10() + 
  scale_x_log10()