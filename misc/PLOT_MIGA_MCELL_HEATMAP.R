
colnames(keep_res_use)[2:3]=c('mIgA','aff_final')

# # Assuming minimum value in data, 1, and maximum value in data
# min_val <- min(keep_res_use$aff_final)
# max_val <- max(keep_res_use$aff_final)
# mid_val <- 1  # The value you want to show as gray95

# Calculate normalized positions for these values
# Need to normalize to 0-1 scale for the values parameter
norm_min <- 0
norm_mid <- (mid_val - min_val) / (max_val - min_val)
norm_max <- 1

affs_plot = round(unique(keep_res_use$mIgA),2)
# Then update your plot
plot_out = ggplot(keep_res_use, aes(x = coeff_binding, y = mIgA)) +
  geom_tile(aes(fill = aff_final), color = '#d1e2e2') +
  scale_fill_gradientn(
    colors = c("#0000CD", "gray95", "#E52B50"),
    values = c(norm_min, norm_mid, norm_max),
    limits = c(min_val, max_val),  # Set the color scale limits with max of 26
    name = "Affinity"
  ) +
  # Rest of your code remains the same
  scale_y_log10(breaks = affs_plot,
                labels = affs_plot) +
  scale_x_log10(breaks = c(1,2,5,10,20,50,100),
                labels = c(1,2,5,10,20,50,100)) +
  theme_minimal() +
  labs(
    x = expression(paste(epsilon^"m", "/", epsilon^"uc")),
    y = "mSIgA Affinity", 
    title = taxon_name
  ) +
  theme(plot.title = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))

# plot_out = ggplot(keep_res_use, aes(x = coeff_binding, y = mIgA)) +
#   geom_tile(aes(fill = aff_final), color = '#d1e2e2') +
#   scale_fill_gradientn(
#     colors = c("#0000CD", "gray95", "#E52B50"),
#     values = c(norm_min, norm_mid, norm_max),
#     limits = c(min_val, max_val),
#     name = "Affinity"
#   )+
#   scale_y_log10(breaks = affs_plot,
#                 labels = affs_plot) +
#   scale_x_log10(breaks = c(1,2,5,10,20,50,100),
#                 labels = c(1,2,5,10,20,50,100))+
#   theme_minimal() +
#   labs(
#     x = expression(paste(epsilon^"m", "/", epsilon^"uc")),
#     y = "mSIgA Affinity", 
#     title = taxon_name
#   ) +
#   theme(plot.title = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))
