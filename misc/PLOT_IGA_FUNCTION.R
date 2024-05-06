if (exists('mIgA_in_keep')){
  keep_miga_vec = mIgA_in_keep
}else{
  keep_miga_vec = data_list_sim$theta_in[(numTaxa + numAgnostic + numGnostic + 1):(numTaxa + numAgnostic + numGnostic + numTaxa)];
}

library(reshape2)

# Create a sequence of x values
x_values <- seq(0.05, 15, by = 0.05)

#diassoctaion rate for omega_i?

d_r = theta[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+17)]

# Calculate the values of the functions
df <- data.frame(
  x = x_values,
  y1 = 1 / (1 + x_values),
  y2 = 1 - 1 / (1 + x_values),
  w  = (1-d_r/(1+x_values))
)

# Reshape the data for plotting with ggplot
df_long <- melt(df, id.vars = 'x')

variable       = 'plasma_affinity'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  dplyr::mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
                taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array), 
                              taxa_array[taxa_index], 
                              taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

out_df_aff = y_out_df_use_l
out_df_e   = y_out_df_use_l %>% filter(taxa=='Enterobacteriaceae')
out_df_ne  = y_out_df_use_l %>% filter(taxa!='Enterobacteriaceae')
out_df_ene = y_out_df_use_l

# Define the x positions for vertical lines
out_df_aff_conv = out_df_aff %>% filter(days==max(out_df_aff$days))
val_high = out_df_aff_conv$median[1]
val_low  = mean(out_df_aff_conv$median[2:4])
vlines <- (c(val_low, val_high))

# Calculate annotation points
annotate_points <- data.frame(
  x = rep(vlines, each = 2),
  y = c(1 / (1 + vlines[1]), 1 - 1 / (1 + vlines[1]), 1 / (1 + vlines[2]), 1 - 1 / (1 + vlines[2])),
  label = sprintf("%.2f", c(1 / (1 + vlines[1]), 1 - 1 / (1 + vlines[1]), 1 / (1 + vlines[2]), 1 - 1 / (1 + vlines[2])))
)

df_long_short_y1      = df_long %>% filter(variable=='y1')
df_long_short_y1_high = df_long_short_y1 %>% filter(x<=1)
df_long_short_y1_low  = df_long_short_y1 %>% filter(x>=1)

df_long_short_y2      = df_long %>% filter(variable=='y2')
df_long_short_y2_high = df_long_short_y2 %>% filter(x<=1)
df_long_short_y2_low  = df_long_short_y2 %>% filter(x>=1)

df_long_short_w       = df_long %>% filter(variable=='w')

annotate_points       = annotate_points[c(1,4),]
annotate_points$label = as.numeric(annotate_points$label)
annotate_points$label = paste0(100*(annotate_points$label), "%")
# annotate_points$y = annotate_points$y + 0.03
annotate_points$y = c(0.95,0.95)
annotate_points$x[1] = annotate_points$x[1]+0.01
annotate_points$x[2] = annotate_points$x[2]-0.2


out_df_ene_last = out_df_ene %>% filter(days==735)

# Define the desired order of taxa
desired_order <- c("Enterobacteriaceae", "Bifidobacteriaceae", "Bacteroidaceae", "Clostridiales")

# Convert the 'taxa' column to a factor with levels in the desired order
out_df_ene_last$taxa <- factor(out_df_ene_last$taxa, levels = desired_order)

# Arrange the dataframe by the ordered factor
out_df_ene_last_ordered <- out_df_ene_last %>%arrange(taxa)

out_df_ene_last = out_df_ene_last_ordered


###### THIS IS sSIgA
# out_df_ene_last$val = 0
# out_df_ene_last[1,]$val = round(1-1/(1+out_df_ene_last[1,]$median),3)  
# out_df_ene_last[2,]$val = round(1/(1+out_df_ene_last[2,]$median),3)  
# out_df_ene_last[3,]$val = round(1/(1+out_df_ene_last[3,]$median),3)  
# out_df_ene_last[4,]$val = round(1/(1+out_df_ene_last[4,]$median),3)  

### PLOT mSIgA!
out_df_ene_last$val = 0
out_df_ene_last[1,]$val = round(1-1/(1+keep_miga_vec[1]),3)  
out_df_ene_last[2,]$val = round(1/(1+keep_miga_vec[2]),3)  
out_df_ene_last[3,]$val = round(1/(1+keep_miga_vec[3]),3)  
out_df_ene_last[4,]$val = round(1/(1+keep_miga_vec[4]),3)  



out_df_ene_last$point_color = c('#f20026','#074fa8','#2bc5d8','#ff9aa1')

# Assuming p_iga_function is your initial ggplot object setup
# p_iga_function_merged = ggplot() +
#   geom_rect(aes(xmin = min(x_values), xmax = 1, ymin = 0.5, ymax = 1), fill = "gray", alpha = 0.1) +
#   geom_rect(aes(xmin = 1, xmax = max(x_values), ymin = 0.5, ymax = 1), fill = "black", alpha = 0.1) +
#   labs(title = "Rate of SIgA Function", x = "Affinity", y = "") +
#   geom_line(data = df_long_short_y1_high, aes(x = x, y = value, color = variable), alpha = 1, size = 1) +
#   geom_line(data = df_long_short_y2_high, aes(x = x, y = value, color = variable), alpha = 0.8, linetype = 'dotted', size = 1) +
#   geom_line(data = df_long_short_y1_low, aes(x = x, y = value, color = variable), alpha = 0.8, linetype = 'dotted', size = 1) +
#   geom_line(data = df_long_short_y2_low, aes(x = x, y = value, color = variable), alpha = 1, size = 1) +
#   geom_line(data = df_long_short_w,  aes(x = x, y = value, color = variable), linetype = 'dashed', alpha = 1, size=1) +
#   scale_x_log10() + 
#   theme_minimal() +  
#   theme(plot.title = element_text(size = 12), legend.position = "none") +
#   scale_color_manual(values = c("#FF6333","gray", "black"), labels = c("" ,"Coating (C)", "Neutrilizing (N)"), name = "Action Type")


# Create a plot
# p_iga_function_merged <- ggplot() +
#   geom_rect(aes(xmin = min(x_values), xmax = 1, ymin = 0.5, ymax = 1), fill = "gray", alpha = 0.1) +
#   geom_rect(aes(xmin = 1, xmax = max(x_values), ymin = 0.5, ymax = 1), fill = "black", alpha = 0.1) +
#   labs(title = "Rate of SIgA Function", x = "Affinity", y = "") +
#   geom_line(data = df_long_short_y1_high, aes(x = x, y = value, color = "Coating (C)"), alpha = 1, size = 1) +
#   geom_line(data = df_long_short_y2_high, aes(x = x, y = value, color = "Neutrilizing (N)"), linetype='dashed',alpha = 0.2, size = 1) +
#   geom_line(data = df_long_short_y1_low, aes(x = x, y = value, color = "Coating (C)"), linetype='dashed',alpha = 0.6, size = 1) +
#   geom_line(data = df_long_short_y2_low, aes(x = x, y = value, color = "Neutrilizing (N)"), alpha = 1, size = 1) +
#   geom_line(data = df_long_short_w,  aes(x = x, y = value, color = "Binding Ability"), linetype = '11', alpha = 1, size=1) +
#   scale_x_log10() + 
#   theme_minimal() +  
#   theme(plot.title = element_text(size = 12)) +
#   scale_color_manual(values = c("Binding Ability" = "#FF6333", "Coating (C)" = "gray", "Neutrilizing (N)" = "black"), name = "") +
#   guides(color = guide_legend(title = NULL, override.aes = list(linetype = c("11","solid", "solid"))))


p_iga_function_merged <- ggplot() +
  geom_rect(aes(xmin = min(x_values), xmax = 1, ymin = 0.5, ymax = 1), fill = "gray", alpha = 0.1) +
  geom_rect(aes(xmin = 1, xmax = max(x_values), ymin = 0.5, ymax = 1), fill = "black", alpha = 0.1) +
  labs(title = "Rate of SIgA Function", x = "Affinity", y = "") +
  geom_line(data = df_long_short_y1_high, aes(x = x, y = value, color = "Coating (C)"), alpha = 1, size = 1) +
  geom_line(data = df_long_short_y2_high, aes(x = x, y = value, color = "Neutrilizing (N)"), linetype='dashed',alpha = 0.2, size = 1) +
  geom_line(data = df_long_short_y1_low, aes(x = x, y = value, color = "Coating (C)"), linetype='dashed',alpha = 0.6, size = 1) +
  geom_line(data = df_long_short_y2_low, aes(x = x, y = value, color = "Neutrilizing (N)"), alpha = 1, size = 1) +
  scale_x_log10() + 
  theme_minimal() +  
  theme(plot.title = element_text(size = 12)) +
  scale_color_manual(values = c("Coating (C)" = "gray", "Neutrilizing (N)" = "black"), name = "") +
  guides(color = guide_legend(title = NULL, override.aes = list(linetype = c("solid", "solid"))))


out_df_ene_last$miga = keep_miga_vec
# Add points with colors defined in the point_color column directly
p_iga_function_merged = p_iga_function_merged +
  # geom_point(data = out_df_ene_last, aes(x = median, y = val), color=out_df_ene_last$point_color, size = 5) #eSIgA
  geom_point(data = out_df_ene_last, aes(x = miga, y = val), color=out_df_ene_last$point_color, size = 5) #mSIgA


# Note: Since we're directly using color outside of aes(), it bypasses ggplot2's scale system, so no scale_color_identity() is needed.

p_iga_function_merged = p_iga_function_merged+ theme(
  #legend.title = element_blank(),
  # legend.text = element_text(size = 10),
  # axis.title = element_text(size = 10),
  # axis.text = element_text(size = 10),
  legend.position = c(0.638, 0.016), # This needs adjustment
  legend.justification = c("right", "bottom"), # Anchor point
  legend.box.just = "right",
  legend.box.background = element_rect(colour = "gray", fill='white'),
  legend.title.align = 0.5,
  # legend.margin = margin(b = 0, t=-4), # Increased bottom margin for padding
  legend.title = element_text(size = 11), # Set the legend title size here
  legend.text = element_text(size = 10)   # Set the legend text size here
) 

annotate_points = c()
# annotate_points$x = out_df_ene_last$median
annotate_points$x = c(12,0.05,0.05,0.05)
annotate_points$y = c(0.56,0.66,0.76,0.56)
annotate_points=as.data.frame(annotate_points)

labels = c()
# annotate_points$label[1] = paste0(out_df_ene_last$taxa[1],', ',round(1-1/(1+out_df_ene_last$median[1]),2))
# annotate_points$label[2] = paste0(out_df_ene_last$taxa[2],', ',round(1/(1+out_df_ene_last$median[2]),2))
# annotate_points$label[3] = paste0(out_df_ene_last$taxa[3],', ',round(1/(1+out_df_ene_last$median[3]),2))
# annotate_points$label[4] = paste0(out_df_ene_last$taxa[4],', ',round(1/(1+out_df_ene_last$median[4]),2))

annotate_points$label[1] = paste0(out_df_ene_last$taxa[1],', ',round(1-1/(1+keep_miga_vec[1]),2))
annotate_points$label[2] = paste0(out_df_ene_last$taxa[2],', ',round(1/(1+keep_miga_vec[2]),2))
annotate_points$label[3] = paste0(out_df_ene_last$taxa[3],', ',round(1/(1+keep_miga_vec[3]),2))
annotate_points$label[4] = paste0(out_df_ene_last$taxa[4],', ',round(1/(1+keep_miga_vec[4]),2))


annotate_points$category[1] = 'E'
annotate_points$category[2] = 'B'
annotate_points$category[3] = 'BC'
annotate_points$category[4] = 'C'

p_iga_function_merged = p_iga_function_merged + 
  geom_text(data = subset(annotate_points, category == 'E'), aes(x = x, y = y, label = label), color = "#f20026",  vjust = "inward", hjust = "inward", fontface = "bold") +
  geom_text(data = subset(annotate_points, category == 'BC'), aes(x = x, y = y, label = label), color = '#22a9ba',  vjust = "inward", hjust = "inward", fontface = "bold") +
  geom_text(data = subset(annotate_points, category == 'B'), aes(x = x, y = y, label = label), color = '#074fa8',  vjust = "inward", hjust = "inward", fontface = "bold") +
  geom_text(data = subset(annotate_points, category == 'C'), aes(x = x, y = y, label = label), color = '#ff717b',  vjust = "inward", hjust = "inward", fontface = "bold") 


