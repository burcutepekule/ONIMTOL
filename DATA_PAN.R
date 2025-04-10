####### VALUES EXTRACTED BY CLAUDE USING AN ENHANCED VERSION OF THE FIGURE 2A IN PAN ET AL. ####### 

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggtext)
library(readxl)

# Create the data frame
microbiota_data <- data.frame(
  Genus = rep(c("Escherichia-Shigella", "Bifidobacterium", "Clostridium", 
                "Lactobacillus", "Bacteroides", "Other"), each = 3),
  FeedingGroup = rep(c("Breastfed", "Partially breastfed", "Formula feeding"), times = 6),
  Proportion = c(12, 22, 32,  # Escherichia-Shigella
                 47, 32, 21,  # Bifidobacterium
                 6, 8, 12,    # Clostridium
                 20, 18, 7,   # Lactobacillus
                 8, 9, 23,    # Bacteroides
                 7, 6, 4)     # Other
)

# Create a data frame for p-values
p_values <- data.frame(
  Genus = c("Escherichia-Shigella", "Bifidobacterium", "Clostridium", 
            "Lactobacillus", "Bacteroides", "Other"),
  P1 = c("P1≈0.0015", "P1≈0.0374", "P1≈0.0150", "P1≈0.0101", "P1≈0.0308", "P1≈0.0415"),
  P2 = c("P2<0.0001", "P2≈0.0008", "P2≈0.0081", "P2<0.0001", "P2≈0.0006", "P2≈0.0201")
)

# Convert Genus to a factor with levels in reverse order (for bottom-to-top plotting)
microbiota_data$Genus <- factor(microbiota_data$Genus, 
                                levels = rev(c("Escherichia-Shigella", "Bifidobacterium", "Clostridium", 
                                               "Lactobacillus", "Bacteroides", "Other")))

# Convert feeding group to a factor to ensure consistent colors
microbiota_data$FeedingGroup <- factor(microbiota_data$FeedingGroup, 
                                       levels = c("Breastfed", "Partially breastfed", "Formula feeding"))

# Create the plot
p <- ggplot(microbiota_data, aes(x = Proportion, y = Genus, fill = FeedingGroup)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(aes(label = paste0(Proportion, "%")), 
            position = position_dodge(width = 0.7), 
            hjust = -0.2, vjust = 0.5, size = 3.5) +
  scale_fill_manual(values = c("Breastfed" = "#FF0000", 
                               "Partially breastfed" = "#0066FF", 
                               "Formula feeding" = "#00CC66")) +
  scale_x_continuous(limits = c(0, 55), breaks = seq(0, 50, by = 5)) +
  labs(title = "Wilcoxon rank-sum test bar plot on Genus level (30-42 days)",
       x = "Proportion (%)",
       y = NULL,
       fill = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10)
  )

# Add p-values to the right side of the plot
# First, convert p_values to match the factor levels
p_values$Genus <- factor(p_values$Genus, levels = levels(microbiota_data$Genus))

# Create a function to add p-values as annotations
add_p_values <- function(plot, p_values_df) {
  max_x <- 55
  
  for (i in 1:nrow(p_values_df)) {
    genus <- as.character(p_values_df$Genus[i])
    p1 <- p_values_df$P1[i]
    p2 <- p_values_df$P2[i]
    
    y_pos <- which(levels(microbiota_data$Genus) == genus)
    
    plot <- plot + 
      annotate("text", x = max_x-4.5, y = y_pos + 0.1, label = p1, hjust = 0, size = 3.5) +
      annotate("text", x = max_x-4.5, y = y_pos - 0.1, label = p2, hjust = 0, size = 3.5)
  }
  
  return(plot)
}

# Add p-values to the plot
final_plot <- add_p_values(p, p_values)

# Display the plot
print(final_plot)

# Save the plot as a PNG file
ggsave("./PAN_DATA/Figure2a_quantified.png", final_plot, width = 12, height = 8, dpi = 300, bg = "white")

############ NOW MATCH THESE VALUES WITH RELATIVE ABUNDANCES IN TSUKUDA ET AL. ###############

# RELATIVE ABUNDANCES AT DOL 30
df_relative_smooth    = read_excel("./TSUKUDA_DATA/RELATIVE_ABUNDANCE_DATA_SMOOTHED.xlsx")
df_relative_smooth_30 = df_relative_smooth %>% filter(day==30)
map_30 = sum(df_relative_smooth_30$smoothed_values)

# MATCH 
microbiota_data_formula = microbiota_data %>% filter(FeedingGroup=="Formula feeding")
microbiota_data_formula = microbiota_data_formula %>% dplyr::filter(Genus %in% c('Escherichia-Shigella', 'Bifidobacterium', 'Clostridium', 'Bacteroides'))
sum_proportion          = sum(microbiota_data_formula$Proportion)
microbiota_data_formula = microbiota_data_formula %>% dplyr::rowwise() %>% dplyr::mutate(`Relative Abundance` = round((1/100)*Proportion*map_30/(sum_proportion),2))

# FILTER THE ONES IN THE MODEL
microbiota_data_formula = microbiota_data_formula[c('Genus','Relative Abundance')]

# PUT BACK THE DAY AND CONDITION
microbiota_data_formula$day = 30
microbiota_data_formula$FeedingGroup = "Formula feeding"

write_xlsx(microbiota_data_formula, "./PAN_DATA/RELATIVE_ABUNDANCE_DATA_FORMULA.xlsx")



