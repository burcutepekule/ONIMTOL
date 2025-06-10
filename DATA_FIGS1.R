# Load necessary libraries
rm(list=ls())
library("bayesplot")
library("ggplot2")
library("rstanarm") 
library('rstan')
library('brms')
library('dplyr')
library(readxl)
library(writexl)
library(deSolve)
library('openxlsx')
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library('data.table')
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

## LOAD MERGED SAMPLES AND SUMMARY FILES OF THE MCMC CHAINS
combined_samples = readRDS('./combined_samples_chains.rds')
summary_df       = readRDS('./summary_df_chains.rds')

# Function to count divergent transitions in combined samples
count_divergent_transitions <- function(combined_samples) {
  # Check if 'divergent__' column exists in the data
  if("divergent__" %in% names(combined_samples)) {
    divergent_count <- sum(combined_samples$divergent__ == 1)
    total_transitions <- nrow(combined_samples)
    divergent_percent <- round(divergent_count / total_transitions * 100, 1)
    
    cat(paste0(divergent_count, " of ", total_transitions, " (", 
               divergent_percent, "%) transitions ended with a divergence.\n"))
    
    return(list(
      divergent_count = divergent_count,
      total_transitions = total_transitions,
      divergent_percent = divergent_percent
    ))
  } else {
    # Try to find diagnostic information in other formats
    # For cmdstanr, the column might be named differently
    potential_columns <- grep("divergent|accept_stat__", names(combined_samples), value = TRUE)
    
    if(length(potential_columns) > 0) {
      cat("Potential diagnostic columns found: ", paste(potential_columns, collapse=", "), "\n")
      # You may need to adjust this based on your specific output format
    } else {
      cat("No divergence information found in the combined samples.\n")
      cat("Diagnostic information might not have been included in the saved samples.\n")
    }
  }
}

# Now call this function for each of your AG_MERGED chain sets
# After reading in your combined samples:

# Check for divergent transitions in the combined chains
cat("\nDiagnostic summary for AG_MERGED chains:\n")
divergence_stats <- count_divergent_transitions(combined_samples)


my_colors = c('#074fa8','#2bc5d8','#ff9aa1','#f20026')
abundanceArray_meanSubjects_maternal = readRDS('./abundanceArray_meanSubjects_maternal.rds')
taxa_array = colnames(abundanceArray_meanSubjects_maternal)
numTaxa    = length(taxa_array)
abundanceArray_meanSubjects_maternal <- abundanceArray_meanSubjects_maternal %>%
  dplyr::mutate(across(c(Enterobacteriaceae, Bifidobacteriaceae, Bacteroidaceae, Clostridiales), 
                       ~pmax(., 0, na.rm = TRUE))) 

# Convert to longer format
abundanceArray_meanSubjects_use     = abundanceArray_meanSubjects_maternal
abundanceArray_meanSubjects_use$day = rownames(abundanceArray_meanSubjects_use)
abundanceArray_meanSubjects_longer  = abundanceArray_meanSubjects_use %>% pivot_longer(!day,names_to = 'taxa', values_to = 'abundance')
abundanceArray_meanSubjects_longer$day=as.numeric(abundanceArray_meanSubjects_longer$day)

################################     EXTACT FROM RSTAN OUTPUT  ###################################################
##################################################################################################################
pattern_maternal <- '^output_mat_total_feces_daily_maternal\\..*$'
summaryTable_maternal <- summary_df[grep(pattern_maternal, summary_df$variable), ]
days_array_sparse <- seq(from = 1, to = dim(summaryTable_maternal)[1]/numTaxa, 1)

# Extracting the time index (i) and taxa index (j) from the 'variable' column
summaryTable_maternal <- summaryTable_maternal %>%
  mutate(
    time = days_array_sparse[as.numeric(gsub("output_mat_total_feces_daily_maternal\\.(\\d+)\\..*", "\\1", variable))],
    taxa = as.numeric(gsub("output_mat_total_feces_daily_maternal\\.\\d+\\.(\\d+)", "\\1", variable))
  ) %>%
  mutate(
    taxa = ifelse(!is.na(taxa), taxa_array[taxa], NA)
  )

# ADD IC
summaryTable_maternal$time = summaryTable_maternal$time+1

# Count taxa and get end day
numTaxa = length(taxa_array)
t_end   = max(as.numeric(unique(abundanceArray_meanSubjects_longer$day)))

graphics.off()

summaryTable_maternal_use = summaryTable_maternal[c('time','taxa','median')]
colnames(summaryTable_maternal_use)[1] ='day'
colnames(summaryTable_maternal_use)[3] ='abundance'
summaryTable_maternal_use$type = 'Model'
abundanceArray_meanSubjects_longer$type = 'Observations'
AG_MERGED_plot_data = rbind(summaryTable_maternal_use,abundanceArray_meanSubjects_longer)

p = ggplot() +
  geom_ribbon(data = summaryTable_maternal, aes(x = time, ymin = q5, ymax = q95, fill = taxa), alpha = .2) +
  geom_line(data = AG_MERGED_plot_data, aes(x = day, y = abundance, color=taxa, linetype = type), size = 0.8) +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=='Enterobacteriaceae'),aes(x = day, y = abundance),  color="#910017",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Bifidobacteriaceae"),aes(x = day, y = abundance),  color="#042f65",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Bacteroidaceae"),aes(x = day, y = abundance),  color="#187783",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Clostridiales"),aes(x = day, y = abundance),  color="#f50011",size = 0.8, linetype = "11") +
  scale_fill_manual(values = c("Enterobacteriaceae" = '#f20026', 
                               "Bifidobacteriaceae" = "#074fa8", 
                               "Bacteroidaceae" = '#2bc5d8',
                               "Clostridiales" ='#ff9aa1',
                               "Observations" = "black")) +
  scale_colour_manual(values = c("Enterobacteriaceae" = '#f20026', 
                                 "Bifidobacteriaceae" = "#074fa8", 
                                 "Bacteroidaceae" = '#2bc5d8',
                                 "Clostridiales" ='#ff9aa1')) +
  scale_linetype_manual(values = c('Model' = "solid",  # Correctly mapping the aesthetic to the line type
                                   'Observations' = "11")) + # Ensure this matches your data's linetype values
  facet_wrap(~ taxa, scales = "free", nrow = 2) +
  theme(strip.text = element_blank()) +
  labs(x = "Days", y = "Relative Abundance in Fecal Samples",
       fill = NULL, # Remove fill legend title
       color = NULL, # Remove color legend title
       linetype = NULL) # Correctly label the legend for linetypes


pattern_maternal <- '^output_mat_total_feces_daily_maternal_bm0\\..*$'
summaryTable_maternal <- summary_df[grep(pattern_maternal, summary_df$variable), ]

# days_array_sparse <- seq(from = 1, to = 735, length.out = dim(summaryTable_maternal)[1]/numTaxa)
days_array_sparse <- seq(from = 1, to = dim(summaryTable_maternal)[1]/numTaxa, 1)

# Extracting the time index (i) and taxa index (j) from the 'variable' column
summaryTable_maternal <- summaryTable_maternal %>%
  mutate(
    time = days_array_sparse[as.numeric(gsub("output_mat_total_feces_daily_maternal_bm0\\.(\\d+)\\..*", "\\1", variable))],
    taxa = as.numeric(gsub("output_mat_total_feces_daily_maternal_bm0\\.\\d+\\.(\\d+)", "\\1", variable))
  ) %>%
  mutate(
    taxa = ifelse(!is.na(taxa), taxa_array[taxa], NA)
  )

# ADD IC
summaryTable_maternal$time = summaryTable_maternal$time+1

# Count taxa and get end day
numTaxa = length(taxa_array)
t_end   = max(as.numeric(unique(abundanceArray_meanSubjects_longer$day)))

summaryTable_maternal_use = summaryTable_maternal[c('time','taxa','median')]
colnames(summaryTable_maternal_use)[1] ='day'
colnames(summaryTable_maternal_use)[3] ='abundance'
summaryTable_maternal_use$type = 'Model'
abundanceArray_meanSubjects_longer$type = 'Observations'
AG_MERGED_plot_data = rbind(summaryTable_maternal_use,abundanceArray_meanSubjects_longer)

p_BM0 = ggplot() +
  geom_ribbon(data = summaryTable_maternal, aes(x = time, ymin = q5, ymax = q95, fill = taxa), alpha = .2) +
  geom_line(data = AG_MERGED_plot_data, aes(x = day, y = abundance, color=taxa, linetype = type), size = 0.8) +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=='Enterobacteriaceae'),aes(x = day, y = abundance),  color="#910017",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Bifidobacteriaceae"),aes(x = day, y = abundance),  color="#042f65",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Bacteroidaceae"),aes(x = day, y = abundance),  color="#187783",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Clostridiales"),aes(x = day, y = abundance),  color="#f50011",size = 0.8, linetype = "11") +
  scale_fill_manual(values = c("Enterobacteriaceae" = '#f20026', 
                               "Bifidobacteriaceae" = "#074fa8", 
                               "Bacteroidaceae" = '#2bc5d8',
                               "Clostridiales" ='#ff9aa1',
                               "Observations" = "black")) +
  scale_colour_manual(values = c("Enterobacteriaceae" = '#f20026', 
                                 "Bifidobacteriaceae" = "#074fa8", 
                                 "Bacteroidaceae" = '#2bc5d8',
                                 "Clostridiales" ='#ff9aa1')) +
  scale_linetype_manual(values = c('Model' = "solid",  # Correctly mapping the aesthetic to the line type
                                   'Observations' = "11")) + # Ensure this matches your data's linetype values
  facet_wrap(~ taxa, scales = "free", nrow = 2) +
  theme(strip.text = element_blank()) +
  labs(x = "Days", y = "Relative Abundance in Fecal Samples",
       fill = NULL, # Remove fill legend title
       color = NULL, # Remove color legend title
       linetype = NULL) # Correctly label the legend for linetypes


######################################################################################################################
######################################################################################################################


abundanceArray_meanSubjects_converged = readRDS('abundanceArray_meanSubjects_converged.rds')
colnames(abundanceArray_meanSubjects_converged)=taxa_array
abundanceArray_meanSubjects_converged <- abundanceArray_meanSubjects_converged %>%
  dplyr::mutate(across(c(Enterobacteriaceae, Bifidobacteriaceae, Bacteroidaceae, Clostridiales), 
                       ~pmax(., 0, na.rm = TRUE))) 

# Convert to longer format
abundanceArray_meanSubjects_use     = abundanceArray_meanSubjects_converged
abundanceArray_meanSubjects_use$day = rownames(abundanceArray_meanSubjects_use)
abundanceArray_meanSubjects_longer  = abundanceArray_meanSubjects_use %>% pivot_longer(!day,names_to = 'taxa', values_to = 'abundance')
abundanceArray_meanSubjects_longer$day=as.numeric(abundanceArray_meanSubjects_longer$day)


# Count taxa and get end day
numTaxa = length(taxa_array)
t_end   = max(as.numeric(unique(abundanceArray_meanSubjects_longer$day)))

pattern_maternal <- '^output_mat_total_feces_daily_converged\\..*$'
summaryTable_maternal <- summary_df[grep(pattern_maternal, summary_df$variable), ]

# days_array_sparse <- seq(from = 1, to = 735, length.out = dim(summaryTable_maternal)[1]/numTaxa)
days_array_sparse <- seq(from = 1, to = dim(summaryTable_maternal)[1]/numTaxa, 1)

# Extracting the time index (i) and taxa index (j) from the 'variable' column
summaryTable_maternal <- summaryTable_maternal %>%
  mutate(
    time = days_array_sparse[as.numeric(gsub("output_mat_total_feces_daily_converged\\.(\\d+)\\..*", "\\1", variable))],
    taxa = as.numeric(gsub("output_mat_total_feces_daily_converged\\.\\d+\\.(\\d+)", "\\1", variable))
  ) %>%
  mutate(
    taxa = ifelse(!is.na(taxa), taxa_array[taxa], NA)
  )

# ADD IC
summaryTable_maternal$time = summaryTable_maternal$time+1

# Count taxa and get end day
numTaxa = length(taxa_array)
t_end   = max(as.numeric(unique(abundanceArray_meanSubjects_longer$day)))

summaryTable_maternal_use = summaryTable_maternal[c('time','taxa','median')]
colnames(summaryTable_maternal_use)[1] ='day'
colnames(summaryTable_maternal_use)[3] ='abundance'
summaryTable_maternal_use$type = 'Model'
abundanceArray_meanSubjects_longer$type = 'Observations'
AG_MERGED_plot_data = rbind(summaryTable_maternal_use,abundanceArray_meanSubjects_longer)

p_conv = ggplot() +
  geom_ribbon(data = summaryTable_maternal, aes(x = time, ymin = q5, ymax = q95, fill = taxa), alpha = .2) +
  geom_line(data = AG_MERGED_plot_data, aes(x = day, y = abundance, color=taxa, linetype = type), size = 0.8) +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=='Enterobacteriaceae'),aes(x = day, y = abundance),  color="#910017",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Bifidobacteriaceae"),aes(x = day, y = abundance),  color="#042f65",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Bacteroidaceae"),aes(x = day, y = abundance),  color="#187783",size = 0.8, linetype = "11") +
  geom_line(data = subset(AG_MERGED_plot_data, type=='Observations' & taxa=="Clostridiales"),aes(x = day, y = abundance),  color="#f50011",size = 0.8, linetype = "11") +
  scale_fill_manual(values = c("Enterobacteriaceae" = '#f20026', 
                               "Bifidobacteriaceae" = "#074fa8", 
                               "Bacteroidaceae" = '#2bc5d8',
                               "Clostridiales" ='#ff9aa1',
                               "Observations" = "black")) +
  scale_colour_manual(values = c("Enterobacteriaceae" = '#f20026', 
                                 "Bifidobacteriaceae" = "#074fa8", 
                                 "Bacteroidaceae" = '#2bc5d8',
                                 "Clostridiales" ='#ff9aa1')) +
  scale_linetype_manual(values = c('Model' = "solid",  # Correctly mapping the aesthetic to the line type
                                   'Observations' = "11")) + # Ensure this matches your data's linetype values
  facet_wrap(~ taxa, scales = "free", nrow = 2) +
  theme(strip.text = element_blank()) +
  labs(x = "Days", y = "Relative Abundance in Fecal Samples",
       fill = NULL, # Remove fill legend title
       color = NULL, # Remove color legend title
       linetype = NULL) # Correctly label the legend for linetypes

spacer_grob = grid::nullGrob()

row_1 = cowplot::plot_grid(p, align='h', rel_heights = c(1),labels=c('A'),label_fontface = "bold")
row_2 = cowplot::plot_grid(p_BM0, align='h', rel_widths = c(1),labels=c('B'),label_fontface = "bold")
row_3 = cowplot::plot_grid(p_conv, align='h', rel_widths = c(1),labels=c('C'),label_fontface = "bold")
nested = cowplot::plot_grid(row_1,row_2,row_3,ncol = 1,align='v')

graphics.off()
png(file =paste0("./FIGURE_S1.png"),   # The directory you want to save the file in
    width     = 8,
    height    = 11,
    units     = "in",
    res       = 300)
nested
dev.off()
