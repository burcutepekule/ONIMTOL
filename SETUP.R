#################################
library(tidyverse)
library(cowplot)  # this is a plotting library but I'm unsure if you want to retain it
library(readxl)
library(icesTAF)
library(plyr)
library(dplyr)
library(writexl)
# library(cmdstanr)
library(loo)
library(mvtnorm)
library(MASS)
library(deSolve)
library(patchwork)

options(width = 10000)

options(xtable.floating = FALSE)
options(xtable.timestamp = "")

todaystr            = format(Sys.Date(), "%d%m%Y");
direc2save          = paste0("OUT/",todaystr,"/")
mkdir(direc2save)

theme_set(theme_bw())
mc.cores = parallel::detectCores()
print(c('cores: ',mc.cores))


# # Function to read a CSV and extract samples
# read_samples <- function(file, wu) {
#   # Read only the header to get column names
#   col_names <- names(read.csv(file, nrows = 0, comment.char = "#"))
#   
#   # Create a colClasses vector, defaulting to "NULL" (skips the column)
#   col_classes <- setNames(rep("NULL", length(col_names)), col_names)
#   
#   # Set the columns you want to keep to their appropriate type (e.g., "numeric")
#   col_classes[keep_cols] <- "numeric"
#   
#   # Read the CSV while skipping the initial metadata and comments, and only selected columns
#   data <- read.csv(file, comment.char = "#", skip = 10, colClasses = col_classes)
#   
#   # Optionally, filter out warmup samples if they are included
#   data <- data[((wu+1):nrow(data)), ]
#   
#   return(data)
# }
bin2months = function(df_in) {
  df_in = df_in %>% mutate(bin_idx = case_when(days <= 30 ~ 1, days <= 60 ~ 2, days <= 90 ~ 3, days <= 120 ~ 4, days <= 150 ~ 5, days <= 180 ~ 6, days <= 210 ~ 7, days <= 240 ~ 8,   days <= 270 ~ 9, days <= 300 ~ 10, days <= 330 ~ 11, days > 330 ~ 12 ))

  return(df_in)
}


bin2weeks = function(df_in) {
  df_in = df_in %>%mutate(bin_idx = case_when(
    days <= 7 ~ 1,
    days <= 14 ~ 2,
    days <= 21 ~ 3,
    days <= 28 ~ 4,
    days <= 35 ~ 5,
    days <= 42 ~ 6,
    days <= 49 ~ 7,
    days <= 56 ~ 8,
    days <= 63 ~ 9,
    days <= 70 ~ 10,
    days <= 77 ~ 11,
    days <= 84 ~ 12,
    days <= 91 ~ 13,
    days <= 98 ~ 14,
    days <= 105 ~ 15,
    days <= 112 ~ 16,
    days <= 119 ~ 17,
    days <= 126 ~ 18,
    days <= 133 ~ 19,
    days <= 140 ~ 20,
    days <= 147 ~ 21,
    days <= 154 ~ 22,
    days <= 161 ~ 23,
    days <= 168 ~ 24,
    days <= 175 ~ 25,
    days <= 182 ~ 26,
    days <= 189 ~ 27,
    days <= 196 ~ 28,
    days <= 203 ~ 29,
    days <= 210 ~ 30,
    days <= 217 ~ 31,
    days <= 224 ~ 32,
    days <= 231 ~ 33,
    days <= 238 ~ 34,
    days <= 245 ~ 35,
    days <= 252 ~ 36,
    days <= 259 ~ 37,
    days <= 266 ~ 38,
    days <= 273 ~ 39,
    days <= 280 ~ 40,
    days <= 287 ~ 41,
    days <= 294 ~ 42,
    days <= 301 ~ 43,
    days <= 308 ~ 44,
    days <= 315 ~ 45,
    days <= 322 ~ 46,
    days <= 329 ~ 47,
    days > 329 ~ 48
  ))
  return(df_in)
}

# Function to convert a number into a vector of its digits
number_to_vector <- function(number) {
  # Convert the number to a character string, then split it into individual characters
  digits_as_chars <- strsplit(as.character(number), "")
  # Convert the character vector to a numeric vector
  digits_as_numbers <- as.numeric(digits_as_chars[[1]])
  return(digits_as_numbers)
}

div0 <- function(x,y){
  if(x==0 && y==0){
    return(0)
  }else{
    return(x/y)
  }
  
}
read_from_data_sim <- function(colname, data, day_start){
  string             = data[37, 1]
  clean_string       = trimws(string)
  columnames         = strsplit(clean_string, ",")[[1]]
  inds               = grep(paste0("^",colname), columnames)
  
  string             = data[38, 1]
  clean_string       = trimws(string)
  numbers            = strsplit(clean_string, ",")[[1]]
  named_array        = as.numeric(numbers[inds])
  names(named_array) = columnames[inds]
  
  # Split the string into individual numbers
  
  # Extract the names and the values
  names_array = names(named_array)
  values_array= as.numeric(named_array) # Ensure the abundance values are numeric
  
  # Split the names to extract time and taxa indices
  split_names = strsplit(names_array, "\\.")
  
  # Extract time and taxa and convert them to the correct type
  times = sapply(split_names, function(x) as.numeric(x[2]))+day_start-1
  taxa  = sapply(split_names, function(x) as.numeric(x[3]))
  
  # Create the dataframe
  df_named <- data.frame(
    median = values_array,
    taxa = taxa_array[taxa],
    time = times
  )
  return(df_named)
}

read_from_data_sim_2 <- function(colname, data, data_list_sim){
  string             = data[37, 1]
  clean_string       = trimws(string)
  columnames         = strsplit(clean_string, ",")[[1]]
  inds               = grep(paste0("^",colname), columnames)
  
  string             = data[38, 1]
  clean_string       = trimws(string)
  numbers            = strsplit(clean_string, ",")[[1]]
  named_array        = as.numeric(numbers[inds])
  names(named_array) = columnames[inds]
  
  # Split the string into individual numbers
  
  # Extract the names and the values
  names_array = names(named_array)
  values_array= as.numeric(named_array) # Ensure the abundance values are numeric
  
  # Split the names to extract time and taxa indices
  split_names = strsplit(names_array, "\\.")
  
  # Extract time and taxa and convert them to the correct type
  times = data_list_sim$ts_pred
  taxa  = sapply(split_names, function(x) as.numeric(x[3]))
  
  # Create the dataframe
  df_named <- data.frame(
    median = values_array,
    taxa = taxa_array[taxa],
    time = times
  )
  return(df_named)
}

read_from_data_sim_6 <- function(colname, data, data_list_sim){
  string             = data[37, 1]
  clean_string       = trimws(string)
  columnames         = strsplit(clean_string, ",")[[1]]
  inds               = grep(paste0("^",colname), columnames)
  
  string             = data[38, 1]
  clean_string       = trimws(string)
  numbers            = strsplit(clean_string, ",")[[1]]
  named_array        = as.numeric(numbers[inds])
  names(named_array) = columnames[inds]
  
  # Split the string into individual numbers
  
  # Extract the names and the values
  names_array = names(named_array)
  values_array= as.numeric(named_array) # Ensure the abundance values are numeric
  
  # Split the names to extract time and taxa indices
  split_names = strsplit(names_array, "\\.")
  
  # Extract time and taxa and convert them to the correct type
  times = data_list_sim$t_data
  taxa  = sapply(split_names, function(x) as.numeric(x[3]))
  
  # Create the dataframe
  df_named <- data.frame(
    median = values_array,
    taxa = taxa_array[taxa],
    time = times
  )
  return(df_named)
}

read_from_data <- function(colname, data, day_start, i){
  string             = data[45, 1]
  clean_string       = trimws(string)
  columnames         = strsplit(clean_string, ",")[[1]]
  inds               = grep(paste0("^",colname), columnames)
  
  string             = data[i, 1]
  clean_string       = trimws(string)
  numbers            = strsplit(clean_string, ",")[[1]]
  named_array        = as.numeric(numbers[inds])
  names(named_array) = columnames[inds]
  
  # Split the string into individual numbers
  
  # Extract the names and the values
  names_array = names(named_array)
  values_array= as.numeric(named_array) # Ensure the abundance values are numeric
  
  # Split the names to extract time and taxa indices
  split_names = strsplit(names_array, "\\.")
  
  # Extract time and taxa and convert them to the correct type
  times = sapply(split_names, function(x) as.numeric(x[2]))+day_start-1
  taxa  = sapply(split_names, function(x) as.numeric(x[3]))
  
  # Create the dataframe
  df_named <- data.frame(
    median = values_array,
    taxa = taxa_array[taxa],
    time = times
  )
  return(df_named)
}

read_from_data_noT <- function(colname, data, day_start, i){
  string             = data[45, 1]
  clean_string       = trimws(string)
  columnames         = strsplit(clean_string, ",")[[1]]
  inds               = grep(paste0("^",colname), columnames)
  
  string             = data[i, 1]
  clean_string       = trimws(string)
  numbers            = strsplit(clean_string, ",")[[1]]
  named_array        = as.numeric(numbers[inds])
  names(named_array) = columnames[inds]
  
  # Split the string into individual numbers
  
  # Extract the names and the values
  names_array = names(named_array)
  values_array= as.numeric(named_array) # Ensure the abundance values are numeric
  
  # Split the names to extract time and taxa indices
  split_names = strsplit(names_array, "\\.")
  
  # Extract time and taxa and convert them to the correct type
  taxa  = sapply(split_names, function(x) as.numeric(x[2]))
  
  # Create the dataframe
  df_named <- data.frame(
    median = values_array,
    taxa = taxa_array[taxa]
  )
  return(df_named)
}

read_from_data_noT_sim <- function(colname, data, day_start){
  string             = data[37, 1]
  clean_string       = trimws(string)
  columnames         = strsplit(clean_string, ",")[[1]]
  inds               = grep(paste0("^",colname), columnames)
  
  string             = data[38, 1]
  clean_string       = trimws(string)
  numbers            = strsplit(clean_string, ",")[[1]]
  named_array        = as.numeric(numbers[inds])
  names(named_array) = columnames[inds]
  
  # Split the string into individual numbers
  
  # Extract the names and the values
  names_array = names(named_array)
  values_array= as.numeric(named_array) # Ensure the abundance values are numeric
  
  # Split the names to extract time and taxa indices
  split_names = strsplit(names_array, "\\.")
  
  # Extract time and taxa and convert them to the correct type
  taxa  = sapply(split_names, function(x) as.numeric(x[2]))
  
  # Create the dataframe
  df_named <- data.frame(
    median = values_array,
    taxa = taxa_array[taxa]
  )
  return(df_named)
}

read_from_data_noT_sim_2 <- function(colname, data, data_list_sim){
  string             = data[37, 1]
  clean_string       = trimws(string)
  columnames         = strsplit(clean_string, ",")[[1]]
  inds               = grep(paste0("^",colname), columnames)
  
  string             = data[38, 1]
  clean_string       = trimws(string)
  numbers            = strsplit(clean_string, ",")[[1]]
  named_array        = as.numeric(numbers[inds])
  names(named_array) = columnames[inds]
  
  # Split the string into individual numbers
  
  # Extract the names and the values
  names_array = names(named_array)
  values_array= as.numeric(named_array) # Ensure the abundance values are numeric
  
  # Split the names to extract time and taxa indices
  split_names = strsplit(names_array, "\\.")
  
  # Extract time and taxa and convert them to the correct type
  times = data_list_sim$ts_pred
  
  # Create the dataframe
  df_named <- data.frame(
    median = values_array,
    taxa = taxa_array[taxa]
  )
  return(df_named)
}

# # Function to read a CSV and extract samples
# read_samples <- function(file, wu) {
#   # Read the CSV while skipping the initial metadata and comments
#   data <- read.csv(file, comment.char = "#")
#   # Optionally, filter out warmup samples if they are included
#   data <- data[((wu+1):nrow(data)), ] # no footer if not ended, so until nrow(data)
#   return(data)
# }

# Function to read a CSV and extract samples
read_samples <- function(file, wu, iter, readall) {
  
  if(readall==1){
    # Explicit column names
    explicit_cols <- c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__",
                       "growthRate_vector.1", "growthRate_vector.2", "growthRate_vector.3", "growthRate_vector.4",
                       "interactionMat_vector_Agnostic.1", "interactionMat_vector_Agnostic.2", 
                       "interactionMat_vector_Agnostic.3", "interactionMat_vector_Agnostic.4",
                       "interactionMat_vector_Gnostic.1", "interactionMat_vector_Gnostic.2",
                       "interactionMat_vector_Gnostic.3", "interactionMat_vector_Gnostic.4", 
                       "interactionMat_vector_Gnostic.5", "interactionMat_vector_Gnostic.6", 
                       "interactionMat_vector_Gnostic.7", "interactionMat_vector_Gnostic.8", 
                       "interactionMat_vector_Gnostic.9", "interactionMat_vector_Gnostic.10", 
                       "interactionMat_vector_Gnostic.11", "interactionMat_vector_Gnostic.12",
                       "mIgA_reactivity_vector_ess", "mIgA_reactivity_vector_lbblcr.1", 
                       "mIgA_reactivity_vector_lbblcr.2", "mIgA_reactivity_vector_lbblcr.3",
                       "O2_degredation", "srate_c", "srate_uc", "srate_k", "kappa", 
                       "add_with_food_bc", "add_with_food_c","target_maternal_add_summed",
                       "target_converged_add_summed","target_day_40","target_iga")
    
    # Patterns for column names
    patterns <- c('^output_mat_total_feces_daily_maternal\\..*$', '^output_mat_total_feces_daily_converged\\..*$', 
                  '^output_mat_coating_ratio_maternal\\..*$','^output_mat_total_lumen_daily_converged\\..*$',
                  '^output_mat_coating_ratio_converged\\..*$', '^output_mat_solid_converged\\..*$',
                  '^output_endogenous_coating_daily\\..*$','^output_endogenous_killing_daily\\..*$',
                  '^output_mat_coated_lumen_converged\\..*$','^output_mat_uncoated_lumen_converged\\..*$',
                  '^output_mat_total_feces_daily_maternal_bm0\\..*$',
                  '^output_mat_total_feces_daily_maternal_bm30\\..*$',
                  '^output_mat_total_feces_daily_converged\\..*$', 
                  '^output_mat_O2_maternal\\..*$', '^output_mat_HMO_maternal\\..*$', 
                  '^output_mat_mIgA_maternal\\..*$', '^output_maternal_coating_daily\\..*$',
                  '^output_maternal_killing_daily\\..*$', '^theta_out\\.[0-9]+$',
                  '^iga_index_maternal\\..*$', '^iga_index_converged\\..*$')
    
  }else{
    # Explicit column names
    explicit_cols <- c("lp__","target_maternal_add_summed","target_converged_add_summed","target_day_40","target_iga")
    
    # Patterns for column names
    patterns <- c('^output_mat_total_feces_daily_maternal\\..*$','^output_mat_total_feces_daily_converged\\..*$',
                  '^output_mat_total_feces_daily_maternal_bm0\\..*$','^output_mat_total_feces_daily_maternal_bm30\\..*$',
                  '^iga_index_maternal\\..*$', '^iga_index_converged\\..*$')
    
  }
  
  col_names <- names(fread(file, nrows = 0))  # Adjust the skip according to your file's structure
  cols_to_read <- col_names %in% explicit_cols
  for (pattern in patterns) {
    cols_to_read <- cols_to_read | grepl(pattern, col_names)
  }
  cols_to_read = which(cols_to_read==TRUE)
  # cols_to_read <- col_names[cols_to_read]  # Keep only the names of the columns to read
  n_lines     = length(readLines(file))
  start_row_1 = 47 # first not-commented line
  if(n_lines>(start_row_1-1+wu)){
    nrows_1     = wu
    start_row_2 = start_row_1-2+nrows_1+5 # commented lines in file
    nrows_2     = iter-wu
    if((nrows_2+start_row_2)>n_lines){
      nrows_2 = n_lines-start_row_2
    }
    data_read_wu_samples   = fread(file, skip = start_row_1-2, nrows = nrows_1, select = cols_to_read) #-2 BECAUSE HEADER
    data_read_iter_samples = fread(file, skip = start_row_2, nrows = nrows_2, select = cols_to_read, data.table = FALSE)
    colnames(data_read_iter_samples) = colnames(data_read_wu_samples)
    return(data_read_iter_samples)
  }else{
    data_read_wu_samples   = fread(file, skip = start_row_1-2, nrows = n_lines, select = cols_to_read) #-2 BECAUSE HEADER
    print('Only WU samples returned')
    return(data_read_wu_samples)
  }
  
}


### Compute summary statistics
### Function to compute summary statistics for a vector
compute_summary <- function(x) {
  data.frame(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q2.5 = quantile(x, probs = 0.025, na.rm = TRUE),
    q5 = quantile(x, probs = 0.050, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q95 = quantile(x, probs = 0.95, na.rm = TRUE),
    q97.5 = quantile(x, probs = 0.975, na.rm = TRUE)
  )
}

qsum = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))

logit = function(x) log(x/(1-x))

inv.logit = function(x) exp(x)/(1+exp(x))

plot_with_threshold <- function(input_array, threshold, title) {
  
  if(length(threshold)==1){
    input_array = input_array[is.finite(input_array)]
    # Determine the range for the y-axis
    y_range <- range(input_array, threshold)
    
    # Plot the array with black points and adjust y-axis limits
    plot(input_array, type = 'p', col = 'black', pch = 16, ylim = y_range, main = title)
    
    # Add a horizontal red line at the threshold value
    abline(h = threshold, col = 'red', lwd = 2)
  }else{
    input_array = input_array[is.finite(input_array)]
    threshold   = threshold[is.finite(threshold)]
    
    # Determine the range for the y-axis
    y_range <- range(input_array, threshold)
    
    # Plot the array with black points and adjust y-axis limits
    plot(input_array, type = 'p', col = 'black', pch = 16, ylim = y_range, main = title)
    
    # Add a horizontal red line at the threshold value
    lines(threshold, col = 'red', lwd = 2)
  }
  
}






