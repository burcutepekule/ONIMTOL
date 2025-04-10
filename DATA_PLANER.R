# Setup
setwd(getwd())
library(zoo)
library(Rcpp)
library(lubridate)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(prodlim)
library(matrixStats)
library(dplyr)
library(plyr)
library(readxl)
library(writexl)
source('SETUP.R')

milk_th_exclusive = 0.75;
milk_th_mixed     = 0.5;

planer_data_iga = read_excel('PLANER_DATA/41586_2016_BFnature17940_MOESM321_ESM_edited.xlsx', sheet='SuppTable_13')
planer_data_iga$`Study ID` = round(as.numeric(planer_data_iga$`Study ID`),1)
write_xlsx(planer_data_iga, "PLANER_DATA/planer_data_iga.xlsx")

planer_data_meta = read_excel('PLANER_DATA/41586_2016_BFnature17940_MOESM321_ESM_edited.xlsx', sheet='SuppTable_2')
# I think column names are mixed here, so gonna fix to make it compatible with the prev. one
colnames_keep = colnames(planer_data_meta)[1:2]
colnames(planer_data_meta)[2:1]=colnames_keep
planer_data_iga$`Study ID` = round(as.numeric(planer_data_iga$`Study ID`),1)
write_xlsx(planer_data_meta, "PLANER_DATA/planer_data_meta.xlsx")

# Extract days and months from the second column
planer_data_iga <- planer_data_iga %>%
  mutate(
    days = as.numeric(str_extract(`Chronologic age of donor at time of sampling, days (months)`, "^[^(]+")),
    months = as.numeric(str_extract(`Chronologic age of donor at time of sampling, days (months)`, "(?<=\\()[^)]+"))
  ) %>% dplyr::select(-`Chronologic age of donor at time of sampling, days (months)`)

# Reorder columns to place "days" and "months" as the 2nd and 3rd columns
planer_data_iga <- planer_data_iga %>%
  dplyr::select(`Study ID`, days, months, everything())

# # Extract the bacteria names inside brackets and use them as new column names for columns 4 and onwards
# new_col_names <- sapply(names(planer_data_iga)[4:ncol(planer_data_iga)], function(col_name) {
#   str_extract(col_name, "(?<=\\().+[^)]")
# })
# names(planer_data_iga)[4:ncol(planer_data_iga)] <- new_col_names

# also change the column name for sampling  "Age at time of fecal sample (days)" to "days"
colnames(planer_data_meta)[3] = 'days'
colnames(planer_data_iga)[1]  = 'ID'
colnames(planer_data_meta)[2] = 'ID'

# Use SuppTable_5 for classification (https://jasbsci.biomedcentral.com/articles/10.1186/s40104-019-0402-1 for clostridium)
colnames_family    = colnames(planer_data_iga)
colnames_family[4] = 'Verrucomicrobiaceae' #"588471 (Akkermansia muciniphila)", NOT from SuppTable_5
colnames_family[5] = 'Ruminococcaceae'     #"C.6 (Ruminococcus torques)",       NOT from SuppTable_5
colnames_family[6] = 'Verrucomicrobiaceae' #"4306262 (Akkermansia muciniphila)",NOT from SuppTable_5
colnames_family[7] = 'Lachnospiraceae'     #"4436046 (Clostridium nexile)",     NOT from SuppTable_5
colnames_family[8] = 'Clostridiaceae'      #"4436046 (Clostridium spiroforme)", NOT from SuppTable_5, also not sure
colnames_family[9] = 'Bifidobacteriaceae'  #"365385 (Bifidobacterium bifidum)", NOT from SuppTable_5
colnames_family[10] = 'Ruminococcaceae'    #"C.4 (Ruminococcus gnavus)",        NOT from SuppTable_5
colnames_family[11] = 'Lachnospiraceae'    #"2555599 (Clostridium nexile)",     NOT from SuppTable_5
colnames_family[12] = 'Clostridiaceae'     #"4372528 (Clostridium glycolicum)", NOT from SuppTable_5
colnames_family[13] = 'Clostridiaceae'     #"4434334 (Clostridium bartlettii)", NOT from SuppTable_5, also not sure
colnames_family[14] = 'Enterobacteriaceae' #"C.3 (Escherichia coli)",           from SuppTable_5
colnames_family[15] = 'Lachnospiraceae'    #"539581 (Eubacterium cylindroides)",NOT from SuppTable_5
colnames_family[16] = 'Veillonellaceae'    #"C.40 (Other Veillonella)",         from SuppTable_5
colnames_family[17] = 'Lachnospiraceae'    #"C.8 (Other Roseburia)",            from SuppTable_5
colnames_family[18] = 'Lachnospiraceae'    #"1751298 (Roseburia intestinalis)", from SuppTable_5
colnames_family[19] = 'Bifidobacteriaceae' #"C.1 (Bifidobacterium longum)",     from SuppTable_5
colnames_family[20] = 'Streptococcaceae'   #"C.44 (Streptococcus thermophilus)",from SuppTable_5
colnames_family[21] = 'Bacteroidaceae'     # "C.15 (Bacteroides vulgatus)",     from SuppTable_5
colnames_family[22] = 'Clostridiaceae'     #"2943548 (Other Clostridiales)",    NOT from SuppTable_5, also not sure
colnames_family[23] = 'Clostridiaceae'     #"606927 (Peptoclostridium difficile)", NOT from SuppTable_5
colnames_family[24] = 'Lachnospiraceae'    #"4469576 (Clostridium bolteae)",       NOT from SuppTable_5, also not sure
colnames_family[25] = 'Clostridiaceae'     #"C.26 (Clostridium clostridioforme)",  from SuppTable_5, also not sure
colnames_family[26] = 'Ruminococcaceae'    #"4396688 (Ruminococcus torques)",      NOT SuppTable_5
colnames_family[27] = 'Clostridiaceae'     # "3709990 (Clostridium innocuum)",     NOT SuppTable_5
colnames_family[28] = 'Clostridiaceae'     # "195937 (Other Clostridiales)",       NOT SuppTable_5
colnames_family[29] = 'Lachnospiraceae'    # "C.12 (Blautia faecis)",              from SuppTable_5
colnames_family[30] = 'Clostridiaceae'     # "4453304 (Other Clostridiales)"       NOT from SuppTable_5
colnames_family[31] = 'Ruminococcaceae'    # "170462 (Ruminococcus sp id8)"        from SuppTable_5
colnames_family[32] = 'Clostridiaceae'     # "C.25 (Clostridium sp ss2 1)"         from SuppTable_5
colnames_family[33] = 'Ruminococcaceae'    # "C.39 (Ruminococcus sp ce2)"          from SuppTable_5
colnames(planer_data_iga)=colnames_family

planer_data_iga_longer = planer_data_iga %>% pivot_longer(!c('ID','days','months'), names_to = "family", values_to = "log_ratio")
planer_data_iga_longer = planer_data_iga_longer %>% filter(log_ratio != "--")
planer_data_iga_longer$log_ratio=as.numeric(planer_data_iga_longer$log_ratio)


planer_data_meta_short             = planer_data_meta[c('ID','days','Ratio of Breast milk : Formula')]
planer_data_iga_longer_meta_merged = merge(planer_data_iga_longer,planer_data_meta_short, by=c('ID','days'))
planer_data_iga_longer_meta_merged_milk_test   = planer_data_iga_longer_meta_merged %>% filter(days<154)
planer_data_iga_longer_meta_merged_milk_test_agg = aggregate(`Ratio of Breast milk : Formula` ~ ID, planer_data_iga_longer_meta_merged_milk_test, FUN=mean)
planer_data_iga_longer_meta_merged_milk_test      = planer_data_iga_longer_meta_merged_milk_test_agg %>% filter(`Ratio of Breast milk : Formula`>=milk_th_exclusive)

planer_data_iga_longer_meta_merged_milk_test_2     = planer_data_iga_longer_meta_merged %>% filter(days>154)
planer_data_iga_longer_meta_merged_milk_test_agg_2 = aggregate(`Ratio of Breast milk : Formula` ~ ID, planer_data_iga_longer_meta_merged_milk_test_2, FUN=mean)
planer_data_iga_longer_meta_merged_milk_test_2     = planer_data_iga_longer_meta_merged_milk_test_agg_2 %>% filter(`Ratio of Breast milk : Formula`<milk_th_mixed)



ids_keep = intersect(unique(planer_data_iga_longer_meta_merged_milk_test$ID),unique(planer_data_iga_longer_meta_merged_milk_test_2$ID))

planer_data_iga_longer_meta_merged_milk   = planer_data_iga_longer_meta_merged %>% filter(ID %in% ids_keep)

# no filtering
# planer_data_iga_longer_meta_merged_milk = planer_data_iga_longer_meta_merged
  
planer_data_iga_longer_ento = planer_data_iga_longer_meta_merged_milk %>% filter(family=='Enterobacteriaceae')
planer_data_iga_longer_bifi = planer_data_iga_longer_meta_merged_milk %>% filter(family=='Bifidobacteriaceae')
planer_data_iga_longer_bact = planer_data_iga_longer_meta_merged_milk %>% filter(family=='Bacteroidaceae')
planer_data_iga_longer_clos = planer_data_iga_longer_meta_merged_milk %>% filter(family=='Clostridiaceae' | family=='Ruminococcaceae' | family == 'Lachnospiraceae')

planer_data_iga_longer_ento = bin2months(planer_data_iga_longer_ento)
planer_data_iga_longer_bifi = bin2months(planer_data_iga_longer_bifi)
planer_data_iga_longer_bact = bin2months(planer_data_iga_longer_bact)
planer_data_iga_longer_clos = bin2months(planer_data_iga_longer_clos)

# Create the boxplot
p_e = ggplot(planer_data_iga_longer_ento, aes(x = factor(bin_idx), y = log_ratio)) +
  geom_boxplot() +
  labs(x = "Bin Idx", y = "Log Ratio", title = "Distribution of Log Ratio by Bin") +
  # scale_x_discrete(labels = c("1" = "1", "3" = "3", "6" = "6", "9" = "9", "12" = "12")) +
  theme_minimal()

p_b = ggplot(planer_data_iga_longer_bifi, aes(x = factor(bin_idx), y = log_ratio)) +
  geom_boxplot() +
  labs(x = "Bin Idx", y = "Log Ratio", title = "Distribution of Log Ratio by Bin") +
  # scale_x_discrete(labels = c("1" = "1", "3" = "3", "6" = "6", "9" = "9", "12" = "12")) +
  theme_minimal()

p_bc = ggplot(planer_data_iga_longer_bact, aes(x = factor(bin_idx), y = log_ratio)) +
  geom_boxplot() +
  labs(x = "Bin Idx", y = "Log Ratio", title = "Distribution of Log Ratio by Bin") +
  # scale_x_discrete(labels = c("1" = "1", "3" = "3", "6" = "6", "9" = "9", "12" = "12")) +
  theme_minimal()

p_c = ggplot(planer_data_iga_longer_clos, aes(x = factor(bin_idx), y = log_ratio)) +
  geom_boxplot() +
  labs(x = "Bin Idx", y = "Log Ratio", title = "Distribution of Log Ratio by Bin") +
  # scale_x_discrete(labels = c("1" = "1", "3" = "3", "6" = "6", "9" = "9", "12" = "12")) +
  theme_minimal()

# planer_data_iga_longer_bifi_agg = aggregate(log_ratio ~ days + family, planer_data_iga_longer_bifi, FUN=mean)
# plot(planer_data_iga_longer_bifi_agg$days,planer_data_iga_longer_bifi_agg$log_ratio)
# 
# planer_data_iga_longer_ento_agg = aggregate(log_ratio ~ days + family, planer_data_iga_longer_ento, FUN=mean)
# plot(planer_data_iga_longer_ento_agg$days,planer_data_iga_longer_ento_agg$log_ratio)

# Extract boxplot values for each month bin
boxplot_values <- lapply(split(planer_data_iga_longer_ento$log_ratio, planer_data_iga_longer_ento$bin_idx), function(x) {
  stats <- boxplot.stats(x)
  list(
    lower_whisker = stats$stats[1],
    first_quartile = stats$stats[2],
    median = stats$stats[3],
    third_quartile = stats$stats[4],
    upper_whisker = stats$stats[5],
    outliers = stats$out
  )
})

# Convert the list to a more readable format (optional)
boxplot_values_df_ento <- do.call(rbind, lapply(names(boxplot_values), function(x) {
  data.frame(
    bin_idx = x,
    lower_whisker = boxplot_values[[x]]$lower_whisker,
    first_quartile = boxplot_values[[x]]$first_quartile,
    median = boxplot_values[[x]]$median,
    third_quartile = boxplot_values[[x]]$third_quartile,
    upper_whisker = boxplot_values[[x]]$upper_whisker
  )
}))
rownames(boxplot_values_df_ento) <- NULL # Clean up row names


boxplot_values <- lapply(split(planer_data_iga_longer_bifi$log_ratio, planer_data_iga_longer_bifi$bin_idx), function(x) {
  stats <- boxplot.stats(x)
  list(
    lower_whisker = stats$stats[1],
    first_quartile = stats$stats[2],
    median = stats$stats[3],
    third_quartile = stats$stats[4],
    upper_whisker = stats$stats[5],
    outliers = stats$out
  )
})
# Convert the list to a more readable format (optional)
boxplot_values_df_bifi <- do.call(rbind, lapply(names(boxplot_values), function(x) {
  data.frame(
    bin_idx = x,
    lower_whisker = boxplot_values[[x]]$lower_whisker,
    first_quartile = boxplot_values[[x]]$first_quartile,
    median = boxplot_values[[x]]$median,
    third_quartile = boxplot_values[[x]]$third_quartile,
    upper_whisker = boxplot_values[[x]]$upper_whisker
  )
}))
rownames(boxplot_values_df_bifi) <- NULL # Clean up row names

# Extract boxplot values for each month bin
boxplot_values <- lapply(split(planer_data_iga_longer_bact$log_ratio, planer_data_iga_longer_bact$bin_idx), function(x) {
  stats <- boxplot.stats(x)
  list(
    lower_whisker = stats$stats[1],
    first_quartile = stats$stats[2],
    median = stats$stats[3],
    third_quartile = stats$stats[4],
    upper_whisker = stats$stats[5],
    outliers = stats$out
  )
})

# Convert the list to a more readable format (optional)
boxplot_values_df_bact <- do.call(rbind, lapply(names(boxplot_values), function(x) {
  data.frame(
    bin_idx = x,
    lower_whisker = boxplot_values[[x]]$lower_whisker,
    first_quartile = boxplot_values[[x]]$first_quartile,
    median = boxplot_values[[x]]$median,
    third_quartile = boxplot_values[[x]]$third_quartile,
    upper_whisker = boxplot_values[[x]]$upper_whisker
  )
}))
rownames(boxplot_values_df_bact) <- NULL # Clean up row names

# Extract boxplot values for each month bin
boxplot_values <- lapply(split(planer_data_iga_longer_clos$log_ratio, planer_data_iga_longer_clos$bin_idx), function(x) {
  stats <- boxplot.stats(x)
  list(
    lower_whisker = stats$stats[1],
    first_quartile = stats$stats[2],
    median = stats$stats[3],
    third_quartile = stats$stats[4],
    upper_whisker = stats$stats[5],
    outliers = stats$out
  )
})

# Convert the list to a more readable format (optional)
boxplot_values_df_clos <- do.call(rbind, lapply(names(boxplot_values), function(x) {
  data.frame(
    bin_idx = x,
    lower_whisker = boxplot_values[[x]]$lower_whisker,
    first_quartile = boxplot_values[[x]]$first_quartile,
    median = boxplot_values[[x]]$median,
    third_quartile = boxplot_values[[x]]$third_quartile,
    upper_whisker = boxplot_values[[x]]$upper_whisker
  )
}))
rownames(boxplot_values_df_clos) <- NULL # Clean up row names

p_b
boxplot_values_df_bifi

p_e
boxplot_values_df_ento

p_bc
boxplot_values_df_bact

p_c
boxplot_values_df_clos

print(c(boxplot_values_df_ento$median[length(boxplot_values_df_ento$median)],boxplot_values_df_bifi$median[length(boxplot_values_df_bifi$median)]))
print(c(boxplot_values_df_ento$median[1],boxplot_values_df_bifi$median[1]))

df_long   = readRDS('df_long.rds')
df_long_e = df_long %>% filter(series=='mIgA levels in the infant gut lumen' & t<360 )
df_long_m = df_long %>% filter(series=='Mcell activation' & t<360 )

# Transform boxplot_values_df_ento and boxplot_values_df_bifi for ggplot
boxplot_values_df_ento$bin_idx_scaled <- 30 * as.numeric(boxplot_values_df_ento$bin_idx)
boxplot_values_df_bifi$bin_idx_scaled <- 30 * as.numeric(boxplot_values_df_bifi$bin_idx)

# Combine the boxplot data frames and add a category column
boxplot_values_df_ento$category <- 'Ento'
boxplot_values_df_bifi$category <- 'Bifi'
combined_boxplot_df = rbind(boxplot_values_df_ento, boxplot_values_df_bifi)

# Plot using ggplot2
ggplot() +
  geom_line(data = df_long_e, aes(x = seq_along(value), y = value), color = 'black') +
  geom_line(data = combined_boxplot_df, aes(x = bin_idx_scaled, y = median, color = category)) +
  scale_color_manual(values = c('Ento' = 'red', 'Bifi' = 'blue')) +
  ylim(-0.5, 1.1) +
  labs(x = 'Index', y = 'Value', color = 'Category') +
  theme_minimal()


ggplot() +
  geom_line(data = df_long_e, aes(x = t, y = value), color = "black") +
  geom_line(data = df_long_m, aes(x = t, y = value), color = "darkgray") + 
  geom_vline(xintercept = 60, colour="magenta", linetype = "longdash")+
  geom_vline(xintercept = 120, colour="magenta", linetype = "longdash")+
  geom_vline(xintercept = 360, colour="magenta", linetype = "longdash")+
  geom_boxplot(data = planer_data_iga_longer_bifi, aes(x = (30*bin_idx), y = log_ratio, group = cut(30*bin_idx, breaks = unique(30*bin_idx))), fill = "blue", alpha = 0.5) +
  geom_boxplot(data = planer_data_iga_longer_ento, aes(x = (30*bin_idx), y = log_ratio, group = cut(30*bin_idx, breaks = unique(30*bin_idx))), fill = "red", alpha = 0.5) +
  labs(x = "Days", y = "Value", title = "IgA Index") +
  theme_minimal()


planer_data_iga_longer_bifi
planer_data_iga_longer_ento


boxplot_values_df_bifi
boxplot_values_df_ento
