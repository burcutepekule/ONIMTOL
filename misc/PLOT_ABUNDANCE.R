library(gridExtra)
library(ggpattern) 
library(patchwork)
library(png)
library(grid)

my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
# BC, B, C, E

colnames(abundanceArray_meanSubjects) = taxa_array
abundanceArray_meanSubjects <- as.data.frame(abundanceArray_meanSubjects) %>%
  dplyr::mutate(across(c(Enterobacteriaceae, Bifidobacteriaceae, Bacteroidaceae, Clostridiales), 
                       ~pmax(., 0, na.rm = TRUE))) 


t_data_full    = seq(data_list_sim$t0+1,dim(abundanceArray_meanSubjects)[1]+data_list_sim$t0,1)
variable       = 'uncoated_lumen'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
         taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array), 
                       taxa_array[taxa_index], 
                       taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_uc_lumen = y_out_df_use_l

variable       = 'coated_lumen'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
         taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array), 
                       taxa_array[taxa_index], 
                       taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_c_lumen = y_out_df_use_l

variable       = 'killed_feces_cumul'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_1 = y_out_df_use[1,]
# Apply differential calculation to each column except 'days'
for (col in names(y_out_df_use)[-1]) {  # Skipping the first column ('days')
  y_out_df_use[[col]] <- c(NA, diff(y_out_df_use[[col]]))
}
y_out_df_use[1,] = y_out_df_use_1
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
         taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array), 
                       taxa_array[taxa_index], 
                       taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_k_feces = y_out_df_use_l

colnames(y_uc_lumen)[3]='median_uc'
colnames(y_c_lumen)[3] ='median_c'
colnames(y_k_feces)[3] ='median_k'

y_lumen = merge(y_k_feces, merge(y_uc_lumen, y_c_lumen, by=c('days','taxa')), by=c('days','taxa'))

##### SAMPLED ANTIGENS IN THE PPs? ######

variable       = 'antigens_sampled_uncoated'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_1 = y_out_df_use[1,]
# Apply differential calculation to each column except 'days'
for (col in names(y_out_df_use)[-1]) {  # Skipping the first column ('days')
  y_out_df_use[[col]] <- c(NA, diff(y_out_df_use[[col]]))
}
y_out_df_use[1,] = y_out_df_use_1
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
         taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array), 
                       taxa_array[taxa_index], 
                       taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_sampled_PP = y_out_df_use_l
##### SAMPLED ANTIGENS IN THE LYMPH NODE? ######


srate_coated   = data_list_sim$srate_c_in
srate_uncoated = data_list_sim$srate_uc_in
srate_killed   = data_list_sim$srate_k_in

y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_tot = srate_coated*median_c+srate_killed*median_k+srate_uncoated*median_uc)
y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_ck_tot = srate_coated*median_c+srate_killed*median_k)
y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_c_tot  = srate_coated*median_c)
y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_k_tot  = srate_killed*median_k)
y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_uc_tot = srate_uncoated*median_uc)

y_lumen_clos = y_lumen %>% filter(taxa=='Clostridiales')

observations_long <- abundanceArray_meanSubjects %>%
  dplyr::mutate(time = row_number()) %>%
  pivot_longer(cols = -time, names_to = "taxa", values_to = "Observation")

feces_daily_df = y_lumen[c('days','taxa','median_tot')]
colnames(feces_daily_df)[3]='median'
feces_daily_df_wide   = feces_daily_df %>% pivot_wider(names_from = taxa, values_from = median)
feces_daily_df_wide   = feces_daily_df_wide[order(feces_daily_df_wide$days),]

# Add a column to distinguish between Model and observation
feces_daily_df$type <- "Model"
observations_long$type <- "Observation"
colnames(observations_long)[1]='days'
colnames(observations_long)[3]='value'
colnames(feces_daily_df)[3]='value'

# Combine the datasets
combined_data <- rbind(feces_daily_df, observations_long)

# PLOT THE RELATIVE ABUNDANCES?
total_abundance_scale = readRDS('total_abundance_scale.rds')
# DILUTE
total_abundance_scale      = total_abundance_scale[data_list_sim$t_data]

total_abundance_scale = c(total_abundance_scale,rep(total_abundance_scale[length(total_abundance_scale)],dim(feces_daily_df_wide)[1]-length(total_abundance_scale)))

feces_daily_df_wide_norm           = feces_daily_df_wide
feces_daily_df_wide_norm_abundance = feces_daily_df_wide_norm[,2:5]
feces_daily_df_wide_norm_abundance$Enterobacteriaceae =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Enterobacteriaceae
feces_daily_df_wide_norm_abundance$Bifidobacteriaceae =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Bifidobacteriaceae
feces_daily_df_wide_norm_abundance$Bacteroidaceae     =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Bacteroidaceae
feces_daily_df_wide_norm_abundance$Clostridiales      =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Clostridiales


feces_daily_df_wide_norm = cbind(feces_daily_df_wide[,1],feces_daily_df_wide_norm_abundance)

t_data_full                      = unique(observations_long$days)
observations_wide                = cbind(t_data_full,abundanceArray_meanSubjects)
observations_wide_norm_abundance = observations_wide
colnames(observations_wide_norm_abundance)[1]= 'time'
total_abundance_scale = readRDS('total_abundance_scale.rds')
total_abundance_scale = total_abundance_scale[3:length(total_abundance_scale)]


observations_wide_norm_abundance$Enterobacteriaceae =  (0.01*total_abundance_scale/rowSums(observations_wide[,2:5]))*observations_wide_norm_abundance$Enterobacteriaceae
observations_wide_norm_abundance$Bifidobacteriaceae =  (0.01*total_abundance_scale/rowSums(observations_wide[,2:5]))*observations_wide_norm_abundance$Bifidobacteriaceae
observations_wide_norm_abundance$Bacteroidaceae     =  (0.01*total_abundance_scale/rowSums(observations_wide[,2:5]))*observations_wide_norm_abundance$Bacteroidaceae
observations_wide_norm_abundance$Clostridiales      =  (0.01*total_abundance_scale/rowSums(observations_wide[,2:5]))*observations_wide_norm_abundance$Clostridiales


feces_daily_df_norm    = feces_daily_df_wide_norm %>% pivot_longer(!days, names_to = "taxa", values_to = "value")
observations_long_norm = observations_wide_norm_abundance %>% pivot_longer(!time, names_to = "taxa", values_to = "value")

feces_daily_df_norm$type     <- "Model"
observations_long_norm$type <- "Observation"
colnames(observations_long_norm)[1]='days'
# Combine the datasets
combined_data_norm <- rbind(feces_daily_df_norm, observations_long_norm)
# Plotting

if (!exists("plot1")) {
  plot1 = 0  # plot1 is for the main paper figure where out-of-sample data is labelled
}

if(plot1 == 1){
  library(magick)
  
  img   = readPNG("~/Dropbox/criticalwindow/code/R/RStan/eye.png")
  
  # Convert the image to a rasterGrob
  img_grob = rasterGrob(img, interpolate = TRUE)
  
  p_rel = ggplot(combined_data_norm, aes(x = days, y = value, color = taxa, linetype = type)) +
    geom_line(size = 0.8) +  
    theme_minimal() +
    labs(title = "Relative abundances in fecal samples", x = "DOL", y = "Ratio") +
    scale_color_manual(values = my_colors) +
    scale_linetype_manual(values = c("Model" = "11", "Observation" = "solid")) +
    theme(legend.title = element_blank(), plot.title = element_text(size = 12))+
    annotate("rect", xmin=274, xmax=720, ymin=0, ymax=Inf, alpha=0.25, fill="gray") +
    annotation_custom(grob = img_grob, xmin = 274+100+215, xmax = 735, ymin = 0.78, ymax = 0.88)
  
}else{
  p_rel=ggplot(combined_data_norm, aes(x = days, y = value, color = taxa, linetype = type)) +
    geom_line(size = 0.8) +  theme_minimal() +
    labs(title = "Relative abundances in fecal samples", x = "DOL", y = "Ratio") +
    scale_color_manual(values =my_colors) +
    scale_linetype_manual(values = c("Model" = "11", "Observation" = "solid")) +
    theme(legend.title = element_blank(), plot.title = element_text(size = 12)) 
}



combined_data_norm     = as.data.frame(combined_data_norm)
lumen_daily_cuc_df_e   = combined_data_norm %>% filter(taxa=="Enterobacteriaceae" & type=='Model' & days<=14)
lumen_daily_cuc_df_d7  = lumen_daily_cuc_df_e %>% filter(days<=7)
lumen_daily_cuc_df_d7  = aggregate(value~taxa,lumen_daily_cuc_df_d7,FUN=sum)

lumen_daily_cuc_df_d14 = lumen_daily_cuc_df_e %>% filter(days<=14)
lumen_daily_cuc_df_d14 = aggregate(value~taxa,lumen_daily_cuc_df_d14,FUN=sum)

lumen_daily_cuc_df_m1 = lumen_daily_cuc_df_e %>% filter(days<=28)

lumen_daily_c_df  = y_lumen[c('days','taxa','median_c')]
lumen_daily_uc_df = y_lumen[c('days','taxa','median_uc')]

lumen_daily_c_df$coating  = 'SIgA+'
lumen_daily_uc_df$coating = 'SIgA-'
colnames(lumen_daily_c_df)[3]='median'
colnames(lumen_daily_uc_df)[3]='median'

lumen_daily_cuc_df = rbind(lumen_daily_c_df,lumen_daily_uc_df)
lumen_daily_cuc_df$median=as.numeric(lumen_daily_cuc_df$median)

seq_ticks = seq(0, 0.15, 0.05)
# if(max(lumen_daily_cuc_df$median)>0.15){
#   seq_ticks = seq(0, 0.30, 0.05)
# }
# y_up = seq_ticks[which(seq_ticks>max(lumen_daily_cuc_df$median))[1]]
y_up = round(max(lumen_daily_cuc_df$median)*1.1,2)

p_abs=ggplot(lumen_daily_cuc_df, aes(x = days, y = median, color = taxa, linetype = coating)) +
  geom_line(size = 0.8) +  theme_minimal() +
  labs(title = "Absolute abundances in gut lumen", x = "DOL", y = "(Cells/g LC) x 1e-11") +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c("SIgA+" = "dotdash", "SIgA-" = "solid")) +
  theme(legend.title = element_blank(), plot.title = element_text(size = 12))+
  scale_y_continuous(breaks = seq_ticks) + ylim(0,y_up)
p_abs

lumen_daily_cuc_df     = as.data.frame(lumen_daily_cuc_df)
lumen_daily_cuc_df_e   = lumen_daily_cuc_df %>% filter(taxa=="Enterobacteriaceae" & coating=='SIgA-' & days<=14)
lumen_daily_cuc_df_d7  = lumen_daily_cuc_df_e %>% filter(days<=7)
lumen_daily_cuc_df_d7  = aggregate(median~taxa,lumen_daily_cuc_df_d7,FUN=sum)

lumen_daily_cuc_df_d14 = lumen_daily_cuc_df_e %>% filter(days<=14)
lumen_daily_cuc_df_d14 = aggregate(median~taxa,lumen_daily_cuc_df_d14,FUN=sum)

lumen_daily_cuc_df_m1 = lumen_daily_cuc_df_e %>% filter(days<=28)

# inflammation
variable       = 'level_inflammation_cumul'
y_out_df_use   = y_out_df[c('days',variable)]
y_out_df_use_1 = y_out_df_use[1,]
# Apply differential calculation to each column except 'days'
for (col in names(y_out_df_use)[-1]) {  # Skipping the first column ('days')
  y_out_df_use[[col]] <- c(NA, diff(y_out_df_use[[col]]))
}
y_out_df_use[1,] = y_out_df_use_1
colnames(y_out_df_use)[2] = 'median'

y_out_inf  = y_out_df_use
y_out_df_use$median = y_out_df_use$median/unique(diff(y_out_df_use$days)) # proper values for differentiation

y_out_df_use_micro_abs = y_out_df_use # keep abs value
y_out_df_use$median = y_out_df_use$median/max(y_out_df_use$median) #normalized by the max value
y_out_df_use$median = y_out_df_use$median #normalized by the max value
micro_values = y_out_df_use$median


# O2
variable       = 'level_O2'
y_out_df_use   = y_out_df[c('days',variable)]
colnames(y_out_df_use)[2] = 'median'
y_O2 = y_out_df_use
O2_values = y_O2$median

# Create a data frame for plotting
df_cal = data.frame(
  t = y_out_df_use$days,
  O2 = O2_values,
  micro = micro_values
)

colnames(df_cal)[2]='O2 concentration'
colnames(df_cal)[3]='Microenvironmental stimulation'

df_long = pivot_longer(df_cal, cols = -t, names_to = "series", values_to = "value")

# Define colors and linetypes
colors = c("O2 concentration" = "#9290C3", "Microenvironmental stimulation" = "gray20")

linetypes =  c("O2 concentration" = "dashed", "Microenvironmental stimulation" = "solid")

df_long$series <- factor(df_long$series, levels = c(
  "O2 concentration",
  "Microenvironmental stimulation"
))
df_long_micro=df_long
p_micro <- ggplot(df_long, aes(x = t, y = value, color = series, linetype = series)) +
  geom_line(data = df_long, size=0.7) +
  theme_minimal() +
  scale_color_manual(values = colors, labels = c("O2 concentration", "Microenvironmental\nstimulation")) +
  scale_linetype_manual(values = linetypes, labels = c("O2 concentration", "Microenvironmental\nstimulation")) +
  labs(title = 'Dynamic variables, gut lumen', x = "DOL", y = "") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) +
  guides(size = "none") + 
  theme(legend.title = element_blank(),
        # legend.text = element_text(size = 12),
        # axis.title = element_text(size = 12),
        # axis.text = element_text(size = 12),
        legend.position = c(1, 0.65),
        legend.justification = c("right", "bottom"), # Anchor point
        legend.box.just = "right",
        # legend.margin = margin(b = 4, t=-4), # Increased bottom margin for padding
        legend.background = element_rect(colour = "gray", fill="white")) 


# HERE THINK OF PLOTS FOR COATED VS UNCOATED DISTRIBUTIONS
# below are all absolute abundances 
feces_daily_uc_df = y_lumen[c('days','taxa','median_uc_tot')]
feces_daily_ck_df = y_lumen[c('days','taxa','median_ck_tot')]
feces_daily_c_df  = y_lumen[c('days','taxa','median_c_tot')]
feces_daily_k_df  = y_lumen[c('days','taxa','median_k_tot')]
lumen_daily_c_df  = y_lumen[c('days','taxa','median_c')]
lumen_daily_uc_df = y_lumen[c('days','taxa','median_uc')]

colnames(feces_daily_uc_df)[3] ='feces_uc'
colnames(feces_daily_ck_df)[3] ='feces_ck'
colnames(feces_daily_c_df)[3]  ='feces_c'
colnames(feces_daily_k_df)[3]  ='feces_k'

colnames(lumen_daily_uc_df)[3] ='lumen_uc'
colnames(lumen_daily_c_df)[3]  ='lumen_c'

df_feces = merge(feces_daily_uc_df,feces_daily_ck_df,by=c('days','taxa'))
df_feces = merge(df_feces,feces_daily_c_df,by=c('days','taxa'))
df_feces = merge(df_feces,feces_daily_k_df,by=c('days','taxa'))

df_lumen = merge(lumen_daily_uc_df,lumen_daily_c_df,by=c('days','taxa'))

df_feces_total_uc = aggregate(feces_uc ~ days, df_feces,FUN=sum)
df_feces_total_ck = aggregate(feces_ck ~ days, df_feces,FUN=sum)
df_feces_total_c  = aggregate(feces_c ~ days, df_feces,FUN=sum)
df_feces_total_k  = aggregate(feces_k ~ days, df_feces,FUN=sum)

df_feces_total    = merge(df_feces_total_uc,df_feces_total_ck,by=c('days'))
df_feces_total    = merge(df_feces_total,df_feces_total_c,by=c('days'))
df_feces_total    = merge(df_feces_total,df_feces_total_k,by=c('days'))

df_feces_total$taxa = 'Pooled'

df_lumen_total_uc = aggregate(lumen_uc ~ days, df_lumen,FUN=sum)
df_lumen_total_c  = aggregate(lumen_c ~ days, df_lumen,FUN=sum)
df_lumen_total    = merge(df_lumen_total_uc,df_lumen_total_c,by=c('days'))
df_lumen_total$taxa = 'Pooled'

df_feces = rbind(df_feces,df_feces_total)
df_lumen = rbind(df_lumen,df_lumen_total)

diassoc_rate_base = data_list_sim$theta_in[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+17)]
out_df_aff        = out_df_aff %>% dplyr::rowwise() %>% dplyr::mutate(binding_ability=(median>0)*(diassoc_rate_base>0)*(1-diassoc_rate_base/(1+median)))
df_feces          = merge(df_feces, out_df_aff, by=c('days','taxa'))

df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(feces_ck=max(0,feces_ck))
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(feces_uc=max(0,feces_uc))
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(feces_c=max(0,feces_c))
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(feces_k=max(0,feces_k))

df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_coated_fraction=feces_ck/(feces_ck+feces_uc))
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_coated_fraction_c=feces_c/(feces_ck+feces_uc))
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_coated_fraction_k=feces_k/(feces_ck+feces_uc))
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_uncoated_fraction=feces_uc/(feces_ck+feces_uc))

df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_plus=binding_ability*IgA_coated_fraction)
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_minus=(1-binding_ability)*IgA_uncoated_fraction)
# df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_plus=IgA_coated_fraction)
# df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(IgA_minus=IgA_uncoated_fraction)
df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(iga_index=-1*(log(IgA_plus)-log(IgA_minus))/(log(IgA_plus)+log(IgA_minus)))


df_feces_e = df_feces %>% filter(taxa=='Enterobacteriaceae')


df_feces$sample = 'feces'

df_iga_index     = df_feces[c('days','taxa','iga_index')]
df_coating_ratio = df_feces[c('days','taxa','IgA_coated_fraction')]
df_coating_ratio_c = df_feces[c('days','taxa','IgA_coated_fraction_c')]
df_coating_ratio_k = df_feces[c('days','taxa','IgA_coated_fraction_k')]



df_iga_index_conv       = df_iga_index %>% filter(days==min(735,max(df_iga_index$days)))
df_coating_ratio_c_conv = df_coating_ratio_c  %>% filter(days==min(735,max(df_iga_index$days)))
df_coating_ratio_k_conv = df_coating_ratio_k  %>% filter(days==min(735,max(df_iga_index$days)))

df_coating_ratio_conv = merge(df_coating_ratio_c_conv,df_coating_ratio_k_conv, by=c('days','taxa'))
df_coating_ratio_conv = merge(df_coating_ratio_conv,df_iga_index_conv, by=c('days','taxa'))

df_long <- pivot_longer(df_coating_ratio_conv, cols = starts_with("IgA_coated_fraction_") | starts_with("iga_index"),
                        names_to = "Measurement", values_to = "Value")

df_long$taxa = factor(df_long$taxa, levels = c("Clostridiales", "Bacteroidaceae",  "Bifidobacteriaceae","Enterobacteriaceae"))

# p_coating = ggplot(df_long, aes(x = taxa, y = Value, fill = Measurement)) +
#   geom_col(position = position_dodge()) +
#   coord_flip() +
#   labs(x = "", y = "Value", title = "") +
#   theme_minimal() +
#   scale_fill_manual(values = c("IgA_coated_fraction_c" = "gray",
#                                "IgA_coated_fraction_k" = "gray10",
#                                "iga_index" = "#FC6736"),
#                     labels = c("IgA+ Fraction (C)", "IgA+ Fraction (N)", "IgA Index")) +
#   geom_hline(yintercept = 0, linetype = "dashed")



z = structure(list(Value= df_long$Value, Taxa= df_long$taxa, Measurement = df_long$Measurement), row.names = c(NA, -12L), class = c("tbl_df", "tbl", "data.frame"))

# # Create the combined plot with corrections - OLD, COMBINED
# p_coating <- ggplot(z, aes(x = Taxa, y = Value, fill = Measurement, pattern = Measurement)) +
#   geom_bar_pattern(
#     stat = "identity", position = "dodge",
#     pattern_spacing = 0.05,
#     pattern_angle = 45
#   ) +
#   coord_flip() +
#   geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Corrected line at y=0
#   scale_fill_manual(values=c('gray', 'black','#FC6736'),labels = c("SIgA+ (M)", "SIgA+ (N)", "IgA Index")) +
#   scale_pattern_manual(values=c('none', 'none', 'stripe'),labels = c("SIgA+ (M)", "SIgA+ (N)", "IgA Index")) +
#   ggpubr::theme_pubr() +
#   theme(
#     legend.position = "right",
#     panel.border = element_blank(), # Ensure the plot border is removed
#     panel.grid.major.x = element_line(color = "lightgray", size = 0.1), # Add major grid lines for the y-axis (vertical due to coord_flip)
#     panel.grid.minor.x = element_line(color = "lightgray", size = 0.1), # Optionally, add minor grid lines for the y-axis (vertical due to coord_flip)
#     panel.grid.major.y = element_line(color = "lightgray", size = 0.1), # Add major grid lines for the y-axis (vertical due to coord_flip)
#     panel.grid.minor.y = element_line(color = "lightgray", size = 0.1), # Optionally, add minor grid lines for the y-axis (vertical due to coord_flip)
#     axis.line.y = element_blank(), # Remove axis lines
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_text(size = 10),
#     panel.grid.major = element_blank(), # Optionally remove major grid lines
#     panel.grid.minor = element_blank(), # Optionally remove minor grid lines
#     plot.title = element_text(size = 12)
#   ) +
#   labs(title  = "Coating fractions and IgA Index", y = "Value") +
#   scale_y_continuous(limits = c(-0.3, 0.65), 
#                      breaks = seq(-0.2, 0.60, by = 0.2),
#                      labels = function(x) sprintf("%.1f", x))

colnames(z)[2] = 'Taxon'
z_coating = z %>% dplyr::filter(Measurement != 'iga_index')
z_iga     = z %>% dplyr::filter(Measurement == 'iga_index')

# Create the combined plot with corrections
p_coating = ggplot(z_coating, aes(x = Taxon, y = Value, fill = Measurement, pattern = Measurement)) +
  geom_bar_pattern(
    stat = "identity", position = "dodge",
    pattern_spacing = 0.05,
    pattern_angle = 45
  ) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Corrected line at y=0
  # scale_fill_manual(values=c('gray', 'black'),labels = c("SIgA+ (M)", "SIgA+ (N)")) +
  # scale_pattern_manual(values=c('none', 'none'),labels = c("SIgA+ (M)", "SIgA+ (N)")) +
  scale_fill_manual(values=c('gray', 'black'),labels = c("(M)", "(N)")) +
  scale_pattern_manual(values=c('none', 'none'),labels = c("(M)", "(N)")) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = c(1.075, 0.89),  # Position legend in top right corner
    legend.justification = c(1, 1),   # Anchor point at top right of legend
    legend.box.just = "right",        # Justify legend box
    legend.margin = ggplot2::margin(6, 6, 6, 6),
    legend.background = element_rect(fill = 'NA', color = 'NA'), # Optional: add background to legend
    legend.title = element_blank(),
    panel.border = element_blank(), # Ensure the plot border is removed
    panel.grid.major.x = element_line(color = "lightgray", size = 0.1), # Add major grid lines for the y-axis (vertical due to coord_flip)
    panel.grid.minor.x = element_line(color = "lightgray", size = 0.1), # Optionally, add minor grid lines for the y-axis (vertical due to coord_flip)
    panel.grid.major.y = element_line(color = "lightgray", size = 0.1), # Add major grid lines for the y-axis (vertical due to coord_flip)
    panel.grid.minor.y = element_line(color = "lightgray", size = 0.1), # Optionally, add minor grid lines for the y-axis (vertical due to coord_flip)
    axis.line.y = element_blank(), # Remove axis lines
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(), # Optionally remove major grid lines
    panel.grid.minor = element_blank(), # Optionally remove minor grid lines
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title  = "Coating fractions", y = "Value") +
  scale_y_continuous(limits = c(-0.01, 0.55), 
                     breaks = seq(-0.2, 0.60, by = 0.2),
                     labels = function(x) sprintf("%.1f", x))

p_iga = ggplot(z_iga, aes(x = Taxon, y = Value, fill = Measurement, pattern = Measurement)) +
  geom_bar_pattern(
    stat = "identity", position = "dodge",
    pattern_spacing = 0.05,
    pattern_angle = 45
  ) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Corrected line at y=0
  scale_fill_manual(values=c('#FC6736'),labels = c("IgA Index")) +
  scale_pattern_manual(values=c('stripe'),labels = c("IgA Index")) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "none",
    panel.border = element_blank(), # Ensure the plot border is removed
    panel.grid.major.x = element_line(color = "lightgray", size = 0.1), # Add major grid lines for the y-axis (vertical due to coord_flip)
    panel.grid.minor.x = element_line(color = "lightgray", size = 0.1), # Optionally, add minor grid lines for the y-axis (vertical due to coord_flip)
    panel.grid.major.y = element_line(color = "lightgray", size = 0.1), # Add major grid lines for the y-axis (vertical due to coord_flip)
    panel.grid.minor.y = element_line(color = "lightgray", size = 0.1), # Optionally, add minor grid lines for the y-axis (vertical due to coord_flip)
    axis.line.y = element_blank(), # Remove axis lines
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(), # Optionally remove major grid lines
    panel.grid.minor = element_blank(), # Optionally remove minor grid lines
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  labs(title  = "IgA index", y = "Value") +
  scale_y_continuous(limits = c(-0.3, 0.65), 
                     breaks = seq(-0.2, 0.60, by = 0.2),
                     labels = function(x) sprintf("%.1f", x))


df_lumen_tot = df_lumen
df_lumen_tot = df_lumen_tot %>% dplyr::rowwise() %>% dplyr::mutate(abundance = lumen_uc+lumen_c)
df_lumen_tot = df_lumen_tot %>% dplyr::rowwise() %>% dplyr::mutate(type = ifelse(taxa=='Enterobacteriaceae','P','C'))
df_lumen_tot_agg = aggregate(abundance ~ days+type, df_lumen_tot, FUN=sum)
df_lumen_tot_agg_wide = df_lumen_tot_agg %>% pivot_wider(names_from = type, values_from = abundance)
df_lumen_tot_agg_wide = df_lumen_tot_agg_wide  %>% dplyr::rowwise() %>% dplyr::mutate(pcr = P/C)

df_lumen_tot_agg_wide_m1 = df_lumen_tot_agg_wide %>% filter(days<=28)
# Plotting
p_pcr=ggplot(df_lumen_tot_agg_wide_m1, aes(x = days, y = pcr), color='gray10')+
  geom_line(size = 0.8) +  theme_minimal() +
  labs(title = "Enterobacteriaceae ratio, first month", x = "Days", y = "") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) 

