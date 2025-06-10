if (exists('mIgA_in_keep')){
  keep_miga_vec = mIgA_in_keep
}else{
  keep_miga_vec = data_list_sim$theta_in[(numTaxa + numAgnostic + numGnostic + 1):(numTaxa + numAgnostic + numGnostic + numTaxa)];
}


library(gridExtra)
my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
# BC, B, C, E


#### PLASMA AFFINITY
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
out_df_b   = y_out_df_use_l %>% filter(taxa=='Bifidobacteriaceae')
out_df_bc  = y_out_df_use_l %>% filter(taxa=='Bacteroidaceae')
out_df_c   = y_out_df_use_l %>% filter(taxa=='Clostridiales')


### if dydt keeps daily
variable       = 'selection_thresholds'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  dplyr::mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
                taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array), 
                              taxa_array[taxa_index], 
                              taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

out_df_sth      = y_out_df_use_l
out_df_sth_bact = out_df_sth %>% filter(taxa=='Bacteroidaceae')
out_df_sth_bifi = out_df_sth %>% filter(taxa=='Bifidobacteriaceae')
out_df_sth_clos = out_df_sth %>% filter(taxa=='Clostridiales')


out_df_sth_e   = out_df_sth %>% filter(taxa=='Enterobacteriaceae' & days>0)
out_df_sth_ne  = out_df_sth %>% filter(taxa!='Enterobacteriaceae' & days>0)
colnames(out_df_sth_e)[3]='selection_threshold'
colnames(out_df_sth_ne)[3]='selection_threshold'
out_df_e  = merge(out_df_e,out_df_sth_e,by=c('days','taxa'))
out_df_ne = merge(out_df_ne,out_df_sth_ne,by=c('days','taxa'))
out_df_ene= rbind(out_df_e,out_df_ne)

graphics.off()

out_df_ene = out_df_ene %>% dplyr::rowwise() %>% dplyr::mutate(miga = ifelse(taxa=='Enterobacteriaceae',keep_miga_vec[1],
                                                                             ifelse(taxa=='Bifidobacteriaceae',keep_miga_vec[2],
                                                                                    ifelse(taxa=='Bacteroidaceae',keep_miga_vec[3],
                                                                                           ifelse(taxa=='Clostridiales',keep_miga_vec[4],NA)))))


out_df_ene_long = out_df_ene %>% pivot_longer(!c('days','taxa'), values_to = 'value', names_to = 'type')
out_df_ene_long = out_df_ene_long %>%  filter(type != "selection_threshold")

if (exists("mv_in")) {
  mv=mv_in
} else {
  mv=6
  if(max(out_df_ene_long$value)>mv){
    mv = round(1.2*max(out_df_ene_long$value))
  }
}

# Ensure 'taxa' and 'type' are factors and set their levels as needed
out_df_ene_long$taxa <- factor(out_df_ene_long$taxa)  # Adjust with desired levels if necessary
out_df_ene_long$type <- factor(out_df_ene_long$type, levels = c("median", "miga"))

##Plot
p_aff_log <- ggplot(out_df_ene_long, aes(x = days, y = value, color = taxa, linetype = type)) +
  geom_line(size = 1) +
  theme_minimal() +
  scale_linetype_manual(values = c("median" = "solid", "miga" = "dotted"),
                        labels = c("median" ="eSIgA", "miga"="mSIgA")) +
  labs(title = "Average eSIgA Affinity", x = "DOL", y = "") +
  scale_color_manual(values = my_colors[1:4]) +  # Ensure my_colors is defined
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) +
  scale_y_log10(minor_breaks = seq(0, 6, .2) * c(1e-2, 1e-1, 1e0))

p_aff_log = p_aff_log + theme(legend.title = element_blank(),
                               legend.position = c(1, 0.13),
                               legend.justification = c("right", "bottom"),
                               legend.box.just = "right",
                               # legend.margin = margin(t = -1, b = 0),
                               legend.box.background = element_rect(colour = NA, fill = NA),
                               plot.background = element_rect(fill = NA, colour = NA),
                               panel.background = element_rect(fill = NA, colour = NA))

# If you need to manually adjust the order of the legends
p_aff_log = p_aff_log + guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

out_df_ene_long_pick = out_df_ene_long %>% filter(type=='median' & value>0)
min(out_df_ene_long_pick$days)

aff_end = out_df_ene_long_pick %>% filter(days==max(days))

p_aff_log = p_aff_log + coord_cartesian(xlim = c(min(out_df_ene_long_pick$days), max(data_list_sim$ts_pred)), ylim = c(0.01, mv))

#### eIgA COUNT
variable       = 'eIgA'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

out_df_eiga    = y_out_df_use_l
y_out_df_use_l = y_out_df_use_l %>%
  dplyr::mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
                taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array), 
                              taxa_array[taxa_index], 
                              taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_out_df_use_l_eiga = y_out_df_use_l

# Plotting

p_eiga = ggplot(y_out_df_use_l, aes(x = days, y = median, color = taxa)) +
  geom_line(size = 1) +  theme_minimal() +
  labs(title = "eSIgA Concentration", x = "DOL", y = "") +
  scale_color_manual(values = my_colors) +  
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) 

p_eiga = p_eiga + theme(legend.title = element_blank(),
                        # legend.text = element_text(size = 10),
                        # axis.title = element_text(size = 10),
                        # axis.text = element_text(size = 10),
                        legend.position = c(0.28, 0.62), # This needs adjustment
                        legend.justification = c("right", "bottom"), # Anchor point
                        legend.box.just = "right",
                        # legend.margin = margin(),
                        legend.box.background = element_rect(colour = NA, fill = NA), # Make legend background transparent
                        # Optional: Make the entire plot background transparent
                        plot.background = element_rect(fill = NA, colour = NA),
                        panel.background = element_rect(fill = NA, colour = NA))


############


out_df_ene_last = out_df_ene %>% filter(days==735)

p_aff_bar = ggplot(out_df_ene_last, aes(x = taxa, y = median, fill = taxa)) +
  geom_col(position = position_dodge()) +
  coord_flip() +
  labs(x = "", y = "", title = "SIgA Affinity") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) +
  scale_fill_manual(values = my_colors[1:4]) +
  theme(legend.position = "none")


#### plasma COUNT

variable       = 'plasma_Bcell_count'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  dplyr::mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
                taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array),
                              taxa_array[taxa_index],
                              taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_out_df_use_plasma = y_out_df_use_l

p_plasma = ggplot(y_out_df_use_l, aes(x = days, y = median, color = taxa)) +
  geom_line(size = 1) +  theme_minimal() +
  labs(title = "IgA+ Cell Count", x = "DOL", y = "") +
  scale_color_manual(values = my_colors) +  
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) 


#### circulating COUNT

variable       = 'circulating_Bcell_count'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  dplyr::mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
                taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array),
                              taxa_array[taxa_index],
                              taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_out_df_use_circulating = y_out_df_use_l

p_circulating = ggplot(y_out_df_use_l, aes(x = days, y = median, color = taxa)) +
  geom_line(size = 1) +  theme_minimal() +
  labs(title = "Circulating B Cell Count", x = "DOL", y = "") +
  scale_color_manual(values = my_colors) +  
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) 


#### naive COUNT

variable       = 'naive_Bcell_count'
y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa)))]
y_out_df_use_l = y_out_df_use %>% pivot_longer(!days, names_to = 'taxa',values_to = 'median')

y_out_df_use_l <- y_out_df_use_l %>%
  dplyr::mutate(taxa_index = as.numeric(str_extract(taxa, "[0-9]+")), # Extract the numeric part
                taxa = ifelse(!is.na(taxa_index) & taxa_index <= length(taxa_array),
                              taxa_array[taxa_index],
                              taxa)) %>%
  dplyr::select(-taxa_index) # Remove the temporary index column

y_out_df_use_naive = y_out_df_use_l

p_naive = ggplot(y_out_df_use_l, aes(x = days, y = median, color = taxa)) +
  geom_line(size = 1) +  theme_minimal() +
  labs(title = "Naive B Cell Count", x = "DOL", y = "") +
  scale_color_manual(values = my_colors) +  
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) 

