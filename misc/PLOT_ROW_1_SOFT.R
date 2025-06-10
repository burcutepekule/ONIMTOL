##### SS Abundances

# Reshape the data from wide to long format

pdc_data_use_local = pdc_data_use %>% filter(bm_cutoff %in% seq(0,720,15) & mf_length %in% seq(0,720,15))
pdc_data_long = pivot_longer(pdc_data_use_local, cols = c(aff_e, aff_bc, aff_b, aff_c), names_to = "variable", values_to = "value")

# Calculate the new variables directly within pdc_data_long
pdc_data_long <- pdc_data_long %>%
  mutate(value_coating = 1 / (1 + value),
         value_killing = 1 - 1 / (1 + value))

# Reshape data: Transforming `value_coating` and `value_killing` into long format
pdc_data_longer <- pdc_data_long %>%
  pivot_longer(cols = c(value_coating, value_killing), names_to = "value_type", values_to = "value_plot")

pdc_data_longer_c = pdc_data_longer %>% dplyr::filter(value_type=='value_coating')
pdc_data_longer_k = pdc_data_longer %>% dplyr::filter(value_type=='value_killing')

pdc_data_longer_b_c  = pdc_data_longer %>% dplyr::filter(variable=='aff_b' & value_type=='value_coating')
pdc_data_longer_b_k  = pdc_data_longer %>% dplyr::filter(variable=='aff_b' & value_type=='value_killing')
pdc_data_longer_bc_c = pdc_data_longer %>% dplyr::filter(variable=='aff_bc' & value_type=='value_coating')
pdc_data_longer_bc_k = pdc_data_longer %>% dplyr::filter(variable=='aff_bc' & value_type=='value_killing')
pdc_data_longer_c_c  = pdc_data_longer %>% dplyr::filter(variable=='aff_c' & value_type=='value_coating')
pdc_data_longer_c_k  = pdc_data_longer %>% dplyr::filter(variable=='aff_c' & value_type=='value_killing')
pdc_data_longer_e_c  = pdc_data_longer %>% dplyr::filter(variable=='aff_e' & value_type=='value_coating')
pdc_data_longer_e_k  = pdc_data_longer %>% dplyr::filter(variable=='aff_e' & value_type=='value_killing')

fk_b  = 1-length(which(is.na(pdc_data_longer_b_k$value_plot)==TRUE))/dim(pdc_data_longer_b_c)[1] # killing
fc_b  = 1-length(which(is.na(pdc_data_longer_b_c$value_plot)==TRUE))/dim(pdc_data_longer_b_c)[1] # coating
fk_bc = 1-length(which(is.na(pdc_data_longer_bc_k$value_plot)==TRUE))/dim(pdc_data_longer_bc_c)[1] # killing
fc_bc = 1-length(which(is.na(pdc_data_longer_bc_c$value_plot)==TRUE))/dim(pdc_data_longer_bc_c)[1] # coating
fk_c  = 1-length(which(is.na(pdc_data_longer_c_k$value_plot)==TRUE))/dim(pdc_data_longer_c_c)[1] # killing
fc_c  = 1-length(which(is.na(pdc_data_longer_c_c$value_plot)==TRUE))/dim(pdc_data_longer_c_c)[1] # coating
fk_e  = 1-length(which(is.na(pdc_data_longer_e_k$value_plot)==TRUE))/dim(pdc_data_longer_e_c)[1] # killing
fc_e  = 1-length(which(is.na(pdc_data_longer_e_c$value_plot)==TRUE))/dim(pdc_data_longer_e_c)[1] # coating

# Plotting
p_bp = ggplot(pdc_data_longer, aes(x = variable, y = value_plot, fill = variable, alpha = value_type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(
    "aff_b" = "#074fa8",
    "aff_bc" = "#2bc5d8",
    "aff_c" = "#ff9aa1",
    "aff_e" = "#f20026"
  )) +
  scale_alpha_manual(values = c(
    "value_coating" = 0.3,  # More transparent
    "value_killing" = 1    # Less transparent (fully opaque)
  )) +
  labs(title = "Distribution of Masking (M) and Neutralizing (N) Rate of eSIgA\nAfter 2 Years (DOL 735)\n", 
       x = "Taxa", 
       y = "Rate of IgA Function ") +
  scale_x_discrete(labels = c("Bifidobacteriaceae", "Bacteroidaceae", "Clostridiales", "Enterobacteriaceae")) +
  theme_minimal() + ylim(0.0,1) +
  theme( plot.title = element_text(size = 12), # Adjust the bottom margin to add padding
         axis.text.x = element_text(angle = 0, hjust = 0.5, size=11), legend.position = "none")

offset=0.20
size_inplot = 4
p_bp = p_bp + 
  geom_text(aes(x = 1-offset, y = .95, label = "(M)"), vjust = -0.5, size = size_inplot) +
  geom_text(aes(x = 1+offset, y = .95, label = "(N)"), vjust = -0.5, size = size_inplot) +
  geom_text(aes(x = 2-offset, y = .95, label = "(M)"), vjust = -0.5, size = size_inplot) +
  geom_text(aes(x = 2+offset, y = .95, label = "(N)"), vjust = -0.5, size = size_inplot) +
  geom_text(aes(x = 3-offset, y = .95, label = "(M)"), vjust = -0.5, size = size_inplot) +
  geom_text(aes(x = 3+offset, y = .95, label = "(N)"), vjust = -0.5, size = size_inplot) +
  geom_text(aes(x = 4-offset, y = .95, label = "(M)"), vjust = -0.5, size = size_inplot) +
  geom_text(aes(x = 4+offset, y = .95, label = "(N)"), vjust = -0.5, size = size_inplot)

pdc_data_keep_final = pdc_data_use[c('aff_b','aff_bc','aff_c')]
pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status = ifelse((aff_b<1 & aff_bc<1 & aff_c<1), 1, 0)) # TOLERANCE

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status_b = ifelse((aff_b>1), 0, 1)) # B-HYPER

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status_bc = ifelse((aff_bc>1), 0, 1)) # BC-HYPER

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status_c = ifelse((aff_c>1 ), 0, 1)) # C-HYPER


# Simplified dataframe creation
df_long = data.frame(
  Value = c(
    100 * sum(pdc_data_keep_final$final_status) / nrow(pdc_data_keep_final),
    100 * (1 - sum(pdc_data_keep_final$final_status) / nrow(pdc_data_keep_final))
  ),
  Status = c("Tolerant", "Hyperreactive")
)

# Set the order of the factor for Status
df_long$Status = factor(df_long$Status, levels = c("Hyperreactive","Tolerant"))

p_final = ggplot(df_long, aes(x = Status, y = Value, fill = Status)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "", y = "Value (%)", title = "") +
  theme_minimal() +
  scale_fill_manual(values = c( "Hyperreactive" = "gray10","Tolerant" = "gray"),
                    guide = guide_legend(reverse = TRUE)) + # Reverse legend order
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  theme(
    legend.text = element_text(size = 12), # Increase legend text size
    axis.title.y = element_text(size = 12), # Increase y-axis title size
    axis.text.y = element_text(size = 12) # Increase y-axis text size
  )+
  scale_y_continuous(breaks = seq(0, 60, by = 10)) # Adjust y-axis ticks


pdc_data_keep_final = pdc_data_use_local[c('bm_cutoff','mf_length','aff_b','aff_bc','aff_c')]

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status = ifelse((aff_b<1 & aff_bc<1 & aff_c<1), 1, 0)) # TOLERANCE

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status_b = ifelse((aff_b>1), 0, 1)) # B-HYPER

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status_bc = ifelse((aff_bc>1), 0, 1)) # BC-HYPER

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>% 
  mutate(final_status_c = ifelse((aff_c>1 ), 0, 1)) # C-HYPER

pdc_data_keep_final = pdc_data_keep_final[c('bm_cutoff','mf_length','final_status','final_status_b','final_status_bc','final_status_c')]

pdc_data_keep_final$bm_cutoff <- as.factor(pdc_data_keep_final$bm_cutoff)
pdc_data_keep_final$mf_length <- as.factor(pdc_data_keep_final$mf_length)

pdc_data_keep_final$mf_length <- as.numeric(as.character(pdc_data_keep_final$mf_length))
pdc_data_keep_final$bm_cutoff <- as.numeric(as.character(pdc_data_keep_final$bm_cutoff))

p_hm_no_label = ggplot(pdc_data_keep_final, aes(x = bm_cutoff, y = mf_length)) +
  geom_tile(aes(fill = factor(final_status)), color = '#d1e2e2') +
  scale_fill_manual(values = c("0" = "gray10", "1" = "#edf3f3"), 
                    name = "Functionality", labels = c("0" = "Hyperreactive", "1" = "Tolerant")) +
  scale_y_continuous(breaks = seq(0, max(pdc_data_keep_final$mf_length, na.rm = TRUE), by = 30),
                     labels = seq(0, max(pdc_data_keep_final$mf_length, na.rm = TRUE), by = 30)) +
  scale_x_continuous(breaks = seq(0, max(pdc_data_keep_final$bm_cutoff, na.rm = TRUE), by = 30),
                     labels = seq(0, max(pdc_data_keep_final$bm_cutoff, na.rm = TRUE), by = 30))+
  theme_minimal() +
  labs(x = "EBF Duration (days)", y = "MF Duration (days)", title = "Functionality of the Endogenous Immune System\nAfter 2 Years (DOL 735)\n") +
  theme(plot.title = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))

# p_hm <- p_hm + theme(plot.margin = unit(c(-100, 0, 0, 30), "pt"))
p_hm = p_hm_no_label + geom_text(aes(x = 120, y = 180, label = paste0(round(df_long$Value[1],1),'%')), color = 'gray10', size = 5, angle = 0)
p_hm = p_hm + geom_text(aes(x = 20, y = 15, label = paste0(round(df_long$Value[2],1),'%')), color = '#edf3f3', size = 5, angle = 0)

############ HEATMAP WITH CTS VALUES

pdc_data_keep_final = pdc_data_use_local[c('bm_cutoff','mf_length','aff_b','aff_bc','aff_c')]

pdc_data_keep_final = pdc_data_keep_final %>% rowwise() %>%
  mutate(final_value_symb = mean(c(aff_b,aff_bc,aff_c))) # CTS

pdc_data_keep_final$bm_cutoff <- as.numeric(as.character(pdc_data_keep_final$bm_cutoff))
pdc_data_keep_final$mf_length <- as.numeric(as.character(pdc_data_keep_final$mf_length))

pdc_data_keep_final = pdc_data_keep_final %>% dplyr::rowwise() %>% dplyr::mutate(value_coating = 1 / (1 + final_value_symb))
pdc_data_keep_final = pdc_data_keep_final %>% dplyr::rowwise() %>% dplyr::mutate(value_killing = 1-value_coating)

p_hm_cts_no_label = ggplot(pdc_data_keep_final, aes(x = bm_cutoff, y = mf_length)) +
  geom_tile(aes(fill = value_coating), color = '#d1e2e2') +
  scale_fill_gradientn(
    colors = c("#E52B50","gray95","#0000CD"),
    values = c(0, 0.5, 1),
    name = ""
  )+
  scale_y_continuous(breaks = seq(0, max(pdc_data_keep_final$mf_length, na.rm = TRUE), by = 30),
                     labels = seq(0, max(pdc_data_keep_final$mf_length, na.rm = TRUE), by = 30)) +
  scale_x_continuous(breaks = seq(0, max(pdc_data_keep_final$bm_cutoff, na.rm = TRUE), by = 30),
                     labels = seq(0, max(pdc_data_keep_final$bm_cutoff, na.rm = TRUE), by = 30))+
  theme_minimal() +
  labs(x = "EBF Duration (days)", y = "MF Duration (days)", title = "Average masking rate for symbiotic commensals\nafter 2 years (DOL 735)\n") +
  theme(plot.title = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1))

# p_hm <- p_hm + theme(plot.margin = unit(c(-100, 0, 0, 30), "pt"))
# p_hm_cts = p_hm_cts_no_label + geom_text(aes(x = 120, y = 180, label = paste0(round(df_long$Value[1],1),'%')), color = 'gray10', size = 5, angle = 0)
# p_hm_cts = p_hm_cts + geom_text(aes(x = 20, y = 15, label = paste0(round(df_long$Value[2],1),'%')), color = '#edf3f3', size = 5, angle = 0)
p_hm_cts = p_hm_cts_no_label

# ########################## HEATMAP WITH TRIANGLE TILES
