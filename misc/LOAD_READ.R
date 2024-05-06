# run PLOT_YAGAHI_FAMILIES_AGG_2USE_DAILY.R
yagahi_data_new        = readRDS('data_smoothed.rds')
family_levels          = c("Enterobacteriaceae","Bifidobacteriaceae","Bacteroidaceae","Clostridiales")
yagahi_data_new$family = factor(yagahi_data_new$family, levels = family_levels)

# Now arrange by 'day' and 'family'
yagahi_data_new_aligned <- yagahi_data_new %>%
  dplyr::select(day, family, smoothed_values, se, xmin, xmax) %>%
  dplyr::arrange(family, day)

yagahi_data = yagahi_data_new_aligned
taxa_array  = unique(yagahi_data$family)
taxa_array  = as.character(taxa_array)

numTaxa           = length(taxa_array)
yagahi_data_keep  = yagahi_data
index_taxa_all    = 1:length(taxa_array)

########### LOAD INTERACTION MASK VECTOR #######################################

interactionMask      = read_excel('family_properties.xlsx',sheet = 'interactions',col_names = TRUE, skip=0) # for LOAD_INIT_1698698657.R
interactionMask_df   = as.data.frame(interactionMask)
interactionMask_df[is.na(interactionMask_df)] = 0
colnames(interactionMask_df)[1] = 'effected'
interactionMask_df_reorder      = interactionMask_df[,c('effected',taxa_array)]

rownames(interactionMask_df_reorder) = interactionMask_df_reorder$effected
interactionMask_df_reorder = interactionMask_df_reorder[taxa_array,]
interactionMask_df_reorder = interactionMask_df_reorder[,2:dim(interactionMask_df_reorder)[2]]
interactionMask_vector     = as.vector(t(interactionMask_df_reorder))

numNegative = length(which(interactionMask_vector==-1))
numPositive = length(which(interactionMask_vector==+1))
numAgnostic = length(which(interactionMask_vector==0))
numGnostic  = numNegative+numPositive

################################################################################
######## FIRST, THE PART FOR MATERNAL ANTIBODY REACTIVITY ######################
milkandsolid_days    = readRDS("save_data_all_milkandsolid.rds")
day_end_milk         = round(mean(milkandsolid_days$day))-1 #172 days (5.73 months)
day_end     = day_end_in
yagahi_data = yagahi_data %>% filter(day<=day_end & day>=day_start)
days_array  = unique(yagahi_data$day)
days_array  = days_array[2:length(days_array)]
yagahi_data_small = yagahi_data[c('day','smoothed_values','family')]
yagahi_data_wider = yagahi_data_small %>% pivot_wider(names_from = family, values_from = smoothed_values)
yagahi_data_wider_keep = yagahi_data_wider
yagahi_data_wider            = yagahi_data_wider[,2:dim(yagahi_data_wider)[2]]
y0_meanSubjects              = unlist(yagahi_data_wider[1,])
abundanceArray_meanSubjects  = yagahi_data_wider
y0_meanSubjects_total        = y0_meanSubjects

####### ONLY FOR OBSERVED TAXA EARLY ON ########################################
mask_array                  = as.integer(taxa_array %in% taxa_array)
indexTaxa                   = index_taxa_all[which(mask_array == 1)]
numTaxa                     = length(taxa_array)
abundanceArray_meanSubjects = abundanceArray_meanSubjects[c(taxa_array)]
y0_meanSubjects_total       = y0_meanSubjects_total[taxa_array];
abundanceArray_meanSubjects = abundanceArray_meanSubjects[2:dim(abundanceArray_meanSubjects)[1],]

###################################################################################################
metabolism = read_excel('family_properties.xlsx',sheet = 'metabolism',col_names = TRUE, skip=0)
metabolism = metabolism %>% filter(Family %in% taxa_array)
O2Dependency_vector    = metabolism$o2
TLR9_vector            = metabolism$TLR9
TLR4_vector            = metabolism$TLR4
gramstain_vector       = metabolism$gr
SCFA_vector            = metabolism$SCFA
addWithFood_vector     = metabolism$withfood

