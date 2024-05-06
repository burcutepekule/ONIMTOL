library(tidyverse)
library(dplyr)
library(readxl)
library(writexl)
sigmoid = function(x) 1 / (1 + exp(-x))
smoothed_data_all     = readRDS('smoothed_data_all.RDS')
all_families          = unique(smoothed_data_all$family)
smoothed_data_all_adj = c()

for(family_pick in all_families){
  print(family_pick)
  
  smoothed_data_pick  = smoothed_data_all %>% filter(family==family_pick)
  smoothed_values_vec = smoothed_data_pick$smoothed_values
  midpoint = 172 # transition to MF
  scale    = 5 # Adjust this to control the smoothness of the transition
  weights  = sigmoid((smoothed_data_pick$day - midpoint) / scale)
  
  if(family_pick %in% c('Bacteroidaceae','Clostridiales')){
    smoothed_data_pick$smoothed_values = (smoothed_data_pick$smoothed_values)*weights
  }
  smoothed_data_all_adj=rbind(smoothed_data_all_adj,smoothed_data_pick)
}

smoothed_data_all_adj    = as.data.frame(smoothed_data_all_adj)
smoothed_data_all_adj_e  = smoothed_data_all_adj %>% dplyr::filter(family=='Enterobacteriaceae')
smoothed_data_all_adj_b  = smoothed_data_all_adj %>% dplyr::filter(family=='Bifidobacteriaceae')
smoothed_data_all_adj_bc = smoothed_data_all_adj %>% dplyr::filter(family=='Bacteroidaceae')
smoothed_data_all_adj_c  = smoothed_data_all_adj %>% dplyr::filter(family=='Clostridiales')
ggp_data_use_all         = smoothed_data_all_adj


ggp_data_use_all_wide = ggp_data_use_all[c('day','family','smoothed_values')] %>% pivot_wider(names_from = family, values_from = smoothed_values)
ggp_data_use_all_wide = ggp_data_use_all_wide[order(ggp_data_use_all_wide$day),]

# new transition to MF
new_172=min(which(ggp_data_use_all_wide[,4]>0.1)[1], which(ggp_data_use_all_wide[,5]>0.1)[1]) 
saveRDS(new_172,'new_172.rds')
ggp_data_use_all_wide[1:new_172,c(4,5)]=0

ggp_data_use_all = ggp_data_use_all %>%
  mutate(smoothed_values = case_when(
    family %in% c("Clostridiales", "Bacteroidaceae") & day < new_172 ~ 0,
    TRUE ~ smoothed_values
  ))


total_abundance_scale   = rowSums(ggp_data_use_all_wide[,2:5]);
saveRDS(total_abundance_scale,'total_abundance_scale.rds')

interpolated_data_extended = readRDS('interpolated_data_extended.rds') 
colnames(interpolated_data_extended)=c('day','family','mult') 
interpolated_data_extended$mult = 1e-11*interpolated_data_extended$mult # normalize with 1e-11

interpolated_data_extended_single = interpolated_data_extended %>% filter(family=='Enterobacteriaceae')# same total scaling for all taxa, names are there to make dataframe merge easier
data_smooth   = interpolated_data_extended_single$mult
fit           = loess(data_smooth ~ seq_along(data_smooth), span = 1)
smoothed_data = predict(fit)

interpolated_data_extended    = as.data.frame(interpolated_data_extended)
interpolated_data_extended_e  = interpolated_data_extended %>% filter(family=='Enterobacteriaceae')
interpolated_data_extended_b  = interpolated_data_extended %>% filter(family=='Bifidobacteriaceae')
interpolated_data_extended_bc = interpolated_data_extended %>% filter(family=='Bacteroidaceae')
interpolated_data_extended_c  = interpolated_data_extended %>% filter(family=='Clostridiales')
interpolated_data_extended_e$mult  = smoothed_data
interpolated_data_extended_b$mult  = smoothed_data
interpolated_data_extended_bc$mult = smoothed_data
interpolated_data_extended_c$mult  = smoothed_data

interpolated_data_extended = rbind(interpolated_data_extended_e,interpolated_data_extended_b,
                                   interpolated_data_extended_bc,interpolated_data_extended_c)
interpolated_data_extended = as.data.frame(interpolated_data_extended)

ggp_data_use_all_keep = ggp_data_use_all
ggp_data_use_all = merge(interpolated_data_extended,ggp_data_use_all_keep,by=c('day','family'))
ggp_data_use_all = ggp_data_use_all %>% dplyr::rowwise() %>% dplyr::mutate(smoothed_values=0.01*smoothed_values*mult) # remove percentage

saveRDS(ggp_data_use_all,'data_smoothed.rds')

yagahi_data = ggp_data_use_all
taxa_array  = unique(yagahi_data$family)
yagahi_data_keep  = yagahi_data
days_array        = unique(yagahi_data$day)
yagahi_data_small = yagahi_data[c('day','smoothed_values','family')]
yagahi_data_wider = yagahi_data_small %>%
  pivot_wider(names_from = family, values_from = smoothed_values)
yagahi_data_wider_keep = yagahi_data_wider
yagahi_data_wider = yagahi_data_wider[order(yagahi_data_wider$day),]
yagahi_data_wider = yagahi_data_wider[,2:dim(yagahi_data_wider)[2]]

abundanceArray_meanSubjects   = yagahi_data_wider
abundanceArray_meanSubjects_c = abundanceArray_meanSubjects

# PLOT EVERYTHING TO SEE FIRST
abundanceArray_meanSubjects_use_c     = abundanceArray_meanSubjects_c
abundanceArray_meanSubjects_use_c$day = rownames(abundanceArray_meanSubjects_use_c)
abundanceArray_meanSubjects_longer_c  = abundanceArray_meanSubjects_use_c %>% pivot_longer(!day,names_to = 'taxa', values_to = 'abundance')
abundanceArray_meanSubjects_longer_c$day=as.numeric(abundanceArray_meanSubjects_longer_c$day)

