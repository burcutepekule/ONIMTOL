rm(list=ls())
library(tidyverse)
library(cowplot) 
library(readxl)
library(icesTAF)
library(plyr)
library(dplyr)
library(writexl)
library(deSolve)
library(patchwork)
library(shiny)

############### LOAD_DATA.R - BEGIN ############################################# 
day_start           = 2 # most Bifidobacteriaceae transfer is assumed to be over at day 3, indexing from 0
fit_span            = 0.5
yagahi_data = readRDS('data_smoothed.rds')
taxa_array  = unique(yagahi_data$family)
yagahi_data_keep  = yagahi_data
days_array        = unique(yagahi_data$day)
yagahi_data_small = yagahi_data[c('day','smoothed_values','family')]
yagahi_data_wider = yagahi_data_small %>%
  pivot_wider(names_from = family, values_from = smoothed_values)
yagahi_data_wider_keep = yagahi_data_wider
yagahi_data_wider = yagahi_data_wider[order(yagahi_data_wider$day),]
yagahi_data_wider = yagahi_data_wider[,2:dim(yagahi_data_wider)[2]]
abundanceArray_meanSubjects              = yagahi_data_wider
abundanceArray_meanSubjects_c            = abundanceArray_meanSubjects
abundanceArray_meanSubjects_use_c        = abundanceArray_meanSubjects_c
abundanceArray_meanSubjects_use_c$day    = rownames(abundanceArray_meanSubjects_use_c)
abundanceArray_meanSubjects_longer_c     = abundanceArray_meanSubjects_use_c %>% pivot_longer(!day,names_to = 'taxa', values_to = 'abundance')
abundanceArray_meanSubjects_longer_c$day = as.numeric(abundanceArray_meanSubjects_longer_c$day)
day_end_in          = 172+563 # max + 563
add_init_estimation = 21
new_172             = readRDS('new_172.rds')
rep_days            = new_172+add_init_estimation-3 # starts day 3
abundance_th        = 1e-4
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

mask_array                  = as.integer(taxa_array %in% taxa_array)
indexTaxa                   = index_taxa_all[which(mask_array == 1)]
numTaxa                     = length(taxa_array)
abundanceArray_meanSubjects = abundanceArray_meanSubjects[c(taxa_array)]
y0_meanSubjects_total       = y0_meanSubjects_total[taxa_array];
abundanceArray_meanSubjects = abundanceArray_meanSubjects[2:dim(abundanceArray_meanSubjects)[1],]
metabolism = read_excel('family_properties.xlsx',sheet = 'metabolism',col_names = TRUE, skip=0)
metabolism = metabolism %>% filter(Family %in% taxa_array)
O2Dependency_vector    = metabolism$o2
TLR9_vector            = metabolism$TLR9
TLR4_vector            = metabolism$TLR4
gramstain_vector       = metabolism$gr
SCFA_vector            = metabolism$SCFA
addWithFood_vector     = metabolism$withfood
############### LOAD_DATA.R - END ############################################# 

############### LOAD_DATA_LIST.R - BEGIN ############################################# 
target_timestamp = 1711378123+10 # NEW
chain_index      = 1234
# LOAD INFERRED-1 RESULTS
summary_df             = readRDS(paste0('summary_df_',target_timestamp,'_',chain_index,'.rds'))
summary_df_theta       = summary_df[grep('^theta_out\\.[0-9]+$', summary_df$variable), ]
theta_in               = summary_df_theta$median

mIgA_reactivity_vector = theta_in[(numTaxa+numAgnostic+numGnostic+1):(numTaxa+numAgnostic+numGnostic+numTaxa)]
# LOAD INFERRED-2 RESULTS
alpha_vector = c(10.062,0.089,0.048,0.040)
kappa_vector = c(314.093,1.264,1.519,2.151)
tau_new       = 180.165
tau_c         = 20.12
tau_delta     = 0.078
C_n           = 0.034
c_n           = 0.000634
C_I           = 3.43

abundanceArray_meanSubjects <- as.data.frame(abundanceArray_meanSubjects) %>%
  dplyr::mutate(across(c(Enterobacteriaceae, Bifidobacteriaceae, Bacteroidaceae, Clostridiales), 
                       ~pmax(., 0, na.rm = TRUE))) 

summary_df_theta = summary_df[grep('^theta_out\\.[0-9]+$', summary_df$variable), ]
summary_df_c     = summary_df[grep('^srate_c$', summary_df$variable), ]
summary_df_uc    = summary_df[grep('^srate_uc$', summary_df$variable), ]
summary_df_k     = summary_df[grep('^srate_k$', summary_df$variable), ]
summary_df_kappa = summary_df[grep('^kappa$', summary_df$variable), ]

theta_in    = summary_df_theta$mean
srate_c_in  = summary_df_c$mean
srate_uc_in = summary_df_uc$mean
srate_k_in  = summary_df_k$mean

mIgA_reactivity_vector = theta_in[(numTaxa+numAgnostic+numGnostic+1):(numTaxa+numAgnostic+numGnostic+numTaxa)]
drate_base             = theta_in[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+17)]

#### INPUTS ARE IN THE UI PART, SO IS data_list_sim

############### LOAD_DATA_LIST.R - END ############################################# 
############### MODEL_BASE_FUNCTIONS.R - BEGIN #####################################

div0 = function(nom, denom) {
  if (nom <= 0 | denom <= 0) {
    return(0)
  } else {
    return(nom / denom)
  }
}

div0_vector = function(nom, denom) {
  len_v = length(nom)
  out = rep(0, len_v)
  
  for (i in 1:len_v) {
    if (nom[i] <= 0 | denom <= 0) {
      out[i] = 0
    } else {
      out[i] = nom[i] / denom
    }
  }
  
  return(out)
}

div0_vector_vector = function(nom, denom) {
  len_v = length(nom)
  out = rep(0, len_v)
  for (i in 1:len_v) {
    nom_temp   = nom[i] 
    denom_temp = denom[i]
    if (nom_temp<= 0 | denom_temp<= 0) {
      out[i] = 0
    } else {
      out[i] = nom_temp/denom_temp
    }
  }
  
  return(out)
}


log10_0 = function(x) {
  if (x <= 0) {
    return(0)
  } else {
    return(log10(x))
  }
}

th_0 = function(x) {
  if (x <= 0) {
    return(0)
  } else {
    return(x)
  }
}

th_0_vector = function(x) {
  len_v = length(x)
  out = rep(0, len_v)
  
  for (i in 1:len_v) {
    if(!is.nan(x[i])){
      if (x[i] <= 0) {
        out[i] = 0
      } else {
        out[i] = x[i]
      }
    }
  }
  
  return(out)
}

sigmoid_0 = function(x) {
  if (x == 0) {
    return(0)
  } else {
    return(1 / (1+exp(-x)))
  }
}

calculate_expected_value_grid_region = function(mean1, std1, lower_bound, upper_bound, num_grid_points) {
  out = rep(0,2)
  
  # Extract lower and upper bounds for x and y
  x_lower = lower_bound[1]
  x_upper = upper_bound[1]
  total_prob = 0
  expected_value = 0
  
  if(!is.nan(mean1)){
    if(x_upper>x_lower & std1>1e-8){
      
      # Create a grid within the specified region
      x_grid = seq(x_lower, x_upper, length.out = num_grid_points)
      # Initialize sum for expected values
      sum_expected_value =0
      total_prob = 0
      
      # Iterate over the grid
      for (i in 1:(length(x_grid) - 1)) {
        # Define the corners of the current grid cell
        x1l = x_grid[i]
        x1u = x_grid[i+1]
        
        cell_mean = (x1l+x1u)/2
        cell_prob = pnorm(x1u, mean = mean1, sd = std1) -  pnorm(x1l, mean = mean1, sd = std1)
        sum_expected_value = sum_expected_value+cell_mean*cell_prob
        total_prob = total_prob+cell_prob
      }
      # Normalize by total probability
      if(total_prob<1e-3){total_prob=0}
      expected_value = div0(sum_expected_value,total_prob)
    }else{
      expected_value = 0
      total_prob = 0
    }
  }
  
  out[1] = total_prob
  out[2] = expected_value
  return(out)
}

buildInteractionMatrix = function(x_i, theta) { # this is not necessary 
  numTaxa = x_i[1]
  numAgnostic = x_i[2]
  numGnostic = x_i[3]
  counter_agnostic = 1
  counter_gnostic = 1
  counter_mask = 1
  
  interactionMat_vector_Agnostic = theta[(1+numTaxa):(numTaxa+numAgnostic)]
  interactionMat_vector_Gnostic = theta[(numTaxa+numAgnostic+1):(numTaxa+numAgnostic+numGnostic)]
  
  interactionMat = matrix(0, nrow = numTaxa, ncol = numTaxa)
  
  for (r in 1:numTaxa) {
    for (c in 1:numTaxa) {
      mask = x_i[3+counter_mask]
      if (mask == 0) {
        interactionMat[r, c] = interactionMat_vector_Agnostic[counter_agnostic]
        counter_agnostic = counter_agnostic+1
      } else {
        interactionMat[r, c] = mask*interactionMat_vector_Gnostic[counter_gnostic]
        counter_gnostic = counter_gnostic+1
      }
      counter_mask = counter_mask+1
    }
  }
  
  interactionMat_flat = as.vector(interactionMat)
  
  return(interactionMat_flat)
}


function_HMO = function(t, EBF_duration, MF_duration, level_hmo, offset) {
  
  V0 = 0.2
  a = 154
  end_time = EBF_duration+MF_duration+1 # +1 in case MF_duration=0
  b = EBF_duration 
  c1 = 40
  c2 = 40
  
  K = 1.28/2 - V0*(1 - exp(-a / c1)) # calculate the meeting point
  
  if(end_time>1){ # any breastfeeding at all - since the system starts from day 3
    if(EBF_duration>0){ # any exclusive breastfeeding
      if(b>a){
        if (t < a) {
          y = K+V0*(1 - exp(-t / c1))  # Charging or growth phase
        } else if (t < b) {
          y = K+V0*(1 - exp(-a / c1))  # Plateau
        } else {
          decline_elapsed = t - b
          total_decline_duration = end_time - b
          mirrored_time = a - (decline_elapsed*a / total_decline_duration)
          K2 = (K+V0*(1 - exp(-a / c1))) / (1 - exp(-a / c2))
          y = K2*(1 - exp(-mirrored_time / c2))  # Mirrored discharging or decline phase
          if (y < 0) y = 0
        }
      }else{
        if (t < b) {
          y = K+V0*(1 - exp(-t / c1))  # Charging or growth phase
        }  else {
          decline_elapsed = t - b
          total_decline_duration = end_time - b
          mirrored_time = b - (decline_elapsed*b / total_decline_duration)
          K2 = (K+V0*(1 - exp(-b / c1))) / (1 - exp(-b / c2))
          y = K2*(1 - exp(-mirrored_time / c2))  # Mirrored discharging or decline phase
          if (y < 0) y = 0
        }
      }
    }else{
      decline_elapsed = t - b
      total_decline_duration = end_time - b
      mirrored_time = a - (decline_elapsed*a / total_decline_duration)
      # K2 = (K+V0*(1 - exp(-b / c1))) / (1 - exp(-b / c2))
      K2 = 0.5*0.44
      y = K2*(1 - exp(-mirrored_time / c2))  # Mirrored discharging or decline phase
      if (y < 0) y = 0
    }
  }else{
    y=0
  }
  y = y + offset
  return(level_hmo*y)
}

logistic_approx = function(x, x0, k){
  return (1 / (1+exp(-k*(x - x0))));
}

function_Mcell = function(t, EBF_duration, t_kickstart, t0_mcell, level_mIgA) {
  t_transition = t0_mcell;
  k = 1;
  if(t_transition<t_kickstart){
    t_transition = t_kickstart;
  }
  return (1 / (1+exp(-k*(t - t_transition))));
}

function_solid = function(t, EBF_duration, MF_duration, t_kickstart, level_hmo, offset) {
  K = 1.280
  r = 0.005
  
  if(t<154){
    total_calories_t = function_HMO(t, 154, 308, 1, offset) # as was fit to data
  }else{
    total_calories_t = K / (1+exp(-r*(t - 154)))
  }
  
  solid_level = (total_calories_t - function_HMO(t, EBF_duration, MF_duration, level_hmo, offset))
  
  if (t < EBF_duration & level_hmo!=0) {
    solid_level = 0
  }
  
  if (t > 720) { 
    solid_level = K / (1+exp(-r*(720 - 154)))
  }
  return(solid_level)
}

function_mIgA = function(t, level_miga) {
  mIgA_level = (1 / 2.67)*(0.4+1.0321*exp(-0.0316*t))  # 2.81 to ensure max is 1
  return(level_miga*mIgA_level)
}

function_steroids = function(t, level_hmo) {
  s_level = (1 / 2.67)*(0.4+1.0321*exp(-0.0316*t))  # 2.81 to ensure max is 1
  return(level_hmo*s_level)
}

add_from_bone_marrow = function(t, C_n, c_n) {
  level_bcell = C_n*exp(-c_n*t)
  return(level_bcell)
}

lambda_activation = function(numTaxa, t, t_kickstart, naive_Bcell_count, antigens_sampled_uncoated, antigens_sampled_coated, baseline_activation, Ag_dependent_activation_c, Ag_dependent_activation_uc) {
  level_activation = rep(0, numTaxa)
  
  level_activation = baseline_activation+Ag_dependent_activation_uc*antigens_sampled_uncoated+Ag_dependent_activation_c*antigens_sampled_coated
  
  return(level_activation*naive_Bcell_count)
}

calculate_rates_onthefly_fast = function(numTaxa, t, selection_thresholds,
                                         proliferation_circulating_base, circulating_Bcell_count, plasma_Bcell_count, 
                                         circulating_affinity, plasma_affinity, influx_circulating, sample_size,
                                         rate_transfer_circulating_base, range_tf, range_apop, 
                                         tau_new, circulating_deviation, tau_c) {
  
  output = rep(0, 5*numTaxa)
  
  for (tx in 1:numTaxa) {
    
    affinity_threshold_circulating = selection_thresholds[tx]
    mean_aff_circulating           = circulating_affinity[tx]
    std_impact                     = tau_c*abs(affinity_threshold_circulating - mean_aff_circulating)
    new_circulating_from_naive     = influx_circulating[tx] # this is equal to rate of activation
    
    if(new_circulating_from_naive>1e-8){
      circulating_count = circulating_Bcell_count[tx]
      plasma_count      = plasma_Bcell_count[tx]
      
      std_aff_circulating  = circulating_deviation[tx]
      
      mean_aff_plasma      = plasma_affinity[tx]
      
      aff_th_diff      = affinity_threshold_circulating
      aff_th_diff_low  = max(0.01,(1-range_tf)*affinity_threshold_circulating)
      aff_th_diff_high = (1+range_tf)*affinity_threshold_circulating
      aff_th_apop_low  = max(0.01,(1-range_apop)*affinity_threshold_circulating)
      aff_th_apop_high = aff_th_diff_high
      
      mean_aff_circulating_darkzone = mean_aff_circulating 
      std_aff_circulating_darkzone  = std_aff_circulating
      
      circulating_darkzone_count    = circulating_count 
      mean_normal                   = mean_aff_circulating_darkzone
      sd_normal                     = std_aff_circulating_darkzone  
      
      # Emprical bounds
      upper_bound_empiric= min(1e5,qnorm(0.999, mean_normal, sd_normal)) # define a practical upperbound
      lower_bound_empiric= max(-1e5,qnorm(0.001, mean_normal, sd_normal)) # define a practical lowerbound
      
      # Calculate the region for transfer to plasma
      lb = aff_th_diff_low
      ub = aff_th_diff_high
      out_transfer   = calculate_expected_value_grid_region(mean_normal, sd_normal, lb, ub, sample_size)
      
      # Calculate the regions for apoptosis - too low affinity
      lb = lower_bound_empiric 
      ub = aff_th_apop_low
      out_apoptosis_low  = calculate_expected_value_grid_region(mean_normal, sd_normal, lb, ub, sample_size)
      
      # Calculate the regions for apoptosis - too high affinity
      lb = aff_th_apop_high
      ub = upper_bound_empiric
      out_apoptosis_high  = calculate_expected_value_grid_region(mean_normal, sd_normal, lb, ub, sample_size)
      
      # Calculate the regions for continue circulating
      lb = aff_th_apop_low
      ub = aff_th_diff_low
      out_circulating  = calculate_expected_value_grid_region(mean_normal, sd_normal, lb, ub, sample_size)
      
      # rate_transfer_circulating          = rate_transfer_circulating_base*(1/(log(750)))*log(1+max(0,(ti-t_first_plasma)))
      rate_transfer_circulating          = rate_transfer_circulating_base
      prob_transfer_circulating          = rate_transfer_circulating*out_transfer[1]
      avg_aff_transfer_circulating       = out_transfer[2]
      prob_apoptosis_low_circulating     = rate_transfer_circulating*out_apoptosis_low[1]
      avg_aff_apoptosis_low_circulating  = out_apoptosis_low[2]
      prob_apoptosis_high_circulating    = rate_transfer_circulating*out_apoptosis_high[1]
      avg_aff_apoptosis_high_circulating = out_apoptosis_high[2]
      
      prob_apoptosis_circulating    = prob_apoptosis_low_circulating+prob_apoptosis_high_circulating
      avg_aff_apoptosis_circulating = div0(prob_apoptosis_low_circulating*avg_aff_apoptosis_low_circulating+prob_apoptosis_high_circulating*avg_aff_apoptosis_high_circulating, prob_apoptosis_low_circulating+prob_apoptosis_high_circulating)
      
      prob_transfer_circulating2plasma = prob_transfer_circulating
      total_transfer2plasma            = prob_transfer_circulating2plasma*circulating_darkzone_count
      
      # this might make it faster at the beginning, but those will live short
      vacancy_mult                     = max(0,(1-plasma_count))
      prob_transfer_circulating2plasma = vacancy_mult*prob_transfer_circulating2plasma
      
      total_transfer2plasma            = prob_transfer_circulating2plasma*circulating_darkzone_count
      total_apoptosis                  = prob_apoptosis_circulating*circulating_darkzone_count
      prob_ct_circulating              = out_circulating[1]
      
      # first, apoptosis and transfer
      circulating_darkzone_count_after_apoptf = circulating_darkzone_count-total_apoptosis-total_transfer2plasma
      
      # then, proliferation
      # proliferation_circulating               = 1e3*proliferation_circulating_base
      proliferation_circulating               = proliferation_circulating_base
      total_proliferation                     = proliferation_circulating*circulating_darkzone_count_after_apoptf
      
      # delta counts - 
      net_delta_circulating                   = new_circulating_from_naive+total_proliferation-total_apoptosis-total_transfer2plasma 
      
      # new affinity before SHM -> mean doesn't change after proliferation, sample size changes
      new_total_affinity_after_proliferation  = out_circulating[2]*circulating_darkzone_count_after_apoptf*(1+proliferation_circulating)
      sd_normal_after_apoptf                  = (aff_th_diff_low-aff_th_apop_low)/6 # upperbound-lowerbound divided by 6, this can be more sophisticated later
      
      # then, SHM - mean affinity doesn't change, standart deviation changes
      new_total_affinity_after_SHM            = new_total_affinity_after_proliferation
      new_standart_dev_after_SHM              = sd_normal_after_apoptf+std_impact
      
      n_x  = circulating_darkzone_count_after_apoptf*(1+proliferation_circulating) 
      mu_x = out_circulating[2]
      std_x= new_standart_dev_after_SHM
      
      n_y  = new_circulating_from_naive
      mu_y = 0 # 0 mean
      std_y= tau_new*new_circulating_from_naive # more added, more the standart deviation
      
      # Calculate the combined mean
      mu_z = (n_x*mu_x+n_y*mu_y)/(n_x+n_y)
      # Calculate the combined std
      std_z = (n_x*(std_x+(mu_x-mu_z)^2)+n_y*(std_y+(mu_y-mu_z)^2))/(n_x+n_y)
      
      final_delta_aff_circulating = mu_z-mean_normal
      final_delta_std_circulating = std_z-sd_normal
      
      if(is.infinite(final_delta_std_circulating)){
        final_delta_std_circulating=1
      }
      
      #### PLASMA
      net_delta_plasma   = total_transfer2plasma
      new_total_affinity = mean_aff_plasma*plasma_count+avg_aff_transfer_circulating*total_transfer2plasma
      old_total_affinity = mean_aff_plasma*plasma_count
      if(is.nan(new_total_affinity)){
        new_total_affinity =0
      }
      final_delta_aff_plasma = div0(new_total_affinity, max(0,plasma_count+net_delta_plasma)) - mean_aff_plasma
      
      output[tx+0*numTaxa] = final_delta_aff_circulating
      output[tx+1*numTaxa] = final_delta_aff_plasma
      output[tx+2*numTaxa] = net_delta_circulating
      output[tx+3*numTaxa] = net_delta_plasma
      output[tx+4*numTaxa] = final_delta_std_circulating
    }
  }
  
  return(output)
}

function_eIgA = function(numTaxa, t, C_I, t_kickstart, number_of_plasma_cells) {
  return(C_I*th_0_vector(number_of_plasma_cells))
}


ODE_MODEL_STEROIDS = function(t, y, parms) { # TO UNDERSTAND WHEN TO OPEN M CELLS
  
  x_i   = parms$x_i
  x_r   = parms$x_r
  theta = parms$theta
  
  # Extract parameters using indices from 'x_i'
  numTaxa       = x_i[1]
  IgA_halflife  = x_r[3*numTaxa+1]
  t_kickstart   = x_r[3*numTaxa+3]
  EBF_duration  = x_r[3*numTaxa+9]
  MF_duration   = x_r[3*numTaxa+10]
  level_miga_in = x_r[3*numTaxa+11]
  level_hmo_in  = x_r[3*numTaxa+12]
  
  dydt        = numeric(1)
  level_s     = y[1]
  
  # Time-dependent parameter calculations
  level_HMO      = function_HMO(t, EBF_duration, MF_duration, level_hmo_in, 0)
  level_s_add    = function_steroids(t, level_miga_in)
  
  dydt[1] = level_HMO*level_s_add - IgA_halflife*level_s
  
  # Return the derivatives
  return(list(dydt))
}

ODE_MODEL = function(t, y, parms) {
  
  x_i   = parms$x_i
  x_r   = parms$x_r
  theta = parms$theta
  
  # Extract parameters using indices from 'x_i'
  numTaxa             = x_i[1]
  numAgnostic         = x_i[2]
  numGnostic          = x_i[3]
  sample_size         = x_i[3+numTaxa*numTaxa+1]
  last_index_theta_in = x_i[3+numTaxa*numTaxa+2]
  
  # Extract parameters using indices from 'x_r'
  O2Dependency_vector             = x_r[(0*numTaxa+1):(1*numTaxa)]
  TLR9_vector                     = x_r[(1*numTaxa+1):(2*numTaxa)]
  addWithFood_vector              = x_r[(2*numTaxa+1):(3*numTaxa)]
  IgA_halflife                    = x_r[3*numTaxa+1]
  bias_coating                    = x_r[3*numTaxa+2]
  t_kickstart                     = x_r[3*numTaxa+3]
  proliferation_circulating       = x_r[3*numTaxa+4]
  tau_new                         = x_r[3*numTaxa+5]
  tau_delta                       = x_r[3*numTaxa+6]
  tau_c                           = x_r[3*numTaxa+7]
  e_on                            = x_r[3*numTaxa+8]
  EBF_duration                    = x_r[3*numTaxa+9]
  MF_duration                     = x_r[3*numTaxa+10]
  level_miga                      = x_r[3*numTaxa+11]
  level_hmo                       = x_r[3*numTaxa+12]
  range_tf                        = x_r[3*numTaxa+13]
  range_apop                      = x_r[3*numTaxa+14]
  
  t0_mcell = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+9]
  
  # Initialize dydt
  dydt = numeric(19*numTaxa+7)
  
  abundancevectemp_uc_lumen      = th_0_vector(y[1:numTaxa])
  abundancevectemp_c_lumen       = th_0_vector(y[(numTaxa+1):(2*numTaxa)])
  abundancevectemp_k_feces_cumul = y[(2*numTaxa+1):(3*numTaxa)]
  
  naive_Bcell_count       = y[(3*numTaxa+1):(4*numTaxa)]
  circulating_Bcell_count = y[(4*numTaxa+1):(5*numTaxa)]
  plasma_Bcell_count      = y[(5*numTaxa+1):(6*numTaxa)]
  
  circulating_affinity = y[(6*numTaxa+1):(7*numTaxa)]
  circulating_std      = pmax(c(0,0,0,0),y[(7*numTaxa+1):(8*numTaxa)])
  plasma_affinity      = y[(8*numTaxa+1):(9*numTaxa)]
  level_eIgA           = y[(10*numTaxa+1):(11*numTaxa)]
  
  antigen_samples_uncoated_cumul = th_0_vector(y[(15*numTaxa+1):(16*numTaxa)])
  antigen_samples_coated_cumul   = th_0_vector(y[(16*numTaxa+1):(17*numTaxa)])
  
  selection_thresholds     = y[(18*numTaxa+1):(19*numTaxa)] 
  selection_thresholds_lag = rep(0,numTaxa)
  
  for(j in 1:numTaxa){
    selection_thresholds_lag[j] = ifelse(t-EBF_duration-1<0, 0, lagvalue(t-1, 18*numTaxa+j))
  }
  
  death_rate_plasma = pmin(rep(1,numTaxa),div0_vector_vector((selection_thresholds-selection_thresholds_lag),selection_thresholds))
  # print(c(t,death_rate_plasma))
  level_O2   = y[19*numTaxa+1]
  level_mIgA = y[19*numTaxa+2]
  level_s    = y[19*numTaxa+6]
  
  # Calculate anaerobe_abundance_total
  # Assuming y[1] and y[1+numTaxa] represent the E. shigella
  anaerobe_abundance_total = y[1]+y[1+numTaxa]
  
  # Ensuring levels are not negative
  if (level_O2 < 0) {
    level_O2 = 0
  }
  if (level_mIgA < 0) {
    level_mIgA = 0
  }
  
  level_Mcell = function_Mcell(t, EBF_duration, t_kickstart, t0_mcell, level_s)
  
  abundancevectemp_total_lumen   = abundancevectemp_uc_lumen+abundancevectemp_c_lumen;
  
  # Extracting sampled parameters using indices from 'theta'
  growthRate_vector               = theta[1:numTaxa]
  mIgA_reactivity_vector          = theta[(numTaxa+numAgnostic+numGnostic+1):(numTaxa+numAgnostic+numGnostic+numTaxa)]
  maternal_mIgA_reactivity_vector = mIgA_reactivity_vector[1:numTaxa]
  O2_degregation                  = theta[numTaxa+numAgnostic+numGnostic+numTaxa+1]
  interactionMat_flat             = theta[(numTaxa+numAgnostic+numGnostic+numTaxa+2):(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa)]
  ################# MODIFICATION FOR E AND B IMPACT  ################# 
  interactionMat                  = matrix(interactionMat_flat, ncol = numTaxa, byrow = TRUE)
  HMODependency_vector            = theta[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+1):(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+4)]
  solidDependency_vector          = theta[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+5):(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+8)]
  diassoc_rate_base               = theta[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+17)]
  kappa_vector                    = theta[(last_index_theta_in+1):(last_index_theta_in+numTaxa)]
  c_n                             = theta[last_index_theta_in+numTaxa+1]
  C_I                             = theta[last_index_theta_in+numTaxa+2]
  alpha_vector                    = theta[(last_index_theta_in+numTaxa+3):(last_index_theta_in+numTaxa+2+numTaxa)]
  sampling_rates_uc               = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+1]
  sampling_rates_c                = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+2]
  baseline_activation             = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+3]
  Ag_dependent_activation_uc      = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+4]
  Ag_dependent_activation_c       = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+5]
  rate_transfer_circulating       = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+6]
  baseline_transfer_circulating   = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+7]
  rate_transfer_plasma            = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+8]
  
  interaction_results     = interactionMat %*% abundancevectemp_total_lumen
  baseline_activation     = level_Mcell*baseline_activation  # it's zero, just to keep the idea here
  
  # Time-dependent parameter calculations
  level_HMO      = function_HMO(t, EBF_duration, MF_duration, level_hmo, 0)
  level_solid    = function_solid(t, EBF_duration, MF_duration, t_kickstart, level_hmo, 0)
  level_mIgA_add = function_mIgA(t, level_miga)
  level_s_add    = function_steroids(t, level_miga)
  level_eIgA_add = function_eIgA(numTaxa, t, C_I, t_kickstart, plasma_Bcell_count)
  
  # Dependencies and net growth rate
  O2Dependency    = O2Dependency_vector[1:numTaxa]*level_O2
  HMODependency   = HMODependency_vector[1:numTaxa]*level_HMO
  solidDependency = solidDependency_vector[1:numTaxa]*level_solid
  growthRate_net  = (rep(1, numTaxa)+O2Dependency+HMODependency+solidDependency)*growthRate_vector
  
  # Inflammation and antigen sampling
  total_inflammation      = sum(kappa_vector*abundancevectemp_uc_lumen)-sum(TLR9_vector*abundancevectemp_total_lumen)
  
  if(total_inflammation<0){total_inflammation=1e-5}
  
  antigens_sampled_uncoated      = rep(0,numTaxa)
  antigens_sampled_uncoated      = (1+total_inflammation)*alpha_vector*(level_Mcell*sampling_rates_uc)*abundancevectemp_uc_lumen
  antigens_sampled_coated        = (level_Mcell*sampling_rates_c)*abundancevectemp_c_lumen
  
  # Naive B cell dynamics
  influx_of_naive_Bcells = rep(add_from_bone_marrow(t, C_n, c_n), numTaxa) # check whether 100* makes sense (it does)
  rate_of_activation     = naive_Bcell_count*lambda_activation(numTaxa, t, t_kickstart, naive_Bcell_count, antigens_sampled_uncoated, antigens_sampled_coated, baseline_activation, Ag_dependent_activation_c, Ag_dependent_activation_uc)
  dy_naive_Bcell_count   = influx_of_naive_Bcells - rate_of_activation
  influx_circulating     = rate_of_activation
  
  antigen_samples_total_cumul= antigen_samples_coated_cumul+antigen_samples_uncoated_cumul
  adj_inf                    = kappa_vector[1]*(growthRate_vector[1]/abs(interactionMat_flat[1]))
  delta_selection_thresholds = rep(0,numTaxa)
  delta_selection_thresholds = level_Mcell*add_from_bone_marrow(t, 1, c_n)*tau_delta*log(1+total_inflammation/adj_inf)*log(1+div0_vector_vector(antigen_samples_uncoated_cumul,antigen_samples_total_cumul))
  
  output      = rep(0,4*numTaxa)
  output      = calculate_rates_onthefly_fast(numTaxa, t, selection_thresholds,
                                              proliferation_circulating, circulating_Bcell_count, plasma_Bcell_count, 
                                              circulating_affinity, plasma_affinity, influx_circulating, sample_size, 
                                              rate_transfer_circulating, range_tf, range_apop,
                                              tau_new, circulating_std, tau_c)
  
  dy_circulating_affinity    = output[(0*numTaxa+1):(1*numTaxa)]
  dy_plasma_affinity         = output[(1*numTaxa+1):(2*numTaxa)]
  dy_circulating_Bcell_count = output[(2*numTaxa+1):(3*numTaxa)] 
  dy_plasma_Bcell_count      = output[(3*numTaxa+1):(4*numTaxa)] - death_rate_plasma*plasma_Bcell_count # longevity depending on the T cell help
  dy_circulating_std         = output[(4*numTaxa+1):(5*numTaxa)]
  
  endogenous_mIgA_reactivity_vector = pmax(c(0,0,0,0),plasma_affinity)
  
  maternal_coating_in   = (maternal_mIgA_reactivity_vector>0)*(rep(1, numTaxa) / (bias_coating*maternal_mIgA_reactivity_vector+rep(1, numTaxa)))
  maternal_killing_in   = (maternal_mIgA_reactivity_vector>0)*(1-maternal_coating_in)
  endogenous_coating_in = (endogenous_mIgA_reactivity_vector>0)*(rep(1, numTaxa) / (bias_coating*endogenous_mIgA_reactivity_vector+rep(1, numTaxa)))
  endogenous_killing_in = (endogenous_mIgA_reactivity_vector>0)*(1-endogenous_coating_in)
  
  # IgA dynamics
  
  maternal_coating_in[which(maternal_coating_in<=0.5)]     = 0
  maternal_killing_in[which(maternal_killing_in<=0.5)]     = 0
  
  binding_ability_m = (maternal_mIgA_reactivity_vector>0)*(1-diassoc_rate_base/(1+maternal_mIgA_reactivity_vector));
  coating_eff_m     = (maternal_coating_in>0)*binding_ability_m;
  killing_eff_m     = (maternal_killing_in>0)*binding_ability_m;
  
  endogenous_coating_in[which(endogenous_coating_in<=0.5)] = 0
  endogenous_killing_in[which(endogenous_killing_in<=0.5)] = 0
  
  binding_ability_e = (endogenous_mIgA_reactivity_vector>0)*(1-diassoc_rate_base/(1+endogenous_mIgA_reactivity_vector));
  coating_eff_e     = (endogenous_coating_in>0)*binding_ability_e;
  killing_eff_e     = (endogenous_killing_in>0)*binding_ability_e;
  
  maternal_coating = level_mIgA*maternal_coating_in
  maternal_killing = level_mIgA*maternal_killing_in
  
  endogenous_coating = level_eIgA*endogenous_coating_in
  endogenous_killing = level_eIgA*endogenous_killing_in
  
  
  k               = (1-IgA_halflife/C_I)
  dydt_eIgA        = level_eIgA_add*(1-k*level_eIgA)-rep(1, numTaxa)*IgA_halflife*level_eIgA
  
  total_coating    = coating_eff_m*maternal_coating+e_on*coating_eff_e*endogenous_coating
  total_killing    = killing_eff_m*maternal_killing+e_on*killing_eff_e*endogenous_killing
  # Community dynamics
  dy_total          = abundancevectemp_total_lumen*(growthRate_net+interaction_results) - total_killing*abundancevectemp_uc_lumen
  dy_uncoated_lumen = dy_total+rep(IgA_halflife, numTaxa)*abundancevectemp_c_lumen - total_coating*abundancevectemp_uc_lumen
  dy_coated_lumen   = -rep(IgA_halflife, numTaxa)*abundancevectemp_c_lumen+total_coating*abundancevectemp_uc_lumen
  dy_killed_feces   = total_killing*abundancevectemp_uc_lumen
  
  # Assigning calculated derivatives to dydt vector
  dydt[(0*numTaxa+1):(1*numTaxa)]  = dy_uncoated_lumen
  dydt[(1*numTaxa+1):(2*numTaxa)]  = dy_coated_lumen
  dydt[(2*numTaxa+1):(3*numTaxa)]  = dy_killed_feces
  
  dydt[(3*numTaxa+1):(4*numTaxa)] = dy_naive_Bcell_count
  dydt[(4*numTaxa+1):(5*numTaxa)] = dy_circulating_Bcell_count
  dydt[(5*numTaxa+1):(6*numTaxa)] = dy_plasma_Bcell_count
  
  dydt[(6*numTaxa+1):(7*numTaxa)]  = dy_circulating_affinity
  dydt[(7*numTaxa+1):(8*numTaxa)]  = dy_circulating_std 
  # print(c(t,dy_circulating_std))
  
  dydt[(8*numTaxa+1):(9*numTaxa)]  = dy_plasma_affinity
  dydt[(9*numTaxa+1):(10*numTaxa)] = rep(0,numTaxa) #empty - use it to keeo track of stg you want!
  
  # IgA dynamics
  dydt[(10*numTaxa+1):(11*numTaxa)] = dydt_eIgA
  
  # IgA coating and killing dynamics
  dydt[(11*numTaxa+1):(12*numTaxa)] = maternal_coating
  dydt[(12*numTaxa+1):(13*numTaxa)] = maternal_killing
  dydt[(13*numTaxa+1):(14*numTaxa)] = endogenous_coating
  dydt[(14*numTaxa+1):(15*numTaxa)] = endogenous_killing
  
  # For observation
  dydt[(15*numTaxa+1):(16*numTaxa)] = antigens_sampled_uncoated
  dydt[(16*numTaxa+1):(17*numTaxa)] = antigens_sampled_coated
  dydt[(17*numTaxa+1):(18*numTaxa)] = rate_of_activation
  dydt[(18*numTaxa+1):(19*numTaxa)] = delta_selection_thresholds
  
  # Environmental and other components
  dydt[19*numTaxa+1] = -O2_degregation*level_O2*anaerobe_abundance_total
  dydt[19*numTaxa+2] = level_HMO*level_mIgA_add - IgA_halflife*level_mIgA
  dydt[19*numTaxa+3] = total_inflammation
  dydt[19*numTaxa+4] = level_HMO
  dydt[19*numTaxa+5] = level_solid
  dydt[19*numTaxa+6] = level_HMO*level_s_add - IgA_halflife*level_s
  
  # Return the derivatives
  return(list(dydt))
}

progressEnv <- new.env()
progressEnv$progress <- 0

# Define a wrapper function for ODE_MODEL
ODE_MODEL_wrapper = function(t, state, parms) {
  progressEnv$progress <- 0
  if (round(t) %% 100 == 0) {
    progressEnv$progress= round(t,0)
    setProgress(progressEnv$progress/840, detail = paste("Computed until day",progressEnv$progress))
    flush.console()
  }  
  return(ODE_MODEL(t, state, parms))
}
############### MODEL_BASE_FUNCTIONS.R - END ################################## 
# 
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      h1 {
        font-size: 22px; /* Adjust font size as needed */
      }
      #shiny-notification-panel {
        position: fixed;
        top: 50%;
        left: 70%;
        transform: translate(-50%, -50%);
        width: 50%;
      }
      .shiny-notification {
        width: 50%;
        border-radius: 0;
      }
    "))
  ),
  titlePanel(tags$h1("The ontogeny of immune tolerance: an integrated model of the early-life gut microbiome and adaptive immunity")),
  sidebarLayout(
    sidebarPanel(
      sliderInput("t_plot", "Plotting Interval (Days)", min = 5, max = 840, value = 840),
      sliderInput("EBF_duration", "Exclusive Breasfeeding Duration (Days)", min = 0, max = 154, value = 154),
      sliderInput("MF_duration", "Mixed Feeding Duration (Days)", min = 0, max = 308, value = 308),
      # actionButton("toggleButton1", "Endogenous SIgA Production: On"),
      actionButton("toggleButton2", "Hyperreactive SIgA Breastmilk: Off"),
      actionButton("toggleButton3", "SIgA Deficient Breastmilk: Off"),
      width = 3  # Adjust width of the sidebar
    ),
    mainPanel(
      fluidRow(plotOutput("odePlot2", height = "300px")),
      fluidRow(plotOutput("odePlot3", height = "300px")),
      fluidRow(plotOutput("odePlot1", height = "300px")),
      # plotOutput("odePlot2", height = "30%"),
      # plotOutput("odePlot3", height = "30%"),
      # plotOutput("odePlot1", height = "30%")
      width = 9  # Adjust width of the main panel
    )
  )
)

server <- function(input, output, session) {
  
  progress <- reactiveVal(0)  # Initialize progress at 0
  
  toggleValue1 = reactiveVal(1)  # Initially On
  toggleValue2 = reactiveVal(0)  # Initially Off
  toggleValue3 = reactiveVal(0)  # Initially Off
  
  # Observe event for the second toggle button
  observeEvent(input$toggleButton1, {
    new_val <- ifelse(toggleValue1() == 1, 0, 1)
    toggleValue1(new_val)
    updateActionButton(session, "toggleButton1", label = paste("Endogenous SIgA Production:", ifelse(new_val == 1, "On", "Off")))
  })
  
  # Observe event for the second toggle button
  observeEvent(input$toggleButton2, {
    new_val <- ifelse(toggleValue2() == 1, 0, 1)
    toggleValue2(new_val)
    updateActionButton(session, "toggleButton2", label = paste("Hyperreactive SIgA BM:", ifelse(new_val == 1, "On", "Off")))
  })
  
  # Observe event for the second toggle button
  observeEvent(input$toggleButton3, {
    new_val <- ifelse(toggleValue3() == 1, 0, 1)
    toggleValue3(new_val)
    updateActionButton(session, "toggleButton3", label = paste("SIgA Deficient BM:", ifelse(new_val == 1, "On", "Off")))
  })
  
  
  odeSolution <- reactive({
    
    withProgress(message = "Running...", value = 0, {
      
      add_pred          = 120
      e_on_in           = toggleValue1()
      ibd_switch        = toggleValue2()
      siga_def_switch   = toggleValue3()
      
      if(siga_def_switch==1){
        level_miga_in     = 0 # Scales the mSIgA levels, 0 to 1 (=0 mimics SIgA deficiency in the breastmilk)
      }else{
        level_miga_in     = 1 # Scales the mSIgA levels, 0 to 1 (=0 mimics SIgA deficiency in the breastmilk)
      }
      
      if(ibd_switch==1){
        mIgA_in           = rep(4.18,4) # allergic/IBD case, hyperreactive against symbiotic commensals (using the affinity value for E for all)
      }else{
        mIgA_in           = NA # mSIgA affinities. NA is the baseline scenario, where mIgA_in vector is c(4.18,0.18,0.13,0.12) for E, B, BC, and C
      }
      
      level_hmo_in      = 1 # Scales the HMOs levels, 0 to 1 (=1, HMOs level in the control case)
      mcell_th          = 0.5 # M cells open when mSIgA drops below (100*mcell_th)% of its maximum value. 
      EBF_duration      = input$EBF_duration # Exclusive Breastfeeding duration (days)
      MF_duration       = input$MF_duration # Mixed Feeding duration (days)
      
      ############# LOAD_DATA_LIST.R (cnt'd) ####################
      data_list_sim = list(
        
        numTaxa                = length(taxa_array),
        numTimeSteps           = length(days_array),
        numTimeSteps_pred      = length(days_array)+add_pred,
        numAgnostic            = numAgnostic,
        numGnostic             = numGnostic,
        last_index_theta_in    = length(theta_in),
        
        interactionMask_vector   = interactionMask_vector,
        y0_meanSubjects          = y0_meanSubjects_total, 
        observations             = tibble(abundanceArray_meanSubjects),
        
        O2Dependency_vector    = O2Dependency_vector,
        TLR4_vector            = TLR4_vector,
        TLR9_vector            = TLR9_vector,
        addWithFood_vector     = addWithFood_vector,
        
        IgA_halflife                   = 1/5, 
        bias_coating                   = 1,
        t_kickstart_min                = 30, # 1 month of silence for sure
        proliferation_circulating      = 1,
        sample_size                    = 10, # number of grid points for integration
        
        theta_in         = theta_in,
        srate_c_in       = srate_c_in,
        srate_uc_in      = srate_uc_in,
        srate_k_in       = srate_k_in,
        
        kappa_e  = kappa_vector[1],
        kappa_b  = kappa_vector[2], 
        kappa_bc = kappa_vector[3], 
        kappa_c  = kappa_vector[4],
        
        alpha_e  = alpha_vector[1],
        alpha_b  = alpha_vector[2],
        alpha_bc = alpha_vector[3],
        alpha_c  = alpha_vector[4],
        
        c_n                        = c_n, #ca
        range_tf                   = .25,
        range_apop                 = .75,
        rate_transfer_circulating  = 1,
        rate_transfer_plasma       = 1,
        C_n                        = C_n,
        C_I                        = C_I,
        Ag_dependent_activation_uc = 1,
        Ag_dependent_activation_c  = 1e-1,
        sampling_rates_uc_in       = 0.005,
        sampling_rates_c_in        = 0.05,
        solid_offset               = 60, 
        tau_new                    = tau_new,
        tau_delta                  = tau_delta,
        
        e_on_in                    = e_on_in,
        EBF_duration               = EBF_duration,
        EBF_duration_int           = EBF_duration,
        MF_duration                = MF_duration,
        level_hmo_in               = level_hmo_in,
        level_miga_in              = level_miga_in,
        tau_c                      = tau_c,
        
        t0         = day_start, #starting time
        t_data_int = as.integer(days_array),
        t_data     = days_array, #time bins of data
        ts_pred    = c(days_array,seq(max(days_array)+1,max(days_array)+add_pred)) #time bins of prediction (not doing prediction currently)
        
      )
      
      
      first_nonzero_values <- apply(data_list_sim$observations, 2, function(x) {
        # Find the index of the first nonzero value
        idx <- which(x > 1e-4)[1]
        # Return the first nonzero value using the index
        if (!is.na(idx)) {
          return(x[idx])
        } else {
          return(NA) # Return NA if there are no nonzero values
        }
      })
      
      abundanceArray_meanSubjects_day155 = abundanceArray_meanSubjects[155-3,]
      data_list_sim$add_e  = as.numeric(abundanceArray_meanSubjects_day155$Enterobacteriaceae)/function_solid(154+60,154,308,30,1,0)
      data_list_sim$add_bc = first_nonzero_values[3]/function_solid(154+60,154,308,30,1,0)
      data_list_sim$add_c  = first_nonzero_values[4]/function_solid(154+60,154,308,30,1,0)
      
      mIgA_in_keep = data_list_sim$theta_in[(numTaxa + numAgnostic + numGnostic + 1):(numTaxa + numAgnostic + numGnostic + numTaxa)]
      if(any(is.na(mIgA_in))){
        mIgA_in = mIgA_in_keep
      }
      
      mIgA_in_name = paste0(round(100*mIgA_in[1]),"_",round(100*mIgA_in[2]),"_",round(100*mIgA_in[3]),"_",round(100*mIgA_in[4]))
      data_list_sim$theta_in[(numTaxa + numAgnostic + numGnostic + 1):(numTaxa + numAgnostic + numGnostic + numTaxa)] = mIgA_in
      
      # Initialize hmo and miga
      hmo_0   = function_HMO(0, data_list_sim$EBF_duration, data_list_sim$MF_duration, data_list_sim$level_hmo_in, 0)
      miga_0  = hmo_0*function_mIgA(0, data_list_sim$level_miga_in)
      solid_0 = function_solid(0, data_list_sim$EBF_duration, data_list_sim$MF_duration, data_list_sim$t_kickstart, data_list_sim$level_hmo_in, 0)
      
      
      # Initialize x_r and x_i arrays
      x_i = integer(3+numTaxa*numTaxa+3)
      x_r = numeric(3*numTaxa+27)
      
      # Fill x_r array
      x_r[(0*numTaxa+1):(1*numTaxa)] = data_list_sim$O2Dependency_vector
      x_r[(1*numTaxa+1):(2*numTaxa)] = data_list_sim$TLR9_vector
      x_r[(2*numTaxa+1):(3*numTaxa)] = data_list_sim$addWithFood_vector
      x_r[(3*numTaxa+1):(3*numTaxa+14)] = c(data_list_sim$IgA_halflife, #1
                                            data_list_sim$bias_coating, #2
                                            data_list_sim$t_kickstart, #3
                                            data_list_sim$proliferation_circulating, #4
                                            data_list_sim$tau_new, #5
                                            data_list_sim$tau_delta,#6
                                            data_list_sim$tau_c,#7
                                            data_list_sim$e_on_in, #8
                                            data_list_sim$EBF_duration, #9
                                            data_list_sim$MF_duration, #10
                                            data_list_sim$level_miga_in, #11
                                            data_list_sim$level_hmo_in, #12
                                            data_list_sim$range_tf, #13
                                            data_list_sim$range_apop)#14
      
      # Fill x_i array
      x_i[1:3]                   = c(data_list_sim$numTaxa, data_list_sim$numAgnostic, data_list_sim$numGnostic)
      x_i[4:(3+numTaxa*numTaxa)] = data_list_sim$interactionMask_vector
      x_i[3+numTaxa*numTaxa+1]   = data_list_sim$sample_size
      x_i[3+numTaxa*numTaxa+2]   = data_list_sim$last_index_theta_in
      
      # Initialize theta vector
      theta = numeric(data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+8+6+1)
      
      # Assign initial values to theta
      theta[1:(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+8)] = data_list_sim$theta_in[1:(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+8)]
      theta[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+17)]  = data_list_sim$theta_in[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+17)] # diassociation rate
      
      drate_base             = data_list_sim$theta_in[(numTaxa+numAgnostic+numGnostic+numTaxa+1+numTaxa*numTaxa+17)]
      mIgA_reactivity_vector = theta[(numTaxa+numAgnostic+numGnostic+1):(numTaxa+numAgnostic+numGnostic+numTaxa)]
      
      init_coating    = as.numeric(mIgA_reactivity_vector>0)*(1/(data_list_sim$bias_coating*mIgA_reactivity_vector+1))
      init_killing    = as.numeric(mIgA_reactivity_vector>0)*(1-1/(data_list_sim$bias_coating*mIgA_reactivity_vector+1))
      coating_eff     = as.numeric(mIgA_reactivity_vector>0)*(1-drate_base/(mIgA_reactivity_vector+1))
      killing_eff     = as.numeric(mIgA_reactivity_vector>0)*(1-drate_base/(mIgA_reactivity_vector+1))
      
      init_coating[init_coating<=0.5]=0
      init_killing[init_killing<=0.5]=0
      coating_eff[init_coating<=0.5]=0
      killing_eff[init_killing<=0.5]=0
      
      init_coating = coating_eff*miga_0*init_coating
      init_killing = killing_eff*miga_0*init_killing
      init_noaction = 1-pmax(init_coating,init_killing)
      
      ### ADJUST BIFIDOBACTERIA FOR EBF DURATION
      k = 0.25 # assume if not breastfed at all, 75% decrease
      a = 1-k
      c = -(1/3)*log(1-(0.999-k)/a)
      m = k+a*(1-exp(-c*data_list_sim$level_hmo_in*data_list_sim$EBF_duration))
      mult = c(1,m,1,1)
      ### ADJUST BIFIDOBACTERIA FOR EBF DURATION
      
      y0 = numeric(19*numTaxa+7)
      y0[(0*numTaxa+1):(1*numTaxa)] = (1/data_list_sim$srate_uc_in)*init_noaction*mult*data_list_sim$y0_meanSubjects[1:numTaxa]
      y0[(1*numTaxa+1):(2*numTaxa)] = (1/data_list_sim$srate_c_in)*init_coating*mult*data_list_sim$y0_meanSubjects[1:numTaxa]
      y0[(2*numTaxa+1):(3*numTaxa)] = (1/data_list_sim$srate_k_in)*init_killing*mult*data_list_sim$y0_meanSubjects[1:numTaxa]
      
      y0[(3*numTaxa+1):(19*numTaxa)] = rep(0, 16*numTaxa)
      y0[(19*numTaxa+1)] = 1
      y0[(19*numTaxa+2)] = 0*miga_0
      y0[(19*numTaxa+3):(19*numTaxa+5)] = 0
      y0[(19*numTaxa+4)] = 0*hmo_0
      y0[(19*numTaxa+5)] = 0*solid_0
      y0[(19*numTaxa+6)] = 0*miga_0
      
      # Update theta with additional values
      theta[(data_list_sim$last_index_theta_in+1):(data_list_sim$last_index_theta_in+numTaxa)] = c(data_list_sim$kappa_e, data_list_sim$kappa_b, data_list_sim$kappa_bc, data_list_sim$kappa_c)
      theta[data_list_sim$last_index_theta_in+numTaxa+1] = data_list_sim$c_n 
      theta[data_list_sim$last_index_theta_in+numTaxa+2] = data_list_sim$C_I
      theta[(data_list_sim$last_index_theta_in+numTaxa+3):(data_list_sim$last_index_theta_in+numTaxa+2+numTaxa)] = c(data_list_sim$alpha_e, data_list_sim$alpha_b, data_list_sim$alpha_bc, data_list_sim$alpha_c)
      theta[(data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+1):(data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa)] = c(0,0,0,0) # empty - in case needed for stg else
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+1] = data_list_sim$sampling_rates_uc
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+2] = data_list_sim$sampling_rates_c
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+3] = 0 # data_list_sim$baseline_activation
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+4] = data_list_sim$Ag_dependent_activation_uc
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+5] = data_list_sim$Ag_dependent_activation_c
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+6] = data_list_sim$rate_transfer_circulating
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+7] = 0 #data_list_sim$baseline_transfer_circulating
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+8] = data_list_sim$rate_transfer_plasma
      theta[data_list_sim$last_index_theta_in+numTaxa+2+numTaxa+numTaxa+9] = 0 # assign 0 for now, this will be the M cell opening time calculated during integration
      
      # Define time vectors for ODE integration
      ts_pred   = seq(0,max(data_list_sim$ts_pred))
      ts_pred_1 = seq(0,data_list_sim$EBF_duration)
      ts_pred_2 = seq(data_list_sim$EBF_duration,max(data_list_sim$ts_pred))
      
      # Define column names here
      column_names = c("days",
                       paste0("uncoated_lumen_", 1:numTaxa),
                       paste0("coated_lumen_", 1:numTaxa),
                       paste0("killed_feces_cumul_", 1:numTaxa),
                       paste0("naive_Bcell_count_", 1:numTaxa),
                       paste0("circulating_Bcell_count_", 1:numTaxa),
                       paste0("plasma_Bcell_count_", 1:numTaxa),
                       paste0("circulating_affinity_", 1:numTaxa),
                       paste0("circulating_std_", 1:numTaxa),
                       paste0("plasma_affinity_", 1:numTaxa),
                       paste0("empty_", 1:numTaxa),
                       paste0("eIgA_", 1:numTaxa),
                       paste0("maternal_coating_", 1:numTaxa),
                       paste0("maternal_killing_", 1:numTaxa),
                       paste0("endogenous_coating_", 1:numTaxa),
                       paste0("endogenous_killing_", 1:numTaxa),
                       paste0("antigens_sampled_uncoated_", 1:numTaxa),
                       paste0("antigens_sampled_coated_", 1:numTaxa),
                       paste0("activated_Bcell_count_", 1:numTaxa),
                       paste0("selection_thresholds_", 1:numTaxa),
                       "level_O2",
                       "level_mIgA",
                       "level_inflammation_cumul",
                       "level_HMO_cumul",
                       "level_solid_cumul",
                       "level_steroids")
      
      mcell_th = 0.5 #
      parms = list(theta = theta, x_r = x_r, x_i = x_i)
      
      y0_s = miga_0
      
      
      # calculate time of m cell opening
      y_out_steroid = deSolve::ode(y = y0_s, times = ts_pred, func = ODE_MODEL_STEROIDS, parms = list(theta = theta, x_r = x_r, x_i = x_i), 
                                   atol = 1.0E-3,  rtol = 1.0E-1, maxsteps = 1000)
      y_out_steroid = as.data.frame(y_out_steroid)
      colnames(y_out_steroid)[2]='value'
      
      if(length(ts_pred_1)>1){
        if(is.na(which(y_out_steroid$value>mcell_th)[1])){ # never goes above mcell_th anyway
          # early m cell opening
          theta[data_list_sim$last_index_theta_in + numTaxa + 2 + numTaxa + numTaxa + 9] = data_list_sim$t_kickstart_min
          y0_in_2 = y0 # Assuming the last row of y_out_1 contains the state at EBF_duration_int, get rid of the first column which is days
        }else{
          y_out_steroid_above_05 = y_out_steroid[which(y_out_steroid$value>mcell_th)[1]:dim(y_out_steroid)[1],]
          if(which(y_out_steroid_above_05$value<mcell_th)[1]>EBF_duration){
            theta[data_list_sim$last_index_theta_in + numTaxa + 2 + numTaxa + numTaxa + 9] = EBF_duration
          }else{
            theta[data_list_sim$last_index_theta_in + numTaxa + 2 + numTaxa + numTaxa + 9] = which(y_out_steroid_above_05$value<mcell_th)[1]
          }
          if(theta[data_list_sim$last_index_theta_in + numTaxa + 2 + numTaxa + numTaxa + 9]<data_list_sim$t_kickstart_min){ #steroid can go below 0.5 before data_list_sim$t_kickstart_min
            theta[data_list_sim$last_index_theta_in + numTaxa + 2 + numTaxa + numTaxa + 9] = data_list_sim$t_kickstart_min
          }
        }
        
        # Run ODE solver for the first segment
        y_out_1 = deSolve::dede(y = y0, times = ts_pred_1, func = ODE_MODEL_wrapper, parms = list(theta = theta, x_r = x_r, x_i = x_i), 
                                atol = 1.0E-3,  rtol = 1.0E-1, maxsteps = 1000, control = list(mxhist = 1e5))
        y_out_1_df = as.data.frame(y_out_1)
        # Assign the column names
        names(y_out_1_df) = column_names
        
        y0_in_2 = y_out_1[nrow(y_out_1), -1] # Assuming the last row of y_out_1 contains the state at EBF_duration_int, get rid of the first column which is days
      }else{
        # early m cell opening
        theta[data_list_sim$last_index_theta_in + numTaxa + 2 + numTaxa + numTaxa + 9] = data_list_sim$t_kickstart_min
        y0_in_2 = y0 # Assuming the last row of y_out_1 contains the state at EBF_duration_int, get rid of the first column which is days
      }
      
      # Update initial state for the second segment
      
      t0_mcell   =  theta[data_list_sim$last_index_theta_in + numTaxa + 2 + numTaxa + numTaxa + 9] 
      solid_init = function_solid(data_list_sim$EBF_duration+data_list_sim$solid_offset, data_list_sim$EBF_duration, data_list_sim$MF_duration, data_list_sim$t_kickstart_min, data_list_sim$level_hmo_in, 0)
      
      # y0_in_2[0*numTaxa+1]  = solid_init*data_list_sim$add_e
      y0_in_2[0*numTaxa+3]  = solid_init*data_list_sim$add_bc
      y0_in_2[0*numTaxa+4]  = solid_init*data_list_sim$add_c
      
      
      y_out_2 = deSolve::dede(y = y0_in_2, times = ts_pred_2, func = ODE_MODEL_wrapper, parms = list(theta = theta, x_r = x_r, x_i = x_i),
                              atol = 1.0E-3,  rtol = 1.0E-1, maxsteps = 1000, control = list(mxhist = 1e5))
      
      y_out_2_df = as.data.frame(y_out_2)
      
      # Assign the column names
      names(y_out_2_df) = column_names
      
      # Combine results
      if(length(ts_pred_1)>1){
        y_out_df = rbind(y_out_1_df, y_out_2_df)
      }else{
        y_out_df = y_out_2_df
        
      }
      # delete the repeated row for the transition day
      y_out_df = y_out_df[-which(y_out_df$days==ts_pred_2[1])[1],]
      setProgress(1, message = "Complete")
      
      return(list(y_out_df=y_out_df, data_list_sim=data_list_sim, t0_mcell=t0_mcell))
    })
  })
  
  output$odePlot1 <- renderPlot({
    
    y_out_df      = odeSolution()$y_out_df
    data_list_sim = odeSolution()$data_list_sim
    t0_mcell      = odeSolution()$t0_mcell
    y_out_df      = as.data.frame(y_out_df)
    
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
    
    
    t_data_full    = seq(2+1,dim(abundanceArray_meanSubjects)[1]+2,1)
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
    total_abundance_scale = total_abundance_scale[seq(3,735,1)]
    
    total_abundance_scale = c(total_abundance_scale,rep(total_abundance_scale[length(total_abundance_scale)],dim(feces_daily_df_wide)[1]-length(total_abundance_scale)))
    
    feces_daily_df_wide_norm           = feces_daily_df_wide
    feces_daily_df_wide_norm_abundance = feces_daily_df_wide_norm[,2:5]
    feces_daily_df_wide_norm_abundance$Enterobacteriaceae =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Enterobacteriaceae
    feces_daily_df_wide_norm_abundance$Bifidobacteriaceae =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Bifidobacteriaceae
    feces_daily_df_wide_norm_abundance$Bacteroidaceae     =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Bacteroidaceae
    feces_daily_df_wide_norm_abundance$Clostridiales      =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Clostridiales
    
    
    feces_daily_df_wide_norm = cbind(feces_daily_df_wide[,1],feces_daily_df_wide_norm_abundance)
    
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
    combined_data_norm = combined_data_norm %>% dplyr::filter(days<=input$t_plot)
    
    p_rel= ggplot(combined_data_norm, aes(x = days, y = value, color = taxa, linetype = type)) +
      geom_line(size = 0.8) +  theme_minimal() +
      labs(title = "Relative abundances in fecal samples", x = "DOL", y = "Ratio") +
      scale_color_manual(values =my_colors) +
      scale_linetype_manual(values = c("Model" = "11", "Observation" = "solid")) +
      theme(legend.title = element_blank(), plot.title = element_text(size = 12), legend.text = element_text(size = 11)) 
    
    combined_data_norm     = as.data.frame(combined_data_norm)
    
    lumen_daily_c_df  = y_lumen[c('days','taxa','median_c')]
    lumen_daily_uc_df = y_lumen[c('days','taxa','median_uc')]
    
    lumen_daily_c_df$coating  = 'SIgA+'
    lumen_daily_uc_df$coating = 'SIgA-'
    colnames(lumen_daily_c_df)[3]='median'
    colnames(lumen_daily_uc_df)[3]='median'
    
    lumen_daily_cuc_df = rbind(lumen_daily_c_df,lumen_daily_uc_df)
    lumen_daily_cuc_df$median=as.numeric(lumen_daily_cuc_df$median)
    
    seq_ticks = seq(0, 0.15, 0.05)
    y_up = round(max(lumen_daily_cuc_df$median)*1.1,2)
    lumen_daily_cuc_df = as.data.frame(lumen_daily_cuc_df)
    lumen_daily_cuc_df = lumen_daily_cuc_df %>% dplyr::filter(days<=input$t_plot)
    
    p_abs=ggplot(lumen_daily_cuc_df, aes(x = days, y = median, color = taxa, linetype = coating)) +
      geom_line(size = 0.8) +  theme_minimal() +
      labs(title = "Absolute abundances in gut lumen", x = "DOL", y = "(Cells/g LC) x 1e-11") +
      scale_color_manual(values = my_colors) +
      scale_linetype_manual(values = c("SIgA+" = "11", "SIgA-" = "solid")) +
      theme(legend.title = element_blank(), plot.title = element_text(size = 12), legend.text = element_text(size = 11))+
      scale_y_continuous(breaks = seq_ticks) + ylim(0,y_up)
    
    plot_grid(p_rel,p_abs,ncol = 2, rel_widths = c(1,1))
    
  })
  
  output$odePlot2 <- renderPlot({
    
    y_out_df      = odeSolution()$y_out_df
    data_list_sim = odeSolution()$data_list_sim
    t0_mcell      = odeSolution()$t0_mcell
    y_out_df      = as.data.frame(y_out_df)
    
    library(gridExtra)
    library(ggpattern) 
    library(patchwork)
    library(png)
    library(grid)
    my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
    
    #### CALORIES
    times = y_out_df$days
    
    # O2
    variable       = 'level_O2'
    y_out_df_use   = y_out_df[c('days',variable)]
    colnames(y_out_df_use)[2] = 'median'
    y_O2 = y_out_df_use
    O2_values = y_O2$median
    
    # HMO
    variable       = 'level_HMO_cumul'
    y_out_df_use   = y_out_df[c('days',variable)]
    y_out_df_use_1 = y_out_df_use[1,]
    # Apply differential calculation to each column except 'days'
    for (col in names(y_out_df_use)[-1]) {  # Skipping the first column ('days')
      y_out_df_use[[col]] = c(NA, diff(y_out_df_use[[col]]))
    }
    y_out_df_use[1,] = y_out_df_use_1
    colnames(y_out_df_use)[2] = 'median'
    y_out_df_use$median = y_out_df_use$median/unique(diff(y_out_df_use$days)) # proper values for differentiation
    y_HMO = y_out_df_use
    HMO_values = y_HMO$median
    
    # SOLID
    variable       = 'level_solid_cumul'
    y_out_df_use   = y_out_df[c('days',variable)]
    y_out_df_use_1 = y_out_df_use[1,]
    # Apply differential calculation to each column except 'days'
    for (col in names(y_out_df_use)[-1]) {  # Skipping the first column ('days')
      y_out_df_use[[col]] = c(NA, diff(y_out_df_use[[col]]))
    }
    y_out_df_use[1,] = y_out_df_use_1
    colnames(y_out_df_use)[2] = 'median'
    y_out_df_use$median = y_out_df_use$median/unique(diff(y_out_df_use$days)) # proper values for differentiation
    y_solid = y_out_df_use
    solid_values = y_solid$median
    
    
    # mSIgA
    variable       = 'level_mIgA'
    y_out_df_use   = y_out_df[c('days',variable)]
    colnames(y_out_df_use)[2] = 'median'
    y_mSIgA = y_out_df_use
    mSIgA_values = y_mSIgA$median
    
    
    EBF_duration = data_list_sim$EBF_duration
    MF_duration  = data_list_sim$MF_duration
    
    # mcell
    Mcell_values = c()
    if(min(times)==0){
      for(t in times){ # t0_mcell comes from integration
        Mcell_values[t+1] = function_Mcell(t, EBF_duration,  data_list_sim$t_kickstart_min, t0_mcell, 1)
      }
    }else{
      for(t in times){ # t0_mcell comes from integration
        Mcell_values[t] = function_Mcell(t, EBF_duration,  data_list_sim$t_kickstart_min, t0_mcell, 1)
      }
    }
    
    # Create a data frame for plotting
    df_cal = data.frame(
      t = times,
      HMO = HMO_values,
      Mcell = Mcell_values,
      Solid = solid_values,
      mSIgA = mSIgA_values
    )
    
    df_cal$t1 = NA
    df_cal$t2 = NA
    
    if(EBF_duration>0){
      df_cal[EBF_duration+0,]$t1 = 0
      df_cal[EBF_duration+1,]$t1 = 1.20
    }
    
    solid_cutoff = EBF_duration + MF_duration
    
    if(solid_cutoff>0){
      df_cal[solid_cutoff+0,]$t2 = 0
      df_cal[solid_cutoff+1,]$t2 = 1.20
    }
    
    # Reshaping the data
    colnames(df_cal)[2]='HMOs calories (x 1e-3)'
    colnames(df_cal)[3]='Endogenous activation'
    colnames(df_cal)[4]='PDPs calories (x 1e-3)'
    colnames(df_cal)[5]='mSIgA levels (lumen)'
    colnames(df_cal)[6]='Transition to MF'
    colnames(df_cal)[7]='Transition to ESF'
    
    df_cal[(EBF_duration+MF_duration+1):dim(df_cal)[1],5] = 0
    
    df_long = pivot_longer(df_cal, cols = -t, names_to = "series", values_to = "value")
    
    # Define colors and linetypes
    colors = c("HMOs calories (x 1e-3)" = "#FF006E", "Endogenous activation" = "gray50", "PDPs calories (x 1e-3)" = "#12db99", 
               "mSIgA levels (lumen)" = "#020887",
               "Transition to MF" = "#FF006E", "Transition to ESF" = "#12db99")
    
    linetypes =  c("HMOs calories (x 1e-3)" = "solid", "Endogenous activation" = "11", "PDPs calories (x 1e-3)" = "solid", 
                   "mSIgA levels (lumen)" = "dotdash",
                   "Transition to MF" = "11", "Transition to ESF" = "11")
    
    df_long$series <- factor(df_long$series, levels = c(
      "Endogenous activation",
      "HMOs calories (x 1e-3)",
      "PDPs calories (x 1e-3)",
      "Transition to MF",
      "Transition to ESF",
      "mSIgA levels (lumen)"
    ))
    
    df_long = as.data.frame(df_long)
    df_long = df_long %>% dplyr::filter(t<=input$t_plot)
    
    p_cals=ggplot(df_long, aes(x = t, y = value, color = series, linetype = series)) +
      geom_line(data = df_long, size=0.7) +
      theme_minimal() +
      scale_color_manual(values = colors) +
      scale_linetype_manual(values = linetypes) +
      labs(title ='Pre-computed model inputs', x = "DOL", y = "") +
      theme(legend.title = element_blank(), plot.title = element_text(size = 12)) +
      guides(size = "none") + theme(legend.title = element_blank(),
                                    legend.text = element_text(size = 11),
                                    # axis.title = element_text(size = 12),
                                    # axis.text = element_text(size = 12),
                                    # legend.position = c(0.95, 0.074), # This needs adjustment
                                    # legend.justification = c("right", "bottom"), # Anchor point
                                    # legend.box.just = "right",
                                    # legend.margin = margin(b = 0, t=-4), # Increased bottom margin for padding
                                    legend.background = element_rect(colour = "gray", fill="white")) 
    
    # Adding lines with arrows on both ends
    yy = 1.25
    # Adding labels for the time periods
    yy_off = 0.15
    if(EBF_duration>0 & MF_duration>0){
      p_cals <- p_cals + 
        annotate("segment", x = 0, xend = EBF_duration, y = yy, yend = yy, colour = "#FF006E", 
                 arrow = arrow(type = "closed", ends = "both",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf) +
        annotate("segment", x = EBF_duration, xend = EBF_duration+MF_duration, y = yy, yend = yy, colour = "#12db99", 
                 arrow = arrow(type = "closed", ends = "both",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf) + 
        annotate("segment", x = EBF_duration+MF_duration, xend = max(df_long$t), y = yy, yend = yy, colour = "#1f6750", 
                 arrow = arrow(type = "closed", ends = "first",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf)
      p_cals <- p_cals + 
        annotate("text", x = EBF_duration/2, y = yy+yy_off, label = "EBF", vjust = 2, size = 4) +
        annotate("text", x = EBF_duration+MF_duration/2, y = yy+yy_off, label = "MF", vjust = 2, size = 4) + 
        annotate("text", x = EBF_duration+MF_duration+(max(df_long$t)-(EBF_duration+MF_duration))/2, y = yy+yy_off, label = "ECF", vjust = 2, size = 4) + ylim(0,1.4)
    }else if(EBF_duration>0 & MF_duration==0){
      p_cals <- p_cals + 
        annotate("segment", x = 0, xend = EBF_duration, y = yy, yend = yy, colour = "#FF006E", 
                 arrow = arrow(type = "closed", ends = "both",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf) +
        annotate("segment", x = EBF_duration+MF_duration, xend = max(df_long$t), y = yy, yend = yy, colour = "#1f6750", 
                 arrow = arrow(type = "closed", ends = "first",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf)
      p_cals <- p_cals + 
        annotate("text", x = EBF_duration/2, y = yy+yy_off, label = "EBF", vjust = 2, size = 4) +
        annotate("text", x = EBF_duration+MF_duration+(max(df_long$t)-(EBF_duration+MF_duration))/2, y = yy+yy_off, label = "ECF", vjust = 2, size = 4) + ylim(0,1.4)
    }else if(EBF_duration==0 & MF_duration>0){
      p_cals <- p_cals + 
        annotate("segment", x = EBF_duration, xend = EBF_duration+MF_duration, y = yy, yend = yy, colour = "#12db99", 
                 arrow = arrow(type = "closed", ends = "both",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf) + 
        annotate("segment", x = EBF_duration+MF_duration, xend = max(df_long$t), y = yy, yend = yy, colour = "#1f6750", 
                 arrow = arrow(type = "closed", ends = "first",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf)
      p_cals <- p_cals + 
        annotate("text", x = EBF_duration+MF_duration/2, y = yy+yy_off, label = "MF", vjust = 2, size = 4) + 
        annotate("text", x = EBF_duration+MF_duration+(max(df_long$t)-(EBF_duration+MF_duration))/2, y = yy+yy_off, label = "ECF", vjust = 2, size = 4) + ylim(0,1.4)
    }else if(EBF_duration==0 & MF_duration==0){
      p_cals <- p_cals + 
        annotate("segment", x = EBF_duration+MF_duration, xend = max(df_long$t), y = yy, yend = yy, colour = "#1f6750", 
                 arrow = arrow(type = "closed", ends = "first",  length = unit(0.15, "inches")), size = 0.5, yintercept = -Inf)
      p_cals <- p_cals + 
        annotate("text", x = EBF_duration+MF_duration+(max(df_long$t)-(EBF_duration+MF_duration))/2, y = yy+yy_off, label = "ECF", vjust = 2, size = 4) + ylim(0,1.4)
    }
    
    p_cals = p_cals + scale_x_continuous(limits = c(0, input$t_plot))
    
    
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
    
    y_out_df_use$median = y_out_df_use$median/max(y_out_df_use$median)
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
    
    p_micro <- ggplot(df_long, aes(x = t, y = value, color = series, linetype = series)) +
      geom_line(data = df_long, size=0.7) +
      theme_minimal() +
      scale_color_manual(values = colors, labels = c("O2 concentration", "Microenvironmental\nstimulation")) +
      scale_linetype_manual(values = linetypes, labels = c("O2 concentration", "Microenvironmental\nstimulation")) +
      labs(title = 'Dynamic variables, gut lumen', x = "DOL", y = "") +
      theme(legend.title = element_blank(), plot.title = element_text(size = 12)) +
      guides(size = "none") + 
      theme(legend.title = element_blank(),
            legend.text = element_text(size = 11),
            # axis.title = element_text(size = 12),
            # axis.text = element_text(size = 12),
            # legend.position = c(1, 0.65),
            # legend.justification = c("right", "bottom"), # Anchor point
            # legend.box.just = "right",
            # legend.margin = margin(b = 4, t=-4), # Increased bottom margin for padding
            legend.background = element_rect(colour = "gray", fill="white")) 
    
    p_micro = p_micro + scale_x_continuous(limits = c(0, input$t_plot))
    
    plot_grid(p_cals,p_micro,ncol = 2,rel_widths = c(1.2,1))
    
  })
  
  output$odePlot3 <- renderPlot({
    
    y_out_df      = odeSolution()$y_out_df
    data_list_sim = odeSolution()$data_list_sim
    t0_mcell      = odeSolution()$t0_mcell
    y_out_df      = as.data.frame(y_out_df)
    
    library(gridExtra)
    library(ggpattern) 
    library(patchwork)
    library(png)
    library(grid)
    my_colors = c('#2bc5d8','#074fa8','#ff9aa1','#f20026')
    
    #### PLASMA AFFINITY
    keep_miga_vec  = c(4.18,0.18,0.13,0.12)
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
    out_df_ene_long = as.data.frame(out_df_ene_long)
    out_df_ene_long = out_df_ene_long %>% dplyr::filter(days<=input$t_plot)
    
    ##Plot
    p_aff_log <- ggplot(out_df_ene_long, aes(x = days, y = value, color = taxa, linetype = type)) +
      geom_line(size = 1) +
      theme_minimal() +
      scale_linetype_manual(values = c("median" = "solid", "miga" = "dotted"),
                            labels = c("median" ="eSIgA", "miga"="mSIgA")) +
      labs(title = "Average eSIgA Affinity", x = "DOL", y = "") +
      scale_color_manual(values = my_colors[1:4]) +  # Ensure my_colors is defined
      theme(legend.title = element_blank(), plot.title = element_text(size = 12), legend.text = element_text(size = 11)) +
      scale_y_log10(minor_breaks = seq(0, 6, .2) * c(1e-2, 1e-1, 1e0))
    
    p_aff_log = p_aff_log + theme(legend.title = element_blank(),
                                  # legend.position = c(1, 0.13),
                                  # legend.justification = c("right", "bottom"),
                                  # legend.box.just = "right",
                                  # legend.margin = margin(t = -1, b = 0),
                                  legend.box.background = element_rect(colour = NA, fill = NA),
                                  plot.background = element_rect(fill = NA, colour = NA),
                                  panel.background = element_rect(fill = NA, colour = NA))
    
    # If you need to manually adjust the order of the legends
    p_aff_log = p_aff_log + guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))
    
    ##################
    out_df_ene_long_pick = out_df_ene_long %>% filter(type=='median' & value>0)
    p_aff_log = p_aff_log + coord_cartesian(xlim = c(min(out_df_ene_long_pick$days), input$t_plot), ylim = c(0.01, mv))
    # p_aff_log = p_aff_log + scale_x_continuous(limits = c(0, input$t_plot))+ scale_y_continuous(limits = c(0.01, mv))
    
    colnames(abundanceArray_meanSubjects) = taxa_array
    abundanceArray_meanSubjects <- as.data.frame(abundanceArray_meanSubjects) %>%
      dplyr::mutate(across(c(Enterobacteriaceae, Bifidobacteriaceae, Bacteroidaceae, Clostridiales),
                           ~pmax(., 0, na.rm = TRUE)))
    
    
    t_data_full    = seq(2+1,dim(abundanceArray_meanSubjects)[1]+2,1)
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
    y_out_df_use   = y_out_df[c('days',paste0(variable,'_',seq(1,numTaxa,1)))]
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
    
    srate_coated   = 0.6975581 # data_list_sim$srate_c_in
    srate_uncoated = 0.8035475 #data_list_sim$srate_uc_in
    srate_killed   = 0.8286261 #data_list_sim$srate_k_in
    
    y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_tot = srate_coated*median_c+srate_killed*median_k+srate_uncoated*median_uc)
    y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_ck_tot = srate_coated*median_c+srate_killed*median_k)
    y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_c_tot  = srate_coated*median_c)
    y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_k_tot  = srate_killed*median_k)
    y_lumen = y_lumen %>% dplyr::rowwise() %>% dplyr::mutate(median_uc_tot = srate_uncoated*median_uc)
    
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
    total_abundance_scale = total_abundance_scale[seq(3,735,1)]
    
    total_abundance_scale = c(total_abundance_scale,rep(total_abundance_scale[length(total_abundance_scale)],dim(feces_daily_df_wide)[1]-length(total_abundance_scale)))
    
    feces_daily_df_wide_norm           = feces_daily_df_wide
    feces_daily_df_wide_norm_abundance = feces_daily_df_wide_norm[,2:5]
    feces_daily_df_wide_norm_abundance$Enterobacteriaceae =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Enterobacteriaceae
    feces_daily_df_wide_norm_abundance$Bifidobacteriaceae =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Bifidobacteriaceae
    feces_daily_df_wide_norm_abundance$Bacteroidaceae     =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Bacteroidaceae
    feces_daily_df_wide_norm_abundance$Clostridiales      =  (0.01*total_abundance_scale/rowSums(feces_daily_df_wide[,2:5]))*feces_daily_df_wide_norm_abundance$Clostridiales
    
    
    feces_daily_df_wide_norm = cbind(feces_daily_df_wide[,1],feces_daily_df_wide_norm_abundance)
    
    observations_wide                = cbind(t_data_full,abundanceArray_meanSubjects)
    observations_wide_norm_abundance = observations_wide
    colnames(observations_wide_norm_abundance)[1]= 'time'
    
    feces_daily_df_norm    = feces_daily_df_wide_norm %>% pivot_longer(!days, names_to = "taxa", values_to = "value")
    observations_long_norm = observations_wide_norm_abundance %>% pivot_longer(!time, names_to = "taxa", values_to = "value")
    
    feces_daily_df_norm$type     <- "Model"
    observations_long_norm$type <- "Observation"
    colnames(observations_long_norm)[1]='days'
    # Combine the datasets
    combined_data_norm <- rbind(feces_daily_df_norm, observations_long_norm)
    # Plotting
    
    combined_data_norm     = as.data.frame(combined_data_norm)
    lumen_daily_c_df  = y_lumen[c('days','taxa','median_c')]
    lumen_daily_uc_df = y_lumen[c('days','taxa','median_uc')]
    
    lumen_daily_c_df$coating  = 'SIgA+'
    lumen_daily_uc_df$coating = 'SIgA-'
    colnames(lumen_daily_c_df)[3]='median'
    colnames(lumen_daily_uc_df)[3]='median'
    
    lumen_daily_cuc_df = rbind(lumen_daily_c_df,lumen_daily_uc_df)
    lumen_daily_cuc_df$median=as.numeric(lumen_daily_cuc_df$median)
    
    seq_ticks = seq(0, 0.15, 0.05)
    y_up = round(max(lumen_daily_cuc_df$median)*1.1,2)
    
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
    
    diassoc_rate_base = 0.7097489
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
    df_feces = df_feces %>% dplyr::rowwise() %>% dplyr::mutate(iga_index=-1*(log(IgA_plus)-log(IgA_minus))/(log(IgA_plus)+log(IgA_minus)))
    
    
    df_feces$sample = 'feces'
    
    df_iga_index     = df_feces[c('days','taxa','iga_index')]
    df_coating_ratio = df_feces[c('days','taxa','IgA_coated_fraction')]
    df_coating_ratio_c = df_feces[c('days','taxa','IgA_coated_fraction_c')]
    df_coating_ratio_k = df_feces[c('days','taxa','IgA_coated_fraction_k')]

    
    df_iga_index_conv       = df_iga_index %>% filter(days==input$t_plot) # was 735 in the main plot in paper
    df_coating_ratio_c_conv = df_coating_ratio_c  %>% filter(days==input$t_plot)
    df_coating_ratio_k_conv = df_coating_ratio_k  %>% filter(days==input$t_plot)
    
    df_coating_ratio_conv = merge(df_coating_ratio_c_conv,df_coating_ratio_k_conv, by=c('days','taxa'))
    df_coating_ratio_conv = merge(df_coating_ratio_conv,df_iga_index_conv, by=c('days','taxa'))
    
    df_long <- pivot_longer(df_coating_ratio_conv, cols = starts_with("IgA_coated_fraction_") | starts_with("iga_index"),
                            names_to = "Measurement", values_to = "Value")
    
    df_long$taxa = factor(df_long$taxa, levels = c("Clostridiales", "Bacteroidaceae",  "Bifidobacteriaceae","Enterobacteriaceae"))
    
    
    z = structure(list(Value= df_long$Value, Taxa= df_long$taxa, Measurement = df_long$Measurement), row.names = c(NA, -12L), class = c("tbl_df", "tbl", "data.frame"))
    
    # Create the combined plot with corrections
    p_coating = ggplot(z, aes(x = Taxa, y = Value, fill = Measurement, pattern = Measurement)) +
      geom_bar_pattern(
        stat = "identity", position = "dodge",
        pattern_spacing = 0.05,
        pattern_angle = 45
      ) +
      coord_flip() +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Corrected line at y=0
      scale_fill_manual(values=c('gray', 'black','#FC6736'),labels = c("SIgA+ (M)", "SIgA+ (N)", "IgA Index")) +
      scale_pattern_manual(values=c('none', 'none', 'stripe'),labels = c("SIgA+ (M)", "SIgA+ (N)", "IgA Index")) +
      ggpubr::theme_pubr() +
      theme(
        legend.position = "right",
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
        plot.title = element_text(size = 12)
      ) +
      labs(title  = "Coating fractions and IgA Index", y = "Value") +
      scale_y_continuous(limits = c(-0.7, 0.7))
    
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
    
    # Plotting
    
    p_eiga = ggplot(y_out_df_use_l, aes(x = days, y = median, color = taxa)) +
      geom_line(size = 1) +  theme_minimal() +
      labs(title = "eSIgA Concentration", x = "DOL", y = "") +
      scale_color_manual(values = my_colors) +  
      theme(legend.title = element_blank(), plot.title = element_text(size = 12)) 
    
    p_eiga = p_eiga + theme(legend.title = element_blank(),
                            legend.text = element_text(size = 11),
                            # axis.title = element_text(size = 10),
                            # axis.text = element_text(size = 10),
                            # legend.position = c(0.28, 0.62), # This needs adjustment
                            # legend.justification = c("right", "bottom"), # Anchor point
                            # legend.box.just = "right",
                            # legend.margin = margin(),
                            legend.box.background = element_rect(colour = NA, fill = NA), # Make legend background transparent
                            # Optional: Make the entire plot background transparent
                            plot.background = element_rect(fill = NA, colour = NA),
                            panel.background = element_rect(fill = NA, colour = NA))
    p_eiga = p_eiga + scale_x_continuous(limits = c(0, input$t_plot))
    plot_grid(p_aff_log,p_eiga,p_coating,ncol = 3, rel_widths = c(1,1,1))
  })
  
}

shinyApp(ui = ui, server = server)

