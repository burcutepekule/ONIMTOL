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
    diff_aff_std                   = round((affinity_threshold_circulating - mean_aff_circulating),6)
    std_impact                     = tau_c*(diff_aff_std^2)
    
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
      
      # # this might make it faster at the beginning, but those will live short
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
      # new_standart_dev_after_SHM              = sd_normal_after_apoptf+std_impact
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
  
  # y = round(y,8) # ADDED 18TH OF JUNE
  # y[y < 1e-8] = 0  # Or another small threshold
  
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
  abundancevectemp_total_lumen_cumul = th_0_vector(y[(9*numTaxa+1):(10*numTaxa)])
  
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
    # print(c(t,t-1))
    selection_thresholds_lag[j] = ifelse(t-EBF_duration-1<0, 0, lagvalue(t-1, 18*numTaxa+j))
  }
  
  death_rate_plasma = pmin(rep(1,numTaxa),div0_vector_vector((selection_thresholds-selection_thresholds_lag),selection_thresholds))
  
  # death_rate_plasma = c(0,0,0,0)
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
  C_n                             = theta[last_index_theta_in+numTaxa+2+numTaxa+numTaxa+10]
  
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
  
  if(total_inflammation<0){total_inflammation=0}
  
  antigens_sampled_uncoated      = rep(0,numTaxa)
  antigens_sampled_uncoated      = (1+total_inflammation)*alpha_vector*(level_Mcell*sampling_rates_uc)*abundancevectemp_uc_lumen
  antigens_sampled_coated        = (level_Mcell*sampling_rates_c)*abundancevectemp_c_lumen
  antigens_sampled_total         = antigens_sampled_uncoated+antigens_sampled_coated
  
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
  
  binding_ability_m = (maternal_mIgA_reactivity_vector>0)*(1-diassoc_rate_base/(1+maternal_mIgA_reactivity_vector));
  coating_eff_m     = binding_ability_m;
  killing_eff_m     = binding_ability_m;
  
  binding_ability_e = (endogenous_mIgA_reactivity_vector>0)*(1-diassoc_rate_base/(1+endogenous_mIgA_reactivity_vector));
  coating_eff_e     = binding_ability_e;
  killing_eff_e     = binding_ability_e;
  
  
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
  dydt[(9*numTaxa+1):(10*numTaxa)] = abundancevectemp_total_lumen 
  
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
  
  return(list(dydt))
}


# Define a wrapper function for ODE_MODEL
ODE_MODEL_wrapper = function(t, state, parms) {
  if (round(t) %% 10 == 0) {
    print(paste("Solving for day:", round(t,3)))
    flush.console()
  }
  
  return(ODE_MODEL(t, state, parms))
}
