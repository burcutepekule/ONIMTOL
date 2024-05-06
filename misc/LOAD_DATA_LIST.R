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