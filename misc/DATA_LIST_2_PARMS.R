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

