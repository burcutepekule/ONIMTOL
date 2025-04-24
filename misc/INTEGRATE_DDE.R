parms = list(theta = theta, x_r = x_r, x_i = x_i)

y0_s = miga_0


# calculate time of m cell opening
y_out_steroid = deSolve::ode(y = y0_s, times = ts_pred, func = ODE_MODEL_STEROIDS, parms = list(theta = theta, x_r = x_r, x_i = x_i), 
                             atol = 1.0E-8,  rtol = 1.0E-5, maxsteps = 1000)
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
  y_out_1 = deSolve::dede(y = y0, times = ts_pred_1, func = ODE_MODEL, parms = list(theta = theta, x_r = x_r, x_i = x_i), 
                          atol = 1.0E-8,  rtol = 1.0E-5, maxsteps = 1000, control = list(mxhist = 1e5))
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

if(!exists('silent')){
  silent=0 # print
}

if(silent==1){
  # Run ODE solver for the second segment without the wrapper function
  y_out_2 = deSolve::dede(y = y0_in_2, times = ts_pred_2, func = ODE_MODEL, parms = list(theta = theta, x_r = x_r, x_i = x_i),
                          atol = 1.0E-8,  rtol = 1.0E-5, maxsteps = 1000, control = list(mxhist = 1e5))
}else{
  # Run ODE solver for the second segment with the wrapper function
  y_out_2 = deSolve::dede(y = y0_in_2, times = ts_pred_2, func = ODE_MODEL_wrapper, parms = list(theta = theta, x_r = x_r, x_i = x_i),
                          atol = 1.0E-8,  rtol = 1.0E-5, maxsteps = 1000, control = list(mxhist = 1e5))
}


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

