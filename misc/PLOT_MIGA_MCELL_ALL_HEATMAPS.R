keep_res_in_e  = keep_res_in[c('coeff_binding','mIgA_e','aff_e_final')]
keep_res_in_b  = keep_res_in[c('coeff_binding','mIgA_b','aff_b_final')]
keep_res_in_bc = keep_res_in[c('coeff_binding','mIgA_bc','aff_bc_final')]
keep_res_in_c  = keep_res_in[c('coeff_binding','mIgA_c','aff_c_final')]
keep_res_in    = keep_res_in %>% dplyr::rowwise() %>% dplyr::mutate(mIgA_ne = mean(mIgA_b,mIgA_bc,mIgA_c))
keep_res_in    = keep_res_in %>% dplyr::rowwise() %>% dplyr::mutate(aff_ne_final = mean(aff_b_final,aff_bc_final,aff_c_final))
keep_res_in_ne = keep_res_in[c('coeff_binding','mIgA_ne','aff_ne_final')]

taxon_name   = taxa_array[1]
keep_res_use = keep_res_in_e
source('./misc/PLOT_MIGA_MCELL_HEATMAP.R')
plot_e       = plot_out

taxon_name   = taxa_array[2]
keep_res_use = keep_res_in_b
source('./misc/PLOT_MIGA_MCELL_HEATMAP.R')
plot_b       = plot_out

taxon_name   = taxa_array[3]
keep_res_use = keep_res_in_bc
source('./misc/PLOT_MIGA_MCELL_HEATMAP.R')
plot_bc      = plot_out

taxon_name   = taxa_array[4]
keep_res_use = keep_res_in_c
source('./misc/PLOT_MIGA_MCELL_HEATMAP.R')
plot_c       = plot_out

taxon_name   = 'Symbiotic Commensals'
keep_res_use = keep_res_in_ne
source('./misc/PLOT_MIGA_MCELL_HEATMAP.R')
plot_ne      = plot_out