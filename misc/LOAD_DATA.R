day_start           = 2 # most Bifidobacteriaceae transfer is assumed to be over at day 3, indexing from 0
fit_span            = 0.5
source('./misc/SMOOTH_DATA.R')
day_end_in          = 172+563 # max + 563
add_init_estimation = 21
rep_days            = new_172+add_init_estimation-3 # starts day 3
abundance_th        = 1e-4
source('./misc/LOAD_READ.R')
abundanceArray_meanSubjects_keep = abundanceArray_meanSubjects
days_array_keep                  = days_array