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

p_cals=ggplot(df_long, aes(x = t, y = value, color = series, linetype = series)) +
  geom_line(data = df_long, size=0.7) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  labs(title ='Pre-computed model inputs', x = "DOL", y = "") +
  theme(legend.title = element_blank(), plot.title = element_text(size = 12)) +
  guides(size = "none") + theme(legend.title = element_blank(),
                                # legend.text = element_text(size = 12),
                                # axis.title = element_text(size = 12),
                                # axis.text = element_text(size = 12),
                                legend.position = c(0.95, 0.074), # This needs adjustment
                                legend.justification = c("right", "bottom"), # Anchor point
                                legend.box.just = "right",
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
# Print the plot
print(p_cals)


# graphics.off()
# png(file =paste0("/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan/PAPERPLOTS/PAPERPLOTS_CALORIES_",EBF_duration,".png"),   # The directory you want to save the file in
#     width     = 7,
#     height    = 3,
#     units     = "in",
#     res       = 600)
# print(p)
# dev.off()
