if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("HDInterval")) install.packages("HDInterval")
library(ggplot2)
library(dplyr)
library(HDInterval)

theta_g_fup_15 = read.csv("Sep-18-01h15-DataV2-mixed-effects-growth-posterior-sample-FUP-15_filtered.csv", header = TRUE)
theta_s_fup_15 = read.csv("Sep-18-01h15-DataV2-mixed-effects-shrink-posterior-sample-FUP-15_filtered.csv", header = TRUE)

theta_g_fup_50 = read.csv("Sep-18-12h14-DataV2-mixed-effects-growth-posterior-sample-FUP-50_filtered.csv", header = TRUE)
theta_s_fup_50 = read.csv("Sep-18-12h14-DataV2-mixed-effects-shrink-posterior-sample-FUP-50_filtered.csv", header = TRUE)

theta_g_fup_100 = read.csv("Sep-19-06h09-DataV2-mixed-effects-growth-posterior-sample-FUP-100_filtered.csv", header = TRUE)
theta_s_fup_100 = read.csv("Sep-19-06h09-DataV2-mixed-effects-shrink-posterior-sample-FUP-100_filtered.csv", header = TRUE)

# jeffreys
theta_g_jeff_15 = read.csv("Nov-18-01h28-DataV2-mixed-effects-growth-posterior-sample-jeffreys15_filtered.csv", header = TRUE)
theta_s_jeff_15 = read.csv("Nov-18-01h28-DataV2-mixed-effects-shrink-posterior-sample-jeffreys15_filtered.csv", header = TRUE)

theta_g_jeff_50 = read.csv("Nov-18-18h01-DataV2-mixed-effects-growth-posterior-sample-jeffreys50_filtered.csv", header = TRUE)
theta_s_jeff_50 = read.csv("Nov-18-18h01-DataV2-mixed-effects-shrink-posterior-sample-jeffreys50_filtered.csv", header = TRUE)

theta_g_jeff_100 = read.csv("Nov-19-14h17-DataV2-mixed-effects-growth-posterior-sample-jeffreys100_filtered.csv", header = TRUE)
theta_s_jeff_100 = read.csv("Nov-19-14h17-DataV2-mixed-effects-shrink-posterior-sample-jeffreys100_filtered.csv", header = TRUE)

# uniform NARROW
theta_g_unif_narrow_15 = read.csv("Nov-03-13h41-DataV2-mixed-effects-growth-posterior-sample-uniform-NARROW-15_filtered.csv", header = TRUE)
theta_s_unif_narrow_15 = read.csv("Nov-03-13h41-DataV2-mixed-effects-shrink-posterior-sample-uniform-NARROW-15_filtered.csv", header = TRUE)

theta_g_unif_narrow_50 = read.csv("Nov-07-00h52-DataV2-mixed-effects-growth-posterior-sample-uniform-NARROW-50_filtered.csv", header = TRUE)
theta_s_unif_narrow_50 = read.csv("Nov-07-00h52-DataV2-mixed-effects-shrink-posterior-sample-uniform-NARROW-50_filtered.csv", header = TRUE)

theta_g_unif_narrow_100 = read.csv("Nov-07-11h41-DataV2-mixed-effects-growth-posterior-sample-uniform-NARROW-100_filtered.csv", header = TRUE)
theta_s_unif_narrow_100 = read.csv("Nov-07-11h41-DataV2-mixed-effects-shrink-posterior-sample-uniform-NARROW-100_filtered.csv", header = TRUE)

# uniform MEDIUM
theta_g_unif_med_15 = read.csv("Nov-09-17h17-DataV2-mixed-effects-growth-posterior-sample-uniform-MEDIUM-15_filtered.csv", header = TRUE)
theta_s_unif_med_15 = read.csv("Nov-09-17h17-DataV2-mixed-effects-shrink-posterior-sample-uniform-MEDIUM-15_filtered.csv", header = TRUE)

theta_g_unif_med_50 = read.csv("Nov-10-01h29-DataV2-mixed-effects-growth-posterior-sample-uniform-MEDIUM-50_filtered.csv", header = TRUE)
theta_s_unif_med_50 = read.csv("Nov-10-01h29-DataV2-mixed-effects-shrink-posterior-sample-uniform-MEDIUM-50_filtered.csv", header = TRUE)

theta_g_unif_med_100 = read.csv("Nov-10-15h31-DataV2-mixed-effects-growth-posterior-sample-uniform-MEDIUM-100_filtered.csv", header = TRUE)
theta_s_unif_med_100 = read.csv("Nov-10-15h31-DataV2-mixed-effects-shrink-posterior-sample-uniform-MEDIUM-100_filtered.csv", header = TRUE)

# uniform WIDE
theta_g_unif_wide_15 = read.csv("Sep-15-22h37-DataV2-mixed-effects-growth-posterior-sample-uniform-WIDE-15_filtered.csv", header = TRUE)
theta_s_unif_wide_15 = read.csv("Sep-15-22h37-DataV2-mixed-effects-shrink-posterior-sample-uniform-WIDE-15_filtered.csv", header = TRUE)

theta_g_unif_wide_50 = read.csv("Nov-07-20h28-DataV2-mixed-effects-growth-posterior-sample-uniform-WIDE-50_filtered.csv", header = TRUE)
theta_s_unif_wide_50 = read.csv("Nov-07-20h28-DataV2-mixed-effects-shrink-posterior-sample-uniform-WIDE-50_filtered.csv", header = TRUE)

theta_g_unif_wide_100 = read.csv("Nov-08-17h17-DataV2-mixed-effects-growth-posterior-sample-uniform-WIDE-100_filtered.csv", header = TRUE)
theta_s_unif_wide_100 = read.csv("Nov-08-17h17-DataV2-mixed-effects-shrink-posterior-sample-uniform-WIDE-100_filtered.csv", header = TRUE)

#################----------------------------------------------------------------

statistics_calc <- function(theta_true, samples, is_list = TRUE) {
  if (is_list)
  {
    # extract theta posteriors
    thetas = as.matrix(samples[["theta_MH"]])
  } else {
    thetas = as.matrix(samples)
  }
  
  hdi_interval90 <- apply(thetas, 2, function(x) hdi(x, .90))
  CP90_in <- hdi_interval90[1,] <= theta_true & hdi_interval90[2,] >= theta_true
  CP90 <- sum(CP90_in)/ncol(thetas)
  hdi_interval95 <- apply(thetas, 2, function(x) hdi(x, .95))
  CP95_in <- hdi_interval95[1,] <= theta_true & hdi_interval95[2,] >= theta_true
  CP95 <- sum(CP95_in)/ncol(thetas)
  
  # plot a chain for sense-check
  plot(thetas[,1], type = 'l', col = 'blue', xlab = 'Iteration', 
       ylab = 'Theta', main = 'Trace plot of Theta', lwd = 0.5, ylim = c(0,2.5))
  lines(thetas[,2], col = "lightblue", lwd=0.5)
  lines(thetas[,3], col = "lightblue1", lwd=0.5)
  lines(thetas[,4], col = "lightblue2", lwd=0.5)
  lines(thetas[,5], col = "lightblue3", lwd=0.5)
  lines(thetas[,6], col = "aquamarine1", lwd=0.5)
  lines(thetas[,7], col = "aquamarine2", lwd=0.5)
  lines(thetas[,8], col = "aquamarine3", lwd=0.5)
  lines(thetas[,9], col = "aquamarine4", lwd=0.5)
  abline(h=hdi_interval90[1,1], col = 'green')
  abline(h=hdi_interval90[2,1], col = 'green')
  abline(h=hdi_interval95[1,1], col = 'red')
  abline(h=hdi_interval95[2,1], col = 'red')
  
  ILE90 <- mean(hdi_interval90[2,] - hdi_interval90[1,])
  ILE95 <- mean(hdi_interval95[2,] - hdi_interval95[1,])
  
  posterior_means <- apply(thetas, 2, mean)
  posterior_median <- apply(thetas, 2, median)
  posterior_biases <- mean((posterior_means-theta_true)) #abs value here? because that is just MAE...
  
  # Bias_squared <- mean((posterior_means-theta_true))
  
  MAE_mean <- mean(abs(posterior_means-theta_true))
  
  MAE_median <- mean(abs(posterior_median-theta_true))
  
  results <- c(CP90, CP95, ILE90, ILE95, MAE_mean, MAE_median)
  return(results)
}

# make the samples matrices

n15_theta_s_samples_unif_narrow= as.matrix(theta_s_unif_narrow_15)
n15_theta_g_samples_unif_narrow= as.matrix(theta_g_unif_narrow_15)
n50_theta_s_samples_unif_narrow= as.matrix(theta_s_unif_narrow_50)
n50_theta_g_samples_unif_narrow= as.matrix(theta_g_unif_narrow_50)
n100_theta_s_samples_unif_narrow= as.matrix(theta_s_unif_narrow_100)
n100_theta_g_samples_unif_narrow= as.matrix(theta_g_unif_narrow_100)

n15_theta_s_samples_unif_med= as.matrix(theta_s_unif_med_15)
n15_theta_g_samples_unif_med= as.matrix(theta_g_unif_med_15)
n50_theta_s_samples_unif_med= as.matrix(theta_s_unif_med_50)
n50_theta_g_samples_unif_med= as.matrix(theta_g_unif_med_50)
n100_theta_s_samples_unif_med= as.matrix(theta_s_unif_med_100)
n100_theta_g_samples_unif_med= as.matrix(theta_g_unif_med_100)

n15_theta_s_samples_unif_wide= as.matrix(theta_s_unif_wide_15)
n15_theta_g_samples_unif_wide= as.matrix(theta_g_unif_wide_15)
n50_theta_s_samples_unif_wide= as.matrix(theta_s_unif_wide_50)
n50_theta_g_samples_unif_wide= as.matrix(theta_g_unif_wide_50)
n100_theta_s_samples_unif_wide= as.matrix(theta_s_unif_wide_100)
n100_theta_g_samples_unif_wide= as.matrix(theta_g_unif_wide_100)

n15_theta_s_samples_jeff = as.matrix(theta_s_jeff_15)
n15_theta_g_samples_jeff = as.matrix(theta_g_jeff_15)
n50_theta_s_samples_jeff = as.matrix(theta_s_jeff_50)
n50_theta_g_samples_jeff = as.matrix(theta_g_jeff_50)
n100_theta_s_samples_jeff = as.matrix(theta_s_jeff_100)
n100_theta_g_samples_jeff = as.matrix(theta_g_jeff_100)

n15_theta_s_samples_fup = as.matrix(theta_s_fup_15)
n15_theta_g_samples_fup = as.matrix(theta_g_fup_15)
n50_theta_s_samples_fup = as.matrix(theta_s_fup_50)
n50_theta_g_samples_fup = as.matrix(theta_g_fup_50)
n100_theta_s_samples_fup = as.matrix(theta_s_fup_100)
n100_theta_g_samples_fup = as.matrix(theta_g_fup_100)

# -------------------------- display shrinkage results ----------------------------------


stats15_s_unif_narrow  = statistics_calc(exp(0.83),  n15_theta_s_samples_unif_narrow, FALSE)
stats50_s_unif_narrow  = statistics_calc(exp(0.83),  n50_theta_s_samples_unif_narrow, FALSE)
stats100_s_unif_narrow = statistics_calc(exp(0.83), n100_theta_s_samples_unif_narrow, FALSE)

stats15_s_unif_med  = statistics_calc(exp(0.83),  n15_theta_s_samples_unif_med, FALSE)
stats50_s_unif_med  = statistics_calc(exp(0.83),  n50_theta_s_samples_unif_med, FALSE)
stats100_s_unif_med = statistics_calc(exp(0.83), n100_theta_s_samples_unif_med, FALSE)

stats15_s_unif_wide  = statistics_calc(exp(0.83),  n15_theta_s_samples_unif_wide, FALSE)
stats50_s_unif_wide  = statistics_calc(exp(0.83),  n50_theta_s_samples_unif_wide, FALSE)
stats100_s_unif_wide = statistics_calc(exp(0.83), n100_theta_s_samples_unif_wide, FALSE)

stats15_s_jeff  = statistics_calc(exp(0.83),  n15_theta_s_samples_jeff, FALSE)
stats50_s_jeff  = statistics_calc(exp(0.83),  n50_theta_s_samples_jeff, FALSE)
stats100_s_jeff = statistics_calc(exp(0.83), n100_theta_s_samples_jeff, FALSE)

stats15_s_fup  = statistics_calc(exp(0.83),  n15_theta_s_samples_fup, FALSE)
stats50_s_fup  = statistics_calc(exp(0.83),  n50_theta_s_samples_fup, FALSE)
stats100_s_fup = statistics_calc(exp(0.83), n100_theta_s_samples_fup, FALSE)

tgi_shrink_statistics_table_unif_narrow = rbind(stats15_s_unif_narrow,
                                                stats50_s_unif_narrow,
                                                stats100_s_unif_narrow)

tgi_shrink_statistics_table_unif_med = rbind(stats15_s_unif_med,
                                             stats50_s_unif_med,
                                             stats100_s_unif_med)

tgi_shrink_statistics_table_unif_wide = rbind(stats15_s_unif_wide,
                                         stats50_s_unif_wide,
                                         stats100_s_unif_wide)

tgi_shrink_statistics_table_jeff = rbind(stats15_s_jeff,
                                         stats50_s_jeff, 
                                         stats100_s_jeff)

tgi_shrink_statistics_table_fup = rbind(stats15_s_fup,
                                         stats50_s_fup, 
                                         stats100_s_fup)

row.names(tgi_shrink_statistics_table_unif_narrow) = c("s unif narrow n=15", "s unif narrow n=50", "s unif narrow n=100")
colnames(tgi_shrink_statistics_table_unif_narrow) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                               "MAE_mean", "MAE_median")

row.names(tgi_shrink_statistics_table_unif_med) = c("s unif med n=15", "s unif med n=50", "s unif med n=100")
colnames(tgi_shrink_statistics_table_unif_med) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                                      "MAE_mean", "MAE_median")

row.names(tgi_shrink_statistics_table_unif_wide) = c("s unif wide n=15", "s unif wide n=50", "s unif wide n=100")
colnames(tgi_shrink_statistics_table_unif_wide) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                                      "MAE_mean", "MAE_median")

row.names(tgi_shrink_statistics_table_jeff) = c("s jeff n=15", "s jeff n=50", "s jeff n=100")
colnames(tgi_shrink_statistics_table_jeff) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                               "MAE_mean", "MAE_median")


row.names(tgi_shrink_statistics_table_fup) = c("s fup n=15", "s fup n=50", "s fup n=100")
colnames(tgi_shrink_statistics_table_fup) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                              "MAE_mean", "MAE_median")


shrink_stats = rbind(round(tgi_shrink_statistics_table_unif_narrow,5),
                     round(tgi_shrink_statistics_table_unif_med,5),
                     round(tgi_shrink_statistics_table_unif_wide,5),
                     round(tgi_shrink_statistics_table_jeff,5),
                     round(tgi_shrink_statistics_table_fup,5))

# -------------------------- display growth results ----------------------------------


stats15_g_unif_narrow  = statistics_calc(exp(-2),  n15_theta_g_samples_unif_narrow, FALSE)
stats50_g_unif_narrow  = statistics_calc(exp(-2),  n50_theta_g_samples_unif_narrow, FALSE)
stats100_g_unif_narrow = statistics_calc(exp(-2), n100_theta_g_samples_unif_narrow, FALSE)

stats15_g_unif_med  = statistics_calc(exp(-2),  n15_theta_g_samples_unif_med, FALSE)
stats50_g_unif_med  = statistics_calc(exp(-2),  n50_theta_g_samples_unif_med, FALSE)
stats100_g_unif_med = statistics_calc(exp(-2), n100_theta_g_samples_unif_med, FALSE)

stats15_g_unif_wide  = statistics_calc(exp(-2),  n15_theta_g_samples_unif_wide, FALSE)
stats50_g_unif_wide  = statistics_calc(exp(-2),  n50_theta_g_samples_unif_wide, FALSE)
stats100_g_unif_wide = statistics_calc(exp(-2), n100_theta_g_samples_unif_wide, FALSE)

stats15_g_jeff  = statistics_calc(exp(-2),  n15_theta_g_samples_jeff, FALSE)
stats50_g_jeff  = statistics_calc(exp(-2),  n50_theta_g_samples_jeff, FALSE)
stats100_g_jeff = statistics_calc(exp(-2), n100_theta_g_samples_jeff, FALSE)

stats15_g_fup  = statistics_calc(exp(-2),  n15_theta_g_samples_fup, FALSE)
stats50_g_fup  = statistics_calc(exp(-2),  n50_theta_g_samples_fup, FALSE)
stats100_g_fup = statistics_calc(exp(-2), n100_theta_g_samples_fup, FALSE)

tgi_growth_statistics_table_unif_narrow = rbind(stats15_g_unif_narrow,
                                                stats50_g_unif_narrow,
                                                stats100_g_unif_narrow)

tgi_growth_statistics_table_unif_med = rbind(stats15_g_unif_med,
                                             stats50_g_unif_med,
                                             stats100_g_unif_med)

tgi_growth_statistics_table_unif_wide = rbind(stats15_g_unif_wide,
                                              stats50_g_unif_wide,
                                              stats100_g_unif_wide)

tgi_growth_statistics_table_jeff = rbind(stats15_g_jeff,
                                         stats50_g_jeff, 
                                         stats100_g_jeff)

tgi_growth_statistics_table_fup = rbind(stats15_g_fup,
                                        stats50_g_fup, 
                                        stats100_g_fup)

row.names(tgi_growth_statistics_table_unif_narrow) = c("g unif narrow n=15", "g unif narrow n=50", "g unif narrow n=100")
colnames(tgi_growth_statistics_table_unif_narrow) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                                      "MAE_mean", "MAE_median")

row.names(tgi_growth_statistics_table_unif_med) = c("g unif med n=15", "g unif med n=50", "g unif med n=100")
colnames(tgi_growth_statistics_table_unif_med) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                                   "MAE_mean", "MAE_median")

row.names(tgi_growth_statistics_table_unif_wide) = c("g unif wide n=15", "g unif wide n=50", "g unif wide n=100")
colnames(tgi_growth_statistics_table_unif_wide) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                                    "MAE_mean", "MAE_median")

row.names(tgi_growth_statistics_table_jeff) = c("g jeff n=15", "g jeff n=50", "g jeff n=100")
colnames(tgi_growth_statistics_table_jeff) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                               "MAE_mean", "MAE_median")


row.names(tgi_growth_statistics_table_fup) = c("g fup n=15", "g fup n=50", "g fup n=100")
colnames(tgi_growth_statistics_table_fup) = c("CP_90", "CP_95", "ILE90", "ILE95",
                                              "MAE_mean", "MAE_median")


growth_stats = rbind(round(tgi_growth_statistics_table_unif_narrow,5),
                     round(tgi_growth_statistics_table_unif_med,5),
                     round(tgi_growth_statistics_table_unif_wide,5),
                     round(tgi_growth_statistics_table_jeff,5),
                     round(tgi_growth_statistics_table_fup,5))
#----------------------- export to csv -----------------------------------------


filename = paste(format(Sys.time(), "%b-%d-%Hh%M"),"-mixed-effects-posterior-statistics.csv",sep="")
write.csv(as.data.frame(rbind(shrink_stats, growth_stats)), filename, row.names = TRUE)

