## Install and load required packages if I don't have them
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("HDInterval")) install.packages("HDInterval")
library(ggplot2)
library(dplyr)
library(HDInterval)

# ---
# SCRIPT CONFIGURATION
# ---
# Set working directory and load data
setwd("~/Desktop/Thesis/code/prime")

# Load the data frame containing actual values and model predictions
# This file should have BIOMVAL, pred_fup_stan, pred_unif_stan, pred_jeff_stan
biom.datf <- readRDS("biom.df.fitted.rds")

# ---
# ASSUMPTION: Load or create Stan summary objects
# ---

# ---
# 1. Predictive Performance Metric Functions
# ---
mae <- function(actual, predicted) {
  mean(abs(actual - predicted), na.rm = TRUE)
}

mse <- function(actual, predicted) {
  mean((actual - predicted)^2, na.rm = TRUE)
}

rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

r_squared <- function(actual, predicted) {
  # Calculate residuals
  ss_res <- sum((actual - predicted)^2, na.rm = TRUE)
  # Calculate total sum of squares
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  # Calculate R-squared
  1 - (ss_res / ss_tot)
}

# ---
# 2. Calculate Predictive Performance Metrics
# ---
# Get the actual observed values
actual_values <- biom.datf$BIOMVAL

# List of prediction column names to evaluate
predictor_cols <- c("pred_fup_stan", "pred_unif_stan", "pred_jeff_stan")
model_names <- c("FUP Prior", "Uniform Prior", "Jeffreys Prior") # Friendly names for output

# Create an empty data.frame to store results
performance_results <- data.frame(
  Model = character(),
  MAE = numeric(),
  MSE = numeric(),
  RMSE = numeric(),
  R_Squared = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each predictor and calculate metrics
for (i in seq_along(predictor_cols)) {
  pred_col_name <- predictor_cols[i]
  model_name <- model_names[i]
  
  # Get the predicted values for the current model
  predicted_values <- biom.datf[[pred_col_name]]
  
  # Calculate metrics
  current_mae <- mae(actual_values, predicted_values)
  current_mse <- mse(actual_values, predicted_values)
  current_rmse <- rmse(actual_values, predicted_values)
  current_r2 <- r_squared(actual_values, predicted_values)
  
  # Add to results data.frame
  performance_results[i, ] <- list(
    model_name, 
    current_mae, 
    current_mse, 
    current_rmse, 
    current_r2
  )
}

# ---
# 3. Print Predictive Performance Results
# ---
cat("--- Predictive Performance Comparison ---\n\n")
print(performance_results, digits = 4)

cat("\n--- Interpretation (Performance) ---\n")
cat("Lower MAE, MSE, and RMSE indicate better fit (less error).\n")
cat("Higher R-Squared (closer to 1) indicates a better fit.\n")


# ---
# 4. Extract Parameter Estimates and Convergence
# ---
# This section assumes the stan summary objects are loaded
# (e.g., stan_fup_summary, stan_unif_summary, stan_jeff_summary)

cat("\n\n--- Parameter Estimates and Convergence ---\n\n")

# Check if all required summary objects exist
if (exists("stan_fup_summary") && 
    exists("stan_unif_summary") && 
    exists("stan_jeff_summary")) {
  
  # List of summary objects
  all_summaries <- list("Uniform Prior" = stan_unif_summary,
                        "Jeffreys Prior" = stan_jeff_summary,
                        "FUP Prior" = stan_fup_summary)
  
  # Parameters and columns to extract
  params_to_extract <- c("theta_g", "theta_s") # Add other params if needed
  cols_to_extract <- c("mean", "se_mean","2.5%", "97.5%", "Rhat")
  
  # Use lapply to extract and format data from each summary
  results_list <- lapply(names(all_summaries), function(model_name) {
    summary_df <- all_summaries[[model_name]]
    
    # Check if parameters exist in summary
    existing_params <- params_to_extract[params_to_extract %in% rownames(summary_df)]
    
    if (length(existing_params) == 0) {
      cat(paste("Warning: None of", 
                paste(params_to_extract, collapse=", "), 
                "found in", model_name, "summary.\n"))
      return(NULL)
    }
    
    # Extract the relevant rows and columns
    # Use 'drop = FALSE' to ensure it stays a data.frame even if 1 param
    param_data_matrix <- summary_df[existing_params, cols_to_extract, drop = FALSE]
    param_data <- as.data.frame(param_data_matrix) # Convert matrix to data.frame
    
    # Add model and parameter columns
    param_data$Model <- model_name
    param_data$Parameter <- rownames(param_data)
    
    # Rename columns for clarity
    colnames(param_data)[colnames(param_data) == "2.5%"] <- "CI_Lower"
    colnames(param_data)[colnames(param_data) == "97.5%"] <- "CI_Upper"
    
    return(param_data)
  })
  
  # Combine all results into one data.frame
  parameter_results <- do.call(rbind, results_list)
  
  # Reorder columns for readability
  if (!is.null(parameter_results)) {
    parameter_results <- parameter_results[, c("Model", "Parameter", "mean", 
                                               "CI_Lower", "CI_Upper", "Rhat")]
    rownames(parameter_results) <- NULL
    
    # Print the combined parameter results
    CI_len <- parameter_results$CI_Upper-parameter_results$CI_Lower
    parameter_results <- cbind(parameter_results, CI_len)
    print(parameter_results, digits = 5)
    
    cat("\n--- Interpretation (Parameters) ---\n")
    cat("mean: Posterior mean estimate\n")
    cat("CI_Lower, CI_Upper: 95% Equal-Tailed Credible Interval\n")
    cat("Rhat: Convergence diagnostic (should be <= 1.01).\n")
  }
  
} else {
  cat("Warning: Could not find one or more Stan summary objects.\n")
  cat("Please ensure 'stan_fup_summary', 'stan_unif_summary', and \n")
  cat("'stan_jeff_summary' objects exist in your R environment.\n")
  cat("Skipping parameter estimate comparison.\n")
}

# --- End of Script ---

stan_unif_summary[c("theta_g","theta_s"),]
stan_jeff_summary[c("theta_g","theta_s"),]
stan_fup_summary[c("theta_g","theta_s"),]

unif_g = unif_posterior_samples$theta_g
unif_s = unif_posterior_samples$theta_s

jeff_g = jeff_posterior_samples$theta_g
jeff_s = jeff_posterior_samples$theta_s

fup_g = fup_posterior_samples$theta_g
fup_s = fup_posterior_samples$theta_s

nadir <- function(thetag, theta_s)
{
  return(log(theta_s/thetag)/(theta_s+thetag))
}

nadir_unif = nadir(unif_g, unif_s)
nadir_jeff = nadir(jeff_g, jeff_s)
nadir_fup = nadir(fup_g, fup_s)

nadir_stats <- function(nadir_samps) {
  hdi_interval90 <- hdi(nadir_samps, .90)
  hdi_interval95 <- hdi(nadir_samps, .95)
  
  posterior_means <- mean(nadir_samps)
  
  results <- data.frame(posterior_mean = posterior_means,HDI_95_lower = hdi_interval95[1], HDI_95_upper = hdi_interval95[2])
  return(results)
}

nadir_stats(nadir_unif)
nadir_stats(nadir_jeff)
nadir_stats(nadir_fup)

#----------------------make plots ----------------------------------------------

df_nadir_unif <- data.frame(
  value = nadir_unif,
  sample_group = "Uniform"
)

df_nadir_jeff <- data.frame(
  value = nadir_jeff,
  sample_group = "Jeffreys"
)

df_nadir_fup <- data.frame(
  value = nadir_fup,
  sample_group = "Functional"
)

# Combine them into one data frame
all_samples_df <- rbind(df_nadir_unif, df_nadir_jeff, df_nadir_fup)

# --- 3. Create the Overlapping Density Plot ---
my_custom_colors <- c("#E69F00", "#6e8a4e", "#0072B2") 

density_plot <- ggplot(all_samples_df, aes(x = value)) +
  geom_density(aes(fill = sample_group, color = sample_group), alpha = 0.3) +
  scale_fill_manual(values = my_custom_colors) +
  scale_color_manual(values = my_custom_colors) +
  
  # (Optional) Add labels and a title
  labs(
    title = "Posterior Densities of Tumour Nadir",
    x = "Time to nadir",
    y = "Posterior density",
    fill = "Prior",
    color = "Prior"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(density_plot)