## Install and load required packages if I don't have them
if (!require("rstan")) install.packages("rstan")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("HDInterval")) install.packages("HDInterval")
if (!require("nlme")) install.packages("nlme")
if (!require("datasets")) install.packages("datasets")

library(datasets)
library(nlme)
library(rstan)
library(ggplot2)
library(dplyr)
library(HDInterval)
options(mc.cores = parallel::detectCores())

# Source the updated colors file
source("colors.R")

############################ FREQUENTIST #######################################

freqmod <- nlme(conc ~ SSfol(Dose, Time, lKe, lKa, lCl),
                data = Theoph,
                fixed = list(lKe + lKa + lCl ~ 1),
                random = lKe + lKa + lCl ~ 1 | Subject,
                start = c(lKe = -2.5, lKa = 0.5, lCl = -3),
                control = nlmeControl(maxIter = 2000, tolerance = 1e-6))
summary(freqmod)

# Fixed effects
fixed.effects(freqmod)

# Random effects
random_effects = ranef(freqmod)
random_effects

# Extract variance components
variance_components <- VarCorr(freqmod)

# Print the variance components
print(variance_components)

Theoph$Subject <- as.numeric(as.factor(Theoph$Subject))

############################## BAYESIAN ########################################

# Set up Stan options for better performance
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load the Theoph dataset
data(Theoph)
Theoph$Subject <- as.factor(Theoph$Subject)

# Prepare data for Stan
theoph_stan_data <- list(
  N = nrow(Theoph),
  J = length(unique(Theoph$Subject)),
  Time = Theoph$Time,
  conc = Theoph$conc,
  Dose = Theoph$Dose,
  Subject = as.numeric(as.factor(Theoph$Subject))
)

# ============================================================================
# Stan Model Definition - One-Compartment with Mixed Effects
# ============================================================================

theoph_stan_model <- "
data {
  int<lower=0> N;              // number of observations
  int<lower=0> J;              // number of subjects
  vector[N] Time;              // time points
  vector[N] Dose;              // dose amounts
  vector[N] conc;              // observed concentrations
  int<lower=1,upper=J> Subject[N]; // subject indicator for each observation
}

parameters {
  // Population parameters (log-scale)
  real beta_lKe;               // population log elimination rate
  real beta_lKa;               // population log absorption rate  
  real beta_lCl;               // population log clearance
  
  // Random effects (non-centered parameterization)
  vector[J] eta_lKe;           // individual deviations for lKe
  vector[J] eta_lKa;           // individual deviations for lKa
  vector[J] eta_lCl;           // individual deviations for lCl
  
  // Variance parameters
  real<lower=0> omega_lKe;     // between-subject SD for lKe
  real<lower=0> omega_lKa;     // between-subject SD for lKa
  real<lower=0> omega_lCl;     // between-subject SD for lCl
  
  real<lower=0> sigma;         // residual error SD
}

transformed parameters {
  vector[J] lKe_i;             // individual log elimination rates
  vector[J] lKa_i;             // individual log absorption rates
  vector[J] lCl_i;             // individual log clearances
  
  vector[J] Ke_i;              // individual elimination rates
  vector[J] Ka_i;              // individual absorption rates
  vector[J] Cl_i;              // individual clearances
  vector[J] V_i;               // individual volumes
  
  vector[N] conc_pred;         // predicted concentrations
  
  // Individual parameters (non-centered parameterization)
  lKe_i = beta_lKe + omega_lKe * eta_lKe;
  lKa_i = beta_lKa + omega_lKa * eta_lKa;
  lCl_i = beta_lCl + omega_lCl * eta_lCl;
  
  // Transform to natural scale
  Ke_i = exp(lKe_i);
  Ka_i = exp(lKa_i);
  Cl_i = exp(lCl_i);
  V_i = Cl_i ./ Ke_i;          // V = Cl/Ke
  
  // Calculate predicted concentrations
  for (n in 1:N) {
    int subj = Subject[n];
    real t = Time[n];
    real dose = Dose[n];
    
    // One-compartment model with first-order absorption
    // C(t) = (Dose * Ka / V / (Ka - Ke)) * (exp(-Ke*t) - exp(-Ka*t))
    if (Ka_i[subj] != Ke_i[subj]) {
      conc_pred[n] = (dose * Ka_i[subj] / V_i[subj] / (Ka_i[subj] - Ke_i[subj])) * 
                     (exp(-Ke_i[subj] * t) - exp(-Ka_i[subj] * t));
    } else {
      // Handle edge case where Ka = Ke (should be rare)
      conc_pred[n] = dose * Ka_i[subj] * t / V_i[subj] * exp(-Ke_i[subj] * t);
    }
  }
}

model {
  // Priors for population parameters (based on frequentist results)
  beta_lKe ~ normal(-2.5, 1);  // log(0.08) ≈ -2.5
  beta_lKa ~ normal(0.5, 1);   // log(1.5) ≈ 0.4
  beta_lCl ~ normal(-1, 1);   // log(1.5) ≈ 0.4
  
  // Priors for between-subject variability
  omega_lKe ~ normal(0, 0.5);
  omega_lKa ~ normal(0, 0.5);
  omega_lCl ~ normal(0, 0.5);
  
  // Prior for residual error
  sigma ~ normal(0, 2);
  
  // Random effects (standard normal for non-centered parameterization)
  eta_lKe ~ std_normal();
  eta_lKa ~ std_normal();
  eta_lCl ~ std_normal();
  
  // Likelihood
  conc ~ normal(conc_pred, sigma);
}

generated quantities {
  // Population-level derived parameters
  real Ke_pop = exp(beta_lKe);
  real Ka_pop = exp(beta_lKa);
  real Cl_pop = exp(beta_lCl);
  real V_pop = Cl_pop / Ke_pop;
  real half_life_pop = log(2) / Ke_pop;
  real tmax_pop = log(Ka_pop/Ke_pop) / (Ka_pop - Ke_pop);
  
  // Log-likelihood for model comparison
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(conc[n] | conc_pred[n], sigma);
  }
  
  // Posterior predictive checks
  vector[N] conc_rep;
  for (n in 1:N) {
    conc_rep[n] = normal_rng(conc_pred[n], sigma);
  }
}
"

# ============================================================================
# Initial Values Function
# ============================================================================

init_fun <- function() {
  list(
    beta_lKe = rnorm(1, -2.5, 0.1),
    beta_lKa = rnorm(1, 0.4, 0.1),
    beta_lCl = rnorm(1, -3, 0.1),
    eta_lKe = rnorm(theoph_stan_data$J, 0, 0.1),
    eta_lKa = rnorm(theoph_stan_data$J, 0, 0.1),
    eta_lCl = rnorm(theoph_stan_data$J, 0, 0.1),
    omega_lKe = abs(rnorm(1, 0, 0.1)),
    omega_lKa = abs(rnorm(1, 0, 0.1)),
    omega_lCl = abs(rnorm(1, 0, 0.1)),
    sigma = abs(rnorm(1, 0.5, 0.1))
  )
}

# ============================================================================
# Fit the Bayesian Model
# ============================================================================

cat("Fitting Bayesian one-compartment model...\n")
cat("This may take several minutes...\n")

# Compile and fit the Stan model
stan_fit <- stan(
  model_code = theoph_stan_model,
  data = theoph_stan_data,
  init = init_fun,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 2,
  cores = 4,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  verbose = FALSE
)

posterior_samples <- rstan::extract(stan_fit, permuted = TRUE)

############################## COMPARE #########################################

# FREQ predictions
Theoph$pred_nlme <- fitted(freqmod)

# BAYES predictions
# Compute posterior mean predictions
pred_stan <- colMeans(posterior_samples$conc_pred)
Theoph$pred_stan <- pred_stan

# Plot 1: Comparison of models on Theoph data
# Using Stellenbosch brand colors
ggplot(Theoph, aes(x = Time)) +
  geom_point(aes(y = conc, color = "Observed"), size = 2, alpha = 0.7) +
  geom_line(aes(y = pred_nlme, color = "Frequentist"), linewidth = 0.8) +
  geom_line(aes(y = pred_stan, color = "Bayesian"), linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~Subject, scales = "free_y") +
  scale_color_manual(
    values = c(
      "Observed" = "#4D4D4D",        # Neutral gray for observed data
      "Frequentist" = "#61223b",      # Confident Maroon (60% tint)
      "Bayesian" = "#caa258"          # Brilliant Gold
    )
  ) +
  labs(
    title = "Comparison of Models on Theoph Data",
    y = "Concentration (mg/L)", 
    x = "Time (hr)",
    color = "Type"
  ) +
  theme_thesis()

# Plot 2: Comparison of Fixed Effects
# Compare parameter estimates
nlme_params <- fixef(freqmod)
stan_params <- summary(stan_fit, pars = c("beta_lKe", "beta_lKa", "beta_lCl"))$summary[, "mean"]

params_df <- data.frame(
  Parameter = c("lKe", "lKa", "lCl"),
  Frequentist = nlme_params,
  Bayesian = stan_params
)

# Reshape for plotting
params_long <- params_df %>%
  tidyr::pivot_longer(cols = c(Frequentist, Bayesian), 
                      names_to = "Model", 
                      values_to = "Estimate")

ggplot(params_long, aes(x = Parameter, y = Estimate, color = Model, group = Model)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linewidth = 0.8) +
  scale_color_manual(
    values = c(
      "Frequentist" = "#61223b",  # Confident Maroon (60% tint)
      "Bayesian" = "#caa258"      # Brilliant Gold
    )
  ) +
  labs(
    title = "Comparison of Fixed Effects",
    y = "Estimate (log scale)", 
    x = "Parameter",
    color = "Model"
  ) +
  theme_thesis()

# Plot 3: Residual diagnostics
Theoph$resid_nlme <- residuals(freqmod)
Theoph$resid_stan <- Theoph$conc - Theoph$pred_stan

ggplot(Theoph) +
  geom_point(aes(y = resid_nlme, x = pred_nlme, color = "Frequentist"), 
             size = 2.5, alpha = 0.6) +
  geom_point(aes(y = resid_stan, x = pred_stan, color = "Bayesian"), 
             size = 2.5, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#1A1A1A", linewidth = 0.6) +
  scale_color_manual(
    values = c(
      "Frequentist" = "#97223b",  # Confident Maroon (60% tint)
      "Bayesian" = "#DEB845"      # Brilliant Gold
    )
  ) +
  labs(
    title = "Model Residuals vs. Fitted Values",
    x = "Fitted Values (mg/L)", 
    y = "Residuals (mg/L)",
    color = "Model"
  ) +
  theme_thesis()
