## Install and load required packages if I don't have them
if (!require("rstan")) install.packages("rstan")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("HDInterval")) install.packages("HDInterval")
library(rstan)
library(ggplot2)
library(dplyr)
library(HDInterval)
options(mc.cores = parallel::detectCores())

#--------------------- process the data to stan format -------------------------

setwd("~/Desktop/Thesis/code/prime")
biom.df.filtered <- readRDS("biom.df.filtered.rds")

dat_long <- biom.df.filtered %>%
  arrange(SUBJID, BIOMYR) %>%
  mutate(subj = as.integer(factor(SUBJID)),
         t = BIOMYR,
         y = BIOMVAL) %>%
  group_by(subj) %>%
  mutate(base = first(BIOMVAL)) %>%  # first after arranging = earliest BIOMVAL
  ungroup()

# One base value per subject (J-length vector)
base_subject <- dat_long %>%
  distinct(subj, base) %>%
  arrange(subj) %>%
  pull(base)

stan_prime_dat <- list(
  N = nrow(dat_long),
  J = length(unique(dat_long$subj)),
  group = dat_long$subj,
  x = dat_long$t,
  y = dat_long$y,
  base = base_subject        # length J
)

str(stan_prime_dat$base)  # check

#------------------------ fit the model using rstan -----------------------------

set.seed(123456)

# Define the Stan model with random effects
tgi_prime_unif_model = "
  // where b_i's ~ N(0, sigma_b^2) and epsilon ~ N(0, sigma_eps^2)
  data {
    int<lower=0> N;              // number of observations
    int<lower=0> J;              // number of groups/subjects
    vector[N] x;                 // predictor (time) variable
    vector[N] y;                 // response variable
    vector[J] base;             // baseline predictor variable
    int<lower=1,upper=J> group[N]; // group/subject indicator for each obs
  }
  
  parameters {
    real<lower=0.01,upper=5> theta_g; // growth parameter
    real<lower=1,upper=10> theta_s; // decay parameter 
    vector[J] b_g_raw;             // non-centered parametstn for rand effs
    vector[J] b_s_raw;             // non-centered parametstn for rand effs
    
    real<lower=0> sigma_epsilon;    // measurement error standard deviation
    real<lower=0> sigma_b_g;             // random effects standard deviation
    real<lower=0> sigma_b_s;             // random effects standard deviation
  }
  
  transformed parameters {
    vector[J] b_g;                 // growth random effects for each subject
    vector[J] b_s;                 // shrink. random effects for each subject
    
    // Apply non-centered parameterization for computational efficiency
    b_g = sigma_b_g * b_g_raw;
    b_s = sigma_b_s * b_s_raw;
  }
  
  model {
    vector[N] mu;                // mean response for each observation
    
    // half-cauchy prior on σ_ε
    sigma_epsilon ~ cauchy(0, 0.1);
    
    // Half-normal prior on σ_b (can change to half_cauchy if preferred)
    sigma_b_g ~ normal(0, 1);      // half-normal due to <lower=0> constraint
    sigma_b_s ~ normal(0, 1);      // half-normal due to <lower=0> constraint
  
    // Normal prior for raw random effects (non-centered parameterization)
    b_g_raw ~ std_normal();
    b_s_raw ~ std_normal();
    
    for (i in 1:N) {
      mu[i] = base[group[i]]*(exp(theta_g * exp(b_g[group[i]]) * x[i]) + exp(-theta_s * exp(b_s[group[i]]) * x[i]) - 1);
    }
    
    // Likelihood
    y ~ normal(mu, sigma_epsilon);
  }
  
  generated quantities {
  // for posterior predictive checks
    vector[N] y_pred;            // posterior predictive samples
    vector[N] mu_pred;
    
    for (i in 1:N) {
      mu_pred[i] = base[group[i]]*(exp(theta_g * exp(b_g[group[i]]) * x[i]) + exp(-theta_s * exp(b_s[group[i]]) * x[i]) - 1);
      y_pred[i] = normal_rng(mu_pred[i], sigma_epsilon);
    }
  }
"

# Write the Stan model to a file
writeLines(tgi_prime_unif_model, "tgi_prime_unif_model.stan")
    
# get the posterior sample for theta
init_fun <- function() {
  list(
    theta_g = 0.5,  # within your bounds [0.135, 2.7]
    theta_s = 1.1   # within your bounds [0.135, 7.5]
  )
}

nchains = 1

init_list <- lapply(1:nchains, function(id) init_fun())

# Compile and fit the Stan model
fit_unif <- stan(
  file = "tgi_prime_unif_model.stan",
  data = stan_prime_dat,
  chains = nchains,               # number of Markov chains
  iter = 10000,              # total number of iterations per chain
  warmup = 5000,            # number of warmup iterations per chain
  thin = 1,                 # period for saving samples
  verbose = FALSE,
  init = init_list
)

rstan::traceplot(fit_unif, pars = c("theta_g", "theta_s"))

stan_unif_summary <- summary(fit_unif)$summary

unif_posterior_samples <- rstan::extract(fit_unif, permuted = TRUE)

theta_rhat <- stan_unif_summary[c("theta_g","theta_s"), "Rhat"]
theta_pointest <- stan_unif_summary[c("theta_g","theta_s"), "mean"]
print(theta_pointest)

biom.datf <- as.data.frame(biom.df.filtered)
biom.datf$pred_unif_stan <- colMeans(unif_posterior_samples$y_pred)

# Comparison plot: Frequentist dashed, Bayesian solid
# Palette & scale helpers ----------------------------------
approach_colors <- c(
  "Observed"    = "#000000",
  "Frequentist" = "#E69F00", #E69F00 # blue
  "Bayes"       = "#0072B2"   # orange
)

# approach_linetypes <- c(
#   "Observed"    = "blank",     # not used for points (ignored)
#   "Frequentist" = "dashed",
#   "Bayes"       = "solid"
# )

scale_color_approach <- function(...) {
  scale_color_manual(values = approach_colors, breaks = c("Observed", "Frequentist", "Bayes"), ...)
}

prime_id_subset <- unique(biom.datf$SUBJID)[1:12]

prime_subset <- biom.datf %>%
  filter(SUBJID %in% prime_id_subset) %>%
  mutate(SUBJID = factor(SUBJID, levels = prime_id_subset))  # preserves facet order

comparison_plot <- ggplot(prime_subset, aes(x = BIOMYR)) +
  geom_point(aes(y = BIOMVAL, color = "Observed"), size = 1) +
  geom_line(aes(y = BIOMVAL, color = "Observed"), linewidth = 0.4) +
  geom_point(aes(y = pred_unif_stan, color = "Bayes"), size = 1) +
  geom_line(aes(y = pred_unif_stan, color = "Bayes"), linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ SUBJID, scales = "free_y") +
  labs(
    title = "Bayesian TGI Model Fit with uniform prior on Prime Data",
    y = "SLD",
    x = "Time (years)",
    color = "Type"
  ) +
  scale_color_approach() +
  theme_minimal()

print(comparison_plot)
    

