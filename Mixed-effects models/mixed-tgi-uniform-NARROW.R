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

set.seed(123456)

# Define the Stan model with random effects
tgi_unif_model_random_effects = "
  // where b_i's ~ N(0, sigma_b^2) and epsilon ~ N(0, sigma_eps^2)
  data {
    int<lower=0> N;              // number of observations
    int<lower=0> J;              // number of groups/subjects
    vector[N] x;                 // predictor variable
    vector[N] y;                 // response variable
    int<lower=1,upper=J> group[N]; // group/subject indicator for each obs
  }
  
  parameters {
    real<lower=0.01,upper=1> theta_g; // growth parameter
    real<lower=1,upper=4> theta_s; // decay parameter 
    vector[J] b_g_raw;             // non-centered parametstn for rand effs
    vector[J] b_s_raw;             // non-centered parametstn for rand effs
    
    real<lower=0> sigma2_epsilon;    // measurement error VARIANCE
    real<lower=0> sigma_b_g;             // random effects standard deviation
    real<lower=0> sigma_b_s;             // random effects standard deviation
  }
  
  transformed parameters {
    vector[J] b_g;                 // growth random effects for each subject
    vector[J] b_s;                 // shrink. random effects for each subject
    real sigma_epsilon;                // measurement error standard deviation
    
    // Apply non-centered parameterization for computational efficiency
    b_g = sigma_b_g * b_g_raw;
    b_s = sigma_b_s * b_s_raw;
    
    // Convert variance to standard deviation
    sigma_epsilon = sqrt(sigma2_epsilon);
    
  }
  
  model {
    vector[N] mu;                // mean response for each observation
    
    // Inverse gamma prior on σ_ε² 
    sigma2_epsilon ~ inv_gamma(0.001, 0.001);
    
    // Half-normal prior on σ_b (can change to half_cauchy if preferred)
    sigma_b_g ~ normal(0, 1);      // half-normal due to <lower=0> constraint
    sigma_b_s ~ normal(0, 1);      // half-normal due to <lower=0> constraint
  
    // Normal prior for raw random effects (non-centered parameterization)
    b_g_raw ~ std_normal();
    b_s_raw ~ std_normal();
    
    for (i in 1:N) {
      mu[i] = exp(theta_g * exp(b_g[group[i]]) * x[i]) + exp(-theta_s * exp(b_s[group[i]]) * x[i]) - 1;
    }
    
    // Uniform prior for thetas (implicit by parameter bounds, no need for explicit statement)
  
    // Likelihood
    y ~ normal(mu, sigma_epsilon);
  }
  
  generated quantities {
  // for posterior predictive checks
    vector[N] y_pred;            // posterior predictive samples
    vector[N] mu_pred;
    
    for (i in 1:N) {
      mu_pred[i] = exp(theta_g * exp(b_g[group[i]]) * x[i]) + exp(-theta_s * exp(b_s[group[i]]) * x[i]) - 1;
      y_pred[i] = normal_rng(mu_pred[i], sigma_epsilon);
    }
  }
"

# Write the Stan model to a file
writeLines(tgi_unif_model_random_effects, "tgi_mixed_model_unif.stan")

generate_tgi_dataset <- function(NumReps = 8, J, 
                                 true_theta_g = exp(-2), sigma_b_g=0.1,
                                 true_theta_s =  exp(0.83), sigma_b_s=0.7, sigma_eps=0.05) {
  # Generate some synthetic data for testing
  
  N <- NumReps*J              # total number of observations
  
  
  # Generate random effects for each group
  b_g <- rnorm(J, 0, sigma_b_g)
  b_s <- rnorm(J, 0, sigma_b_s)
  
  # Generate balanced design with equal observations per group
  group <- rep(1:J, each = N/J)
  
  # Create sequence of x values (could be different for each group too)
  x <- rep(seq(0, 2, length.out = NumReps), J)
  
  # Generate y values according to the model with random effects and noise
  y <- numeric(N)
  for (i in 1:N) {
    y0 = 1
    yg = (exp(true_theta_g * exp(b_g[group[i]]) * x[i])-1)
    ys = (exp(-true_theta_s * exp(b_s[group[i]]) * x[i]))
    y[i] <- y0*(yg+ys) + rnorm(1, 0, sigma_eps)
  }
  
  # Create data frame for ggplot
  sim_data <- data.frame(
    x = x,
    y = y,
    group = factor(group)
  )
  
  data_list <- list(
    N = N,
    J = J,
    x = x,
    y = y,
    group = group
  )
  
  return(data_list)
}

simulate_tgi_posterior_mixed <- function(nreps, nsubjects, 
                                         true_theta_g , sigma_b_g,
                                         true_theta_s, sigma_b_s, 
                                         sigma_eps, M, nchains) {
  
  theta_g_storage <- matrix(nrow = M*nchains, ncol = nreps)
  theta_s_storage <- matrix(nrow = M*nchains, ncol = nreps)
  
  rhat_storage <- matrix(nrow = nreps, ncol = 2)
  
  theta_pointest_storage <- matrix(nrow = nreps, ncol = 2)
  
  data_storage <- matrix(nrow = nsubjects*8, ncol = nreps)
  
  for (repl in 1:nreps) {
    print(paste("---------------Starting STAN for TGI sample", repl,"--------------------------"))
    # generate a dataset
    stan_data_list <- generate_tgi_dataset(NumReps = 8, J = nsubjects, 
                                           true_theta_g = true_theta_g, sigma_b_g=sigma_b_g,
                                           true_theta_s =  true_theta_s, sigma_b_s=sigma_b_s, 
                                           sigma_eps=sigma_eps)
    data_storage[,repl] <- stan_data_list$y
    
    # get the posterior sample for theta
    init_fun <- function() {
      list(
        theta_g = 0.5,  # within your bounds [0.135, 2.7]
        theta_s = 1.5   # within your bounds [0.135, 7.5]
      )
    }
    init_list <- lapply(1:nchains, function(id) init_fun())
    
    # Compile and fit the Stan model
    fit <- stan(
      file = "tgi_mixed_model_unif.stan",
      data = stan_data_list,
      chains = nchains,               # number of Markov chains
      iter = M+3000,              # total number of iterations per chain
      warmup = 3000,            # number of warmup iterations per chain
      thin = 1,                 # period for saving samples
      verbose = FALSE,
      init = init_list,
      control = list(
        adapt_delta = 0.95        # increased from 0.95
      )
    )
    
    if (repl<10) {
      rstan::traceplot(fit, pars = c("theta_g", "theta_s"))
    }
    # store the sample and Rhat values
    stan_summary <- summary(fit)$summary
    posterior_samples <- as.array(fit, pars = c("theta_g", "theta_s"))
    
    theta_rhat <- stan_summary[c("theta_g","theta_s"), "Rhat"]
    theta_pointest <- stan_summary[c("theta_g","theta_s"), "mean"]
    print(theta_pointest)
    
    theta_g_storage[,repl] <- as.vector(posterior_samples[,,"theta_g"])
    theta_s_storage[,repl] <- as.vector(posterior_samples[,,"theta_s"])
    
    rhat_storage[repl,] <- t(as.matrix(theta_rhat))
    theta_pointest_storage[repl,] <- t(as.matrix(theta_pointest))
  }
  
  colnames(theta_g_storage) <- paste(rep("sample", nreps), 1:nreps, sep='')
  filename = paste(format(Sys.time(), "%b-%d-%Hh%M"),"-DataV2-mixed-effects-growth-posterior-sample-uniform-NARROW-",nsubjects,".csv",sep="")
  write.csv(as.data.frame(theta_g_storage), filename, row.names = FALSE)
  
  colnames(theta_s_storage) <- paste(rep("sample", nreps), 1:nreps, sep='')
  filename = paste(format(Sys.time(), "%b-%d-%Hh%M"),"-DataV2-mixed-effects-shrink-posterior-sample-uniform-NARROW-",nsubjects,".csv",sep="")
  write.csv(as.data.frame(theta_s_storage), filename, row.names = FALSE)
  
  return(list(theta_g = theta_g_storage, theta_s = theta_s_storage, rhats = rhat_storage, pointests = theta_pointest_storage, datasets = data_storage))
}

thetas_n15 <- simulate_tgi_posterior_mixed(nreps = 250, nsubjects=15,
                                           true_theta_g = exp(-2), sigma_b_g=0.1,
                                           true_theta_s =  exp(0.83), sigma_b_s=0.7,
                                           sigma_eps=0.05, M = 7000, nchains =1)

thetas_n50 <- simulate_tgi_posterior_mixed(nreps = 250, nsubjects=50,
                                           true_theta_g = exp(-2), sigma_b_g=0.1,
                                           true_theta_s =  exp(0.83), sigma_b_s=0.7,
                                           sigma_eps=0.05, M = 7000, nchains =1)

thetas_n100 <- simulate_tgi_posterior_mixed(nreps = 250, nsubjects=100,
                                            true_theta_g = exp(-2), sigma_b_g=0.1,
                                            true_theta_s =  exp(0.83), sigma_b_s=0.7,
                                            sigma_eps=0.05, M = 7000, nchains =1)