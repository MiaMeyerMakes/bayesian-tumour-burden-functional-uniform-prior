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

set.seed(12345)

# Define the Stan model with random effects
shrink_stan_jeff = "
  // Stan model for y = exp(-theta*exp(b_i)*x) + epsilon, 
  // where b_i ~ N(0, sigma_b^2) and epsilon ~ N(0, sigma_eps^2)
  data {
    int<lower=0> N;              // number of observations
    int<lower=0> J;              // number of groups/subjects
    vector[N] x;                 // predictor variable
    vector[N] y;                 // response variable
    int<lower=1,upper=J> group[N]; // group/subject indicator for each obs
  }
  
  parameters {
    real<lower=0,upper=5> theta; // decay parameter (with your bounds)
    vector[J] b_raw;             // non-centered parameterization for random effects
    
    real<lower=0> sigma2_epsilon;    // measurement error VARIANCE
    real<lower=0> sigma_b;             // random effects standard deviation
  }
  
  transformed parameters {
    vector[J] b;                 // random effects for each group/subject
    real sigma_epsilon;                // measurement error standard deviation
    
    // Apply non-centered parameterization for computational efficiency
    b = sigma_b * b_raw;
    
    // Convert variance to standard deviation
    sigma_epsilon = sqrt(sigma2_epsilon);
    
  }
  
  model {
    vector[N] mu;                // mean response for each observation
    
    // half-cauchy prior on σ_ε 
    sigma_epsilon ~ cauchy(0,1);
    
    // Half-normal prior on σ_b (can change to half_cauchy if preferred)
    sigma_b ~ normal(0, 2.5);      // half-normal due to <lower=0> constraint
  
    // Normal prior for raw random effects (non-centered parameterization)
    b_raw ~ std_normal();
    
    for (n in 1:N) {
      mu[n] = exp(-theta * exp(b[group[n]]) * x[n]);
    }
    
    // Apply jeffreys prior
    {
      real fisher_info_sum = 0;
      
      // Sum over all subjects i
      for (j in 1:J) {
        real subject_sum = 0;
        
        // Sum over all time points j for subject i
        for (n in 1:N) {
          if (group[n] == j) {  // if observation n belongs to subject i
            real t_ij = x[n];   // time point
            real exp_bi = exp(b[j]);
            real exp_term = exp(-2 * theta * t_ij * exp_bi);
            
            subject_sum += square(t_ij) * square(exp_bi) * exp_term;
          }
        }
        
        fisher_info_sum += subject_sum;
      }
      
      // Add log of Jeffreys prior: log(sqrt(fisher_info_sum)) = 0.5 * log(fisher_info_sum)
      if (fisher_info_sum > 0) {
        target += 0.5 * log(fisher_info_sum);
      }
    }
  
    // Likelihood
    y ~ normal(mu, sigma_epsilon);
  }
  
  generated quantities {
  // for posterior predictive checks
    vector[N] y_pred;            // posterior predictive samples
    vector[N] mu_pred;
    
    for (i in 1:N) {
      mu_pred[i] = exp(-theta * exp(b[group[i]])* x[i]);
      y_pred[i] = normal_rng(mu_pred[i], sigma_epsilon);
    }
  }
"


# Write the Stan model to a file
writeLines(shrink_stan_jeff, "random_shrink_jeff.stan")

generate_dataset <- function(NumReps = 7, J, true_theta=0.18, sigma_b=0.3, sigma_eps=0.03) {
  # Generate some synthetic data for testing
  
  N <- NumReps*J              # total number of observations
  
  
  # Generate random effects for each group
  b <- rnorm(J, 0, sigma_b)
  
  # Generate balanced design with equal observations per group
  group <- rep(1:J, each = N/J)
  
  # Create sequence of x values (could be different for each group too)
  x <- rep(seq(0, 10, length.out = N/J), J)
  
  # Generate y values according to the model with random effects and noise
  y <- numeric(N)
  for (i in 1:N) {
    y[i] <- exp(-true_theta * exp(b[group[i]]) * x[i]) + rnorm(1, 0, sigma_eps)
  }
  
  # Create data frame for ggplot
  sim_data <- data.frame(
    x = x,
    y = y,
    group = factor(group),
    true_curve = NA  # Will be filled later
  )
  
  # Add the true curves for each x value and group
  for (i in 1:nrow(sim_data)) {
    group_i <- sim_data$group[i]
    sim_data$true_curve[i] <- exp(-true_theta * exp(b[as.numeric(group_i)]) * sim_data$x[i])
  }
  
  data_list <- list(
    N = N,
    J = J,
    x = x,
    y = y,
    group = group
  )
  
  return(data_list)
}

simulate_theta_posterior_mixed_jeffreys <- function(nreps = 200, 
                                           nsubjects, 
                                           true_theta=0.18, 
                                           sigma_b, 
                                           sigma_eps, M = 5000) {
  theta_storage <- matrix(nrow = M, ncol = nreps)
  rhat_storage <- numeric(nreps)
  theta_pointest_storage <- numeric(nreps)
  data_storage <- matrix(nrow = nsubjects*7, ncol = nreps)
  
  for (repl in 1:nreps) {
    print(paste("---------------Starting STAN for sample", repl,"-------------------",Sys.time()))
    # generate a dataset
    stan_data_list <- generate_dataset(NumReps = 7, J = nsubjects, true_theta=0.18, sigma_b, sigma_eps)
    data_storage[,repl] <- stan_data_list$y
      
    # get the posterior sample for theta
    # Compile and fit the Stan model
    fit <- stan(
      file = "random_shrink_jeff.stan",
      data = stan_data_list,
      chains = 1,               # number of Markov chains
      iter = M+1000,              # total number of iterations per chain
      warmup = 1000,            # number of warmup iterations per chain
      thin = 1,                 # period for saving samples
      verbose = FALSE,
      control = list(adapt_delta = 0.95) # increase adaptation parameter for better sampling
    )
    
    # store the sample and Rhat values
    stan_summary <- summary(fit)$summary
    theta_rhat <- stan_summary["theta", "Rhat"]
    theta_pointest <- stan_summary["theta", "mean"]
    
    theta_sample <- fit@sim[["samples"]][[1]][["theta"]]
    theta_sample <- theta_sample[1001:(M+1000)]
    theta_storage[,repl] <- theta_sample
    rhat_storage[repl] <- theta_rhat
    theta_pointest_storage[repl] <- theta_pointest
  }
  
  colnames(data_storage) <- paste(rep("dataset", nreps), 1:nreps, sep='')
  filename_data =   paste(format(Sys.time(), "%b-%d-%Hh%M"),"-mixed-effects-shrink-data-JEFF-",nsubjects,".csv",sep="")
  write.csv(as.data.frame(data_storage), filename_data, row.names = FALSE)
  
  colnames(theta_storage) <- paste(rep("sample", nreps), 1:nreps, sep='')
  filename = paste(format(Sys.time(), "%b-%d-%Hh%M"),"-mixed-effects-shrink-posterior-sample-JEFF-",nsubjects,".csv",sep="")
  write.csv(as.data.frame(theta_storage), filename, row.names = FALSE)
  
  return(list(thetas = theta_storage, rhats = rhat_storage, pointests = theta_pointest_storage, datasets = data_storage))
}

thetas_jeffreys_n15 <- simulate_theta_posterior_mixed_jeffreys(nreps=200,nsubjects = 15, M=5000, sigma_b=0.3, sigma_eps = 0.03)
thetas_jeffreys_n50 <- simulate_theta_posterior_mixed_jeffreys(nreps=200,nsubjects = 50, M=5000, sigma_b=0.3, sigma_eps = 0.03)
thetas_jeffreys_n100 <- simulate_theta_posterior_mixed_jeffreys(nreps=100,nsubjects = 100, M=5000, sigma_b=0.3, sigma_eps = 0.03)
