library(MASS)  # For multivariate normal
library(RColorBrewer)
library(HDInterval)

# Set seed for reproducibility
set.seed(123)

### Create the exponential sample datasets ######################################

# function to create datasets
sim_emax_datasets <- function(numsets, n, theta0_true, theta1_true, theta2_true, sigma2_true) {
  data <- matrix(nrow = n, ncol = (1+numsets))
  
  x <- seq(0, 4, length.out = n)  # Random x values
  data[,1] <- x
  
  noise_matrix <- matrix(data = rnorm(n*numsets, 0, sqrt(sigma2_true)), nrow = n, ncol = numsets)
  y <- theta0_true + theta1_true*x/(theta2_true+x) + noise_matrix
  data[,2:(numsets+1)] <- y
  colnames(data) <- c('x',apply(matrix(c(rep('sample',numsets), 1:numsets), ncol = 2, byrow = FALSE), 1, function(x) {paste(x[1], x[2], sep="")}))
  data
  return(data)
}

# plot and test on small sample size and numsets size for sense check
testdata <- sim_emax_datasets(numsets=50, n=30, theta0_true=0.2, theta1_true=0.6, theta2_true=0.05, sigma2_true=0.05)
cols = colorRampPalette(brewer.pal(11, "Spectral"))(20)

plot(x = testdata[,1], y = 0.2+0.6*testdata[,1]/(testdata[,1]+ 0.05), type = 'l', col = 'black',
     ylim = c(-2,3), lwd = 3, main = "Simulated datasets from hyperbolic emax model",
     xlab = "x", ylab = "hyberbolic emax")
for (i in 2:20)
{
  lines(x = testdata[,1], y = testdata[,i], col = cols[i])
}

# clear history before we start simulation - keep the data simulating function
rm(list = setdiff(ls(), c("sim_emax_datasets")))

# Simulate data for testing with n = 15, 50 and 100
numsets <- 1000
theta0_true=0.2
theta1_true=0.6
theta2_true=0.05
sigma2_true=0.05

emax_data15 <- sim_emax_datasets(numsets, n=15, theta0_true, theta1_true, theta2_true, sigma2_true)
emax_data50 <- sim_emax_datasets(numsets, n=50, theta0_true, theta1_true, theta2_true, sigma2_true)
emax_data100 <- sim_emax_datasets(numsets, n=100, theta0_true, theta1_true, theta2_true, sigma2_true)

#### MH Sampler #################################################################

# functions to compute priors, posterior and likelihood
log_prior_sigma2 <- function(sig2) {-log(sig2)}

log_prior_theta <- function(theta) {
  theta_prior = sqrt(1/((theta+4)^3 *theta))
  return(log(theta_prior))
}

log_likelihood <- function(theta0=0.2, theta1=0.6, theta2, sigma2, x, y) {
  n = length(x)
  y_obs = theta0+theta1*x/(x+ theta2)
  log_lik = -0.5*n*log(2*pi*sigma2) + ((-0.5/sigma2)*sum((y-y_obs)^2))
  return(log_lik)
}

log_posterior <- function (theta0=0.2, theta1=0.6, theta2, sig2, x, y) {
  log_post = log_likelihood(theta0, theta1, theta2, sig2, x, y) + log_prior_theta(theta2) + log_prior_sigma2(sig2)
  return(log_post)
}

# Metropolis-Hastings algorithm
emax_metropolis_hastings <- function(dataset, num_iterations, burn_in, sigma_theta, sigma_sigma, initial_theta, initial_sigma2) {
  n <- nrow(dataset)
  numsets <- ncol(dataset)-1
  
  # Storage for posterior samples
  samples_sigma2 <- matrix(NA, nrow = num_iterations, ncol = numsets)
  samples_theta <- matrix(NA, nrow = num_iterations, ncol = numsets)
  
  x <- dataset[,1] # x is always the first column of the given dataset
  
  theta_acceptance <- numeric()
  sigma_acceptance <- numeric()
  
  # Generate new chain for each dataset
  for (setnr in 2:ncol(dataset))
  {
    print(paste("Starting MH for sample", setnr))
    # Initialize the chain
    theta_current <- initial_theta
    sigma2_current <- initial_sigma2
    
    # extract the simulated y from dataset
    y <- dataset[,setnr]
    
    theta_acceptance_counter <- 0
    sigma_acceptance_counter <- 0
    # Metropolis-Hastings loop for that sample
    for (iteration in 1:num_iterations) {
      # print(paste("iter", iteration, "theta_curr:", theta_current, ", sig2_curr:", sigma2_current))
      # THETA SAMPLER
      
      # Step 1: Propose a new values for theta
      theta_star <- abs(rnorm(1, theta_current, sigma_theta))
      
      # Step 2: Compute acceptance ratio for theta
      alpha_theta <- log_posterior(theta2=theta_star, sig2=sigma2_current, x=x, y=y)-
        log_posterior(theta2=theta_current, sig2=sigma2_current, x=x, y=y)
      
      # Step 3: Accept or reject theta
      u <- runif(1)
      if (u < exp(alpha_theta))
      {
        theta_current <- theta_star
        theta_acceptance_counter = theta_acceptance_counter + 1
      }
      
      # Step 4: Propose a new value for sigma
      sigma2_star <-  rlnorm(0, 0.1)
      
      # Step 5: Compute acceptance ratio for sigma2
      alpha_sigma_num <- log_posterior(theta2=theta_current, sig2=sigma2_star, x=x, y=y) + log(dlnorm(sigma2_star,0, 0.1))
      alpha_sigma_denom <- log_posterior(theta2=theta_current, sig2=sigma2_current, x=x, y=y)  + log(dlnorm(sigma2_current,0, 0.1))
      alpha_sigma <- alpha_sigma_num - alpha_sigma_denom
      
      # Step 6: Accept or reject sigma2
      u <- runif(1)
      if (u < exp(alpha_sigma)) {
        sigma2_current <- sigma2_star
        sigma_acceptance_counter = sigma_acceptance_counter + 1
      }
      
      # Store the values
      samples_sigma2[iteration,setnr-1] <- sigma2_current
      samples_theta[iteration,setnr-1] <- theta_current
    }
    theta_acceptance <- append(theta_acceptance, sum(theta_acceptance_counter)/num_iterations)
    sigma_acceptance <- append(sigma_acceptance, sum(sigma_acceptance_counter)/num_iterations)
  }
  
  # Plot the first sample's chain for sense check (plot with burn in)
  plot(samples_theta[,1], type = 'l', col = 'blue', xlab = 'Iteration', ylab = 'Theta', main = 'Trace plot of Theta')
  plot(samples_sigma2[,1], type = 'l', col = 'red', xlab = 'Iteration', ylab = 'Sigma^2', main = 'Trace plot of Sigma^2')
  
  print(paste("Theta acceptance:", mean(theta_acceptance)))
  print(paste("Sigma acceptance:", mean(sigma_acceptance)))
  
  # Return MH samples after removing burn-in
  return(list(sigma2_MH = samples_sigma2[(burn_in+1):(num_iterations),], theta_MH = samples_theta[(burn_in+1):(num_iterations),]))
}

# Parameters for the Metropolis-Hastings algorithm
num_iterations <- 10000  # Total number of iterations
burn_in <- 1000          # Burn-in period
sigma_theta <- 0.05      # Step size for theta
sigma_sigma <- 0.1      # Step size for log(sigma2)
initial_theta <- 1    # Initial value for theta
initial_sigma2 <- 1     # Initial value for sigma^2

# Run the Metropolis-Hastings sampler
emax_n15_samples_list <- emax_metropolis_hastings(emax_data15, num_iterations, burn_in, sigma_theta, sigma_sigma, initial_theta, initial_sigma2)
emax_n50_samples_list <- emax_metropolis_hastings(emax_data50, num_iterations, burn_in, sigma_theta, sigma_sigma, initial_theta, initial_sigma2)
emax_n100_samples_list <- emax_metropolis_hastings(emax_data100, num_iterations, burn_in, sigma_theta, sigma_sigma, initial_theta, initial_sigma2)

# # Save samples to prevent unnecessary re-runs
# write.csv(emax_n15_samples_list[["theta_MH"]], "emax_theta_posterior_n15.csv", row.names = FALSE, col.names = FALSE)
# write.csv(emax_n50_samples_list[["theta_MH"]], "emax_theta_posterior_n50.csv", row.names = FALSE, col.names = FALSE)
# write.csv(emax_n100_samples_list[["theta_MH"]], "emax_theta_posterior_n100.csv", row.names = FALSE, col.names = FALSE)
# 
# n15_samplez <- read.csv("emax_theta_posterior_n15.csv", header = TRUE)
# n50_samplez <- read.csv("emax_theta_posterior_n50.csv", header = TRUE)
# n100_samplez <- read.csv("emax_theta_posterior_n100.csv", header = TRUE)

# Function to calculate the coverage probabilities
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
  plot(thetas[,1], type = 'l', col = 'blue', xlab = 'Iteration', ylab = 'Theta', main = 'Trace plot of Theta')
  abline(h=hdi_interval90[1,1], col = 'green')
  abline(h=hdi_interval90[2,1], col = 'green')
  abline(h=hdi_interval95[1,1], col = 'red')
  abline(h=hdi_interval95[2,1], col = 'red')
  
  ILE90 <- mean(hdi_interval90[2,] - hdi_interval90[1,])
  ILE95 <- mean(hdi_interval95[2,] - hdi_interval95[1,])
  
  posterior_means <- apply(thetas, 2, mean)
  posterior_median <- apply(thetas, 2, median)
  
  MAE_mean <- mean(abs(posterior_means-theta_true))
  MAE_median <- mean(abs(posterior_median-theta_true))
  
  results <- c(CP90, CP95, ILE90, ILE95, MAE_mean, MAE_median)
  return(results)
}

emax_stats15 = statistics_calc(theta2_true, emax_n15_samples_list)
emax_stats50 = statistics_calc(theta2_true, emax_n50_samples_list)
emax_stats100 = statistics_calc(theta2_true, emax_n100_samples_list)

emax_statistics_table = rbind(emax_stats15, emax_stats50, emax_stats100)
row.names(emax_statistics_table) = c("n=15", "n=50", "n=100")
colnames(emax_statistics_table) = c("CP_90", "CP_95", "ILE90", "ILE95", "MAE_mean", "MAE_median")
emax_statistics_table
#  emax_f  
#       CP_90 CP_95      ILE90     ILE95   MAE_mean MAE_median
# n=15  0.973 0.991 0.28730129 0.3762861 0.08593298 0.06441821
# n=50  0.877 0.952 0.12467992 0.1526932 0.03316077 0.03289978
# n=100 0.851 0.929 0.09163402 0.1098929 0.02419880 0.02423747

