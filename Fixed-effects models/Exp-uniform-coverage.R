library(pracma) 
library(ggplot2)
library(tidyr)
library(dplyr)


# Set seed for reproducibility
set.seed(123)

prior_density <- function(theta, xupper=10) {
  outside <- exp(-2*theta*xupper)/(4 * theta^3)
  inside <- -2*((theta*xupper)^2)-2*xupper*theta-1+exp(2*theta*xupper)

  # compute and return the density
  density <-sqrt(outside*inside)
  return(density)
}

# Compute the normalizing constant over the valid range 0 < theta < 5
normalizing_constant <- integrate(prior_density, lower = 0.0001, upper = 5)$value

# Define the normalized prior
normalized_prior <- function(theta, xupper=10) {
  prior_density(theta, xupper) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.0001, 5, length.out =1000)

# compute the density values for different x support
prior_sim2 <- normalized_prior(thetas, 2)
prior_sim4 <- normalized_prior(thetas, 4)
prior_sim6 <- normalized_prior(thetas, 6)
prior_sim8 <- normalized_prior(thetas, 8)
prior_sim10 <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, 
                       prior_sim2 = prior_sim2,
                       prior_sim4 = prior_sim4,
                       prior_sim6 = prior_sim6,
                       prior_sim8 = prior_sim8,
                       prior_sim10 = prior_sim10
                       )


# Reshape data to long format (using dplyr::select explicitly)
plot_dat_long <- plot_dat %>%
  dplyr::select(thetas, prior_sim2, prior_sim4, prior_sim6, prior_sim8, prior_sim10) %>%
  pivot_longer(
    cols = -thetas,
    names_to = "prior_type",
    values_to = "density"
  )

# Custom labels for the legend
prior_labels <- c(
  prior_sim2 = "[0,2] support",
  prior_sim4 = "[0,4] support",
  prior_sim6 = "[0,6] support",
  prior_sim8 = "[0,8] support",
  prior_sim10 = "[0,10] support"
)

# Create the plot with legend
ggplot(plot_dat_long, aes(x = thetas, y = density, color = prior_type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    name = NULL,  # Remove legend title if desired
    values = c(
      "prior_sim2" = "maroon",
      "prior_sim4" = "red",
      "prior_sim6" = "coral",
      "prior_sim8" = "orange",
      "prior_sim10" = "yellow2"
    ),
    labels = prior_labels
  ) +
  labs(
    title = "Function uniform prior sensitivity- exponential function",
    x = expression(theta),
    y = expression(p(theta))
  ) +
  theme_minimal()

################################### x in 0,10 #############################################


# Set seed for reproducibility
set.seed(123)

xupper = 10

prior_density_flex <- function(theta, xupper=10) {
  outside <- exp(-2*theta*xupper)/(4 * theta^3)
  inside <- -2*((theta*xupper)^2)-2*xupper*theta-1+exp(2*theta*xupper)
  
  # compute and return the density
  density <-sqrt(outside*inside)
  return(density)
}

# Compute the normalizing constant over the valid range 0 < theta < 5
normalizing_constant <- integrate(prior_density_flex, lower = 0.0001, upper = 5)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex(theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.0001, 5, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.0001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.0001, upper = 5) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

exp_func <- function(x, theta) {exp(-theta*x)}

x <- seq(0, xupper, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- exp_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Regression function for different" ~ theta ~ "with x in [0,10]"),
    x = "x",
    y = expression("exp(-" ~ theta ~ "x)"),
    color = expression(theta)
  ) +
  theme_minimal() +
  scale_color_viridis_d() +  # Use a color scale for better distinction
  theme(legend.position = "bottom")  # Move legend to the bottom


################################### x in 0,7 #############################################


# Set seed for reproducibility
set.seed(123)

xupper = 7

prior_density_flex <- function(theta, xupper=7) {
  outside <- exp(-2*theta*xupper)/(4 * theta^3)
  inside <- -2*((theta*xupper)^2)-2*xupper*theta-1+exp(2*theta*xupper)
  
  # compute and return the density
  density <-sqrt(outside*inside)
  return(density)
}

# Compute the normalizing constant over the valid range 0 < theta < 5
normalizing_constant <- integrate(prior_density_flex, lower = 0.0001, upper = 5)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex(theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.0001, 5, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.0001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.0001, upper = 5) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

exp_func <- function(x, theta) {exp(-theta*x)}

x <- seq(0, xupper, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- exp_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Regression function for different" ~ theta ~ "with x in [0,7]"),
    x = "x",
    y = expression("exp(-" ~ theta ~ "x)"),
    color = expression(theta)
  ) +
  theme_minimal() +
  scale_color_viridis_d() +  # Use a color scale for better distinction
  theme(legend.position = "bottom")  # Move legend to the bottom


################################### x in 0,5 #############################################


# Set seed for reproducibility
set.seed(123)

xupper = 5

prior_density_flex <- function(theta, xupper=5) {
  outside <- exp(-2*theta*xupper)/(4 * theta^3)
  inside <- -2*((theta*xupper)^2)-2*xupper*theta-1+exp(2*theta*xupper)
  
  # compute and return the density
  density <-sqrt(outside*inside)
  return(density)
}

# Compute the normalizing constant over the valid range 0 < theta < 5
normalizing_constant <- integrate(prior_density_flex, lower = 0.0001, upper = 5)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex(theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.0001, 5, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.0001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.0001, upper = 5) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

exp_func <- function(x, theta) {exp(-theta*x)}

x <- seq(0, xupper, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- exp_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Regression function for different" ~ theta ~ "with x in [0,5]"),
    x = "x",
    y = expression("exp(-" ~ theta ~ "x)"),
    color = expression(theta)
  ) +
  theme_minimal() +
  scale_color_viridis_d() +  # Use a color scale for better distinction
  theme(legend.position = "bottom")  # Move legend to the bottom


################################### x in 0,3 #############################################

# Set seed for reproducibility
set.seed(123)

xupper = 3

prior_density_flex <- function(theta, xupper=3) {
  outside <- exp(-2*theta*xupper)/(4 * theta^3)
  inside <- -2*((theta*xupper)^2)-2*xupper*theta-1+exp(2*theta*xupper)
  
  # compute and return the density
  density <-sqrt(outside*inside)
  return(density)
}

# Compute the normalizing constant over the valid range 0 < theta < 5
normalizing_constant <- integrate(prior_density_flex, lower = 0.0001, upper = 5)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex(theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.0001, 5, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.0001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.0001, upper = 5) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

exp_func <- function(x, theta) {exp(-theta*x)}

x <- seq(0, xupper, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- exp_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Regression function for different" ~ theta ~ "with x in [0,3]"),
    x = "x",
    y = expression("exp(-" ~ theta ~ "x)"),
    color = expression(theta)
  ) +
  theme_minimal() +
  scale_color_viridis_d() +  # Use a color scale for better distinction
  theme(legend.position = "bottom")  # Move legend to the bottom

####################################################################################

