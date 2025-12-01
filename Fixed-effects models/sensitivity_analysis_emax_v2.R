
# setwd("\\\\sunrga.stb.sun.ac.za/home/ebw/22675760/Downloads/THESIS")

library(pracma)
library(ggplot2)
library(tidyr)
library(dplyr)

set.seed(123)

gradient_factor_palette <- function(fct,
                                    low = "#0072B2",
                                    high = "#E69F00",
                                    space = c("rgb", "hcl")) {
  space <- match.arg(space)
  lvls <- levels(fct)
  n <- length(lvls)
  if (space == "hcl") {
    if (!requireNamespace("colorspace", quietly = TRUE)) {
      warning("colorspace not installed; using RGB interpolation instead.")
      cols <- grDevices::colorRampPalette(c(low, high))(n)
    } else {
      cols <- colorspace::colorRampPalette(c(low, high), space = "HCL")(n)
    }
  } else {
    cols <- grDevices::colorRampPalette(c(low, high))(n)
  }
  setNames(cols, lvls)
}


# prior_density <- function(theta, xupper) {
#   A_val <- (1/(3*theta))-(3*xupper^2+3*xupper*theta + theta^2)/(3*(xupper+theta)^3)
# 
#   # compute and return the density
#   density <-sqrt(A_val)
#   return(density)
# }
# 
# # Compute the normalizing constant over the valid range 0.004 < theta < 6
# normalizing_constant <- integrate(prior_density, lower = 0.004001, upper = 6)$value
# 
# # Define the normalized prior
# normalized_prior <- function(theta, xupper) {
#   prior_density(theta, xupper) / normalizing_constant
# }
# 
# # generate values for plotting
# thetas <- seq(0.004001, 6, length.out =1000)
# 
# # compute the density values for different x support
# prior_sim1 <- normalized_prior(thetas, 1)
# prior_sim2 <- normalized_prior(thetas, 2)
# prior_sim3 <- normalized_prior(thetas, 3)
# prior_sim4 <- normalized_prior(thetas, 4)
# 
# plot_dat <- data.frame(thetas = thetas, 
#                        prior_sim1 = prior_sim1,
#                        prior_sim2 = prior_sim2,
#                        prior_sim3 = prior_sim3,
#                        prior_sim4 = prior_sim4
#                        )
# 
# # plot the density function using ggplot2
# ggplot(plot_dat, aes(x = thetas)) +
#   geom_line(aes(y = prior_sim1), color = "maroon", linewidth = 1) +
#   geom_line(aes(y = prior_sim2), color = "red", linewidth = 1) +
#   geom_line(aes(y = prior_sim3), color = "coral", linewidth = 1) +
#   geom_line(aes(y = prior_sim4), color = "goldenrod", linewidth = 1) +
#   labs(
#     title = "Function uniform prior sensitivity- hyperbolic emax-type function",
#     x = expression(theta),
#     y = expression(p(theta))
#   ) +
#   theme_minimal()
# 
# # First, reshape your data to long format
# library(tidyr)
# library(dplyr)
# 
# # Reshape data to long format (using dplyr::select explicitly)
# plot_dat_long <- plot_dat %>%
#   dplyr::select(thetas, prior_sim1, prior_sim2, prior_sim3, prior_sim4) %>%
#   pivot_longer(
#     cols = -thetas,
#     names_to = "prior_type",
#     values_to = "density"
#   )
# 
# # Custom labels for the legend
# prior_labels <- c(
#   prior_sim1 = "1 support upper limit",
#   prior_sim2 = "2 support upper limit",
#   prior_sim3 = "3 support upper limit",
#   prior_sim4 = "4 support upper limit"
# )
# 
# # Create the plot with legend
# ggplot(plot_dat_long, aes(x = thetas, y = density, color = prior_type)) +
#   geom_line(linewidth = 1) +
#   scale_color_manual(
#     name = NULL,  # Remove legend title if desired
#     values = c(
#       "prior_sim1" = "maroon",
#       "prior_sim2" = "red",
#       "prior_sim3" = "coral",
#       "prior_sim4" = "goldenrod"
#     ),
#     labels = prior_labels
#   ) +
#   labs(
#     title = "Function uniform prior sensitivity- hyperbolic emax-type function",
#     x = expression(theta),
#     y = expression(p(theta))
#   ) +
#   theme_minimal()
# 
# prior_density_flex <- function(theta, xupper) {
#   A_val <- (1/(3*theta))-(3*xupper^2+3*xupper*theta + theta^2)/(3*(xupper+theta)^3)
#   
#   # compute and return the density
#   density <-sqrt(A_val)
#   return(density)
# }

################################### x in 0,4 #############################################

# Set seed for reproducibility
set.seed(123)

xupper4 = 4

prior_density_flex4 <- function(theta, xupper=4) {
  A_val <- (1/(3*theta))-(3*xupper^2+3*xupper*theta + theta^2)/(3*(xupper+theta)^3)
  
  # compute and return the density
  density <-sqrt(A_val)
  return(density)
}

# Compute the normalizing constant over the valid range 0.004 < theta < 6
normalizing_constant <- integrate(prior_density_flex4, lower = 0.004001, upper = 6)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex4(theta=theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.004001, 6, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.004001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.004001, upper = 6) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

emax_func <- function(x, theta) {0.2+0.6*x/(x+theta)}

x <- seq(0, xupper4, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- emax_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)


theta_colors <- gradient_factor_palette(plot_data$theta,
                                        low = "#0072B2",
                                        high = "#E69F00",
                                        space = "rgb")

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Emax regression function for different" ~ theta ~ "with x in [0,4]"),
    x = "x",
    y = expression("0.2+0.6x/(x+" ~ theta ~ ")"),
    color = expression(theta)
  ) +
  scale_color_manual(values = theta_colors) +
  theme_minimal()


################################### x in 0,3 #############################################

# Set seed for reproducibility
set.seed(123)

xupper3 = 3

prior_density_flex3 <- function(theta, xupper=3) {
  A_val <- (1/(3*theta))-(3*xupper^2+3*xupper*theta + theta^2)/(3*(xupper+theta)^3)
  
  # compute and return the density
  density <-sqrt(A_val)
  return(density)
}

# Compute the normalizing constant over the valid range 0.004 < theta < 6
normalizing_constant <- integrate(prior_density_flex3, lower = 0.004001, upper = 6)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex3(theta=theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.004001, 6, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.004001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.004001, upper = 6) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

emax_func <- function(x, theta) {0.2+0.6*x/(x+theta)}

x <- seq(0, xupper3, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- emax_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)


theta_colors <- gradient_factor_palette(plot_data$theta,
                                        low = "#0072B2",
                                        high = "#E69F00",
                                        space = "rgb")

# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Emax regression function for different" ~ theta ~ "with x in [0,3]"),
    x = "x",
    y = expression("0.2+0.6x/(x+" ~ theta ~ ")"),
    color = expression(theta)
  ) +
  scale_color_manual(values = theta_colors) +
  theme_minimal()


################################### x in 0,2 #############################################

# Set seed for reproducibility
set.seed(123)

xupper2 = 2

prior_density_flex2 <- function(theta, xupper=2) {
  A_val <- (1/(3*theta))-(3*xupper^2+3*xupper*theta + theta^2)/(3*(xupper+theta)^3)
  
  # compute and return the density
  density <-sqrt(A_val)
  return(density)
}

# Compute the normalizing constant over the valid range 0.004 < theta < 6
normalizing_constant <- integrate(prior_density_flex2, lower = 0.004001, upper = 6)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex2(theta=theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.004001, 6, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.004001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.004001, upper = 6) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

emax_func <- function(x, theta) {0.2+0.6*x/(x+theta)}

x <- seq(0, xupper2, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- emax_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)


theta_colors <- gradient_factor_palette(plot_data$theta,
                                        low = "#0072B2",
                                        high = "#E69F00",
                                        space = "rgb")
# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Emax regression function for different" ~ theta ~ "with x in [0,2]"),
    x = "x",
    y = expression("0.2+0.6x/(x+" ~ theta ~ ")"),
    color = expression(theta)
  ) +
  scale_color_manual(values = theta_colors) +
  theme_minimal()




################################### x in 0,1 #############################################

# Set seed for reproducibility
set.seed(123)

xupper1 = 1

prior_density_flex1 <- function(theta, xupper=1) {
  A_val <- (1/(3*theta))-(3*xupper^2+3*xupper*theta + theta^2)/(3*(xupper+theta)^3)
  
  # compute and return the density
  density <-sqrt(A_val)
  return(density)
}

# Compute the normalizing constant over the valid range 0.004 < theta < 6
normalizing_constant <- integrate(prior_density_flex1, lower = 0.004001, upper = 6)$value

# Define the normalized prior
normalized_prior <- function(theta) {
  prior_density_flex1(theta=theta) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.004001, 6, length.out =1000)

# compute the density values
prior_sim <- normalized_prior(thetas)

plot_dat <- data.frame(thetas = thetas, prior_sim = prior_sim)


# Integrate to get the cdf
prior_cdf <- function(x) {integrate(normalized_prior, lower = 0.004001, upper = x)$value}

# get the percentiles:
prior_percentile <- function(p, prior_cdf, lower = 0.004001, upper = 6) 
{uniroot(function(x) prior_cdf(x) - p, lower = lower, upper = upper)$root}

# calculate plot specific percentiles
percs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9999999)  

percentile_values <- sapply(percs, prior_percentile, prior_cdf = prior_cdf)
percentile_values
# 0.00010 0.02220 0.04857 0.08105 0.12277 0.17966 0.26462 0.40936 0.70536 1.49104 5.00000

# Define my nonlinear function

emax_func <- function(x, theta) {0.2+0.6*x/(x+theta)}

x <- seq(0, xupper1, length.out = 100)

# Create an empty data frame to store all plotting data
plot_data <- data.frame(x = numeric(), y = numeric(), theta = numeric())

# Loop through each theta value and compute the corresponding y values
for (theta in percentile_values) {
  y <- emax_func(x, theta)
  plot_data <- rbind(plot_data, data.frame(x = x, y = y, theta = theta))
}

# Round and convert theta to a factor for better legend handling
plot_data$theta <- round(plot_data$theta, 4)
plot_data$theta <- as.factor(plot_data$theta)


theta_colors <- gradient_factor_palette(plot_data$theta,
                                        low = "#0072B2",
                                        high = "#E69F00",
                                        space = "rgb")
# Plot using ggplot2
ggplot(plot_data, aes(x = x, y = y, color = theta)) +
  geom_line(linewidth = 1) +
  labs(
    title = expression("Emax regression function for different" ~ theta ~ "with x in [0,1]"),
    x = "x",
    y = expression("0.2+0.6x/(x+" ~ theta ~ ")"),
    color = expression(theta)
  ) +
  scale_color_manual(values = theta_colors) +
  theme_minimal()

