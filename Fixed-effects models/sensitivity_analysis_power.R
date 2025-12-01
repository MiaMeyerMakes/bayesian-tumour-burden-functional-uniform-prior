library(pracma) 
library(ggplot2)
library(tidyr)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

prior_density <- function(theta, xupper=1) {
  outside <- xupper^(2*theta + 1)
  bottom <- 1/(2*theta + 1)
  lnlam <- log(xupper)
  inner <- bottom*(lnlam^2) + 2*lnlam*bottom^2-xupper*bottom^3

  # compute and return the density
  density <-sqrt(outside*inner)
  return(density)
}

# Compute the normalizing constant over the valid range 0 < theta < 5
normalizing_constant <- integrate(prior_density, lower = 0.0501, upper = 20)$value

# Define the normalized prior
normalized_prior <- function(theta, xupper=1) {
  prior_density(theta, xupper) / normalizing_constant
}

# generate values for plotting
thetas <- seq(0.0001, 5, length.out =1000)

# compute the density values for different x support
prior_sim2 <- normalized_prior(thetas, 0.2)
prior_sim4 <- normalized_prior(thetas, 0.4)
prior_sim6 <- normalized_prior(thetas, 0.6)
prior_sim8 <- normalized_prior(thetas, 0.8)
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
  prior_sim2 = "[0,0.2] support",
  prior_sim4 = "[0,0.4] support",
  prior_sim6 = "[0,0.6] support",
  prior_sim8 = "[0,0.8] support",
  prior_sim10 = "[0,1] support"
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
    title = "Function uniform prior sensitivity- power function",
    x = expression(theta),
    y = expression(p(theta))
  ) +
  theme_minimal()