setwd("~/Desktop/Thesis/code/prime")
biom.df.filtered <- readRDS("biom.df.filtered.rds")
biom.datf <- readRDS("biom.df.fitted.rds")

# Comparison plot: Frequentist dashed, Bayesian solid
# Palette & scale helpers ----------------------------------
prior_colors <- c(
  "Observed"    = "#000000",
  "FUP" = "#E69F00", #E69F00 # blue
  "Uniform"= "#718a52",
  "Jeffreys"= "#0072B2"   # orange

)

scale_color_prior <- function(...) {
  scale_color_manual(values = prior_colors, breaks = c("Observed", "FUP", "Uniform", "Jeffreys"), ...)
}

prime_id_subset <- unique(biom.datf$SUBJID)[1:16]

prime_subset <- biom.datf %>%
  filter(SUBJID %in% prime_id_subset) %>%
  mutate(SUBJID = factor(SUBJID, levels = prime_id_subset))  # preserves facet order

comparison_plot <- ggplot(prime_subset, aes(x = BIOMYR)) +
  geom_point(aes(y = BIOMVAL, color = "Observed"), size = 0.75) +
  geom_line(aes(y = BIOMVAL, color = "Observed"), linewidth = 0.4) +
  geom_point(aes(y = pred_fup_stan, color = "FUP"), size = 0.75) +
  geom_line(aes(y = pred_fup_stan, color = "FUP"), linewidth = 0.4) +
  
  geom_point(aes(y = pred_unif_stan, color = "Uniform"), size = 0.75) +
  geom_line(aes(y = pred_unif_stan, color = "Uniform"), linewidth = 0.4) +

  geom_point(aes(y = pred_jeff_stan, color = "Jeffreys"), size = 0.75) +
  geom_line(aes(y = pred_jeff_stan, color = "Jeffreys"), linewidth = 0.4) +
  
  facet_wrap(~ SUBJID, scales = "free_y") +
  labs(
    title = "",
    y = "SLD",
    x = "Time (years)",
    color = "Data"
  ) +
  scale_color_prior() +
  theme_minimal() +
  theme(legend.position = "bottom")

print(comparison_plot)

#------------------------ print the data profiles-------------------------------

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

# Round and convert theta to a factor for better legend handling

plot_data = biom.df.filtered

plot_data$SUBJID <- as.factor(plot_data$SUBJID)


theta_colors <- gradient_factor_palette(plot_data$SUBJID,
                                        low = "#0072B2",
                                        high = "#E69F00",
                                        space = "rgb")

# ggplot(plot_data, aes(x = x, y = y, color = theta)) +
#   geom_line(linewidth = 1) +
#   labs(
#     title = expression("Regression function for different" ~ theta),
#     x = "x",
#     y = expression("exp(-" ~ theta ~ "x)"),
#     color = expression(theta)
#   ) +
#   scale_color_manual(values = theta_colors) +
#   theme_minimal() +
#   theme(legend.position = "bottom")

ggplot(plot_data, aes(x = BIOMYR, y = BIOMVAL, group = SUBJID, color= SUBJID)) +
  geom_point(colour = "grey33", alpha = 0.35, size = 0.9) +
  geom_line(alpha = 0.6,) +
  scale_x_continuous("Year", breaks = 0.5 * (0:5)) +
  scale_y_continuous("SLD (mm)", breaks = c(3, 30, 100, 300)) +
  scale_color_manual(values = theta_colors) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("SLD Over Time: Panitumumab + FOLFOX")
