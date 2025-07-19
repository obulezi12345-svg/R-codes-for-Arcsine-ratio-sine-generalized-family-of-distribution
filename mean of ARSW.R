

setwd("C:\\Users\\1012 G2\\Documents\\R files\\ARS-W")
library(e1071)
library(plotly)
library(fields)
library(Cairo)
library(RColorBrewer)
library(scatterplot3d)
library(viridis)


dARSW <- function(par, x) {
  theta=par[1]
  lambda=par[2]
  omega=par[3]
  if (theta <= 0 || lambda <= 0 || omega <= 0) {
    stop("Parameters theta, lambda, and omega must be positive.")
  }
  if (any(x < 0)) {
    warning("x must be non-negative. Returning 0 for negative x values.")
  }
  
  # Handle x < 0 separately to ensure PDF is 0
  pdf_values <- numeric(length(x))
  valid_indices <- which(x >= 0)
  
  if (length(valid_indices) > 0) {
    x_valid <- x[valid_indices]
    
    G_x <- 1 - exp(-lambda * x_valid^omega)
    # Baseline Weibull PDF g(x)
    g_x <- ifelse(x_valid == 0 & omega < 1, Inf, # Handle g(0) for omega < 1
                  lambda * omega * x_valid^(omega - 1) * exp(-lambda * x_valid^omega))
    
    sin_G_x <- sin(pi / 2 * G_x)
    cos_G_x <- cos(pi / 2 * G_x)
    
    numerator <- pi * theta * cos_G_x * g_x
    denominator <- (theta + (2 - theta) * sin_G_x)^2
    
    pdf_values[valid_indices] <- numerator / denominator
  }
  
  # Ensure PDF is 0 for x = 0 when omega > 1 (due to x^(omega-1) term)
  # and handle potential NaN/Inf from division by zero or invalid inputs
  pdf_values[is.nan(pdf_values) | is.infinite(pdf_values)] <- 0
  
  return(pdf_values)
}






# Define the function for skewness calculation
skewness_calculation <- function(data) {
  mean_val <- mean(data)
  sd_val <- sd(data)
  skew <- mean((data - mean_val)^3) / sd_val^3
  return(skew)
}

# Function to generate data from the distribution
generate_data <- function(n, par) {
  set.seed(123)  # For reproducibility
  x_vals <- runif(n, min = 0.1, max = 1)  # Generate x values within a range
  data <- sapply(x_vals, function(x) dARSW(par, x))  # Calculate density values
  return(data)
}


# Create grids of parameter values
theta_vals <- seq(0.05, 1.5, length.out = 50)
lambda_vals <- seq(0.05, 1.5, length.out = 50)
omega_vals <- seq(0.05, 1.5, length.out = 50)


Z_mean <- matrix(NA, nrow = length(theta_vals), ncol = length(lambda_vals))

# Loop through parameter combinations and compute mean for each combination
for (i in 1:length(theta_vals)) {
  for (j in 1:length(lambda_vals)) {
    # Generate data from the distribution for each parameter combination
    simulated_data <- generate_data(1000, c(theta_vals[i], lambda_vals[j], omega_vals[i]))
    # Compute mean for the simulated data
    Z_mean[i, j] <- mean(simulated_data)
  }
}
#alpha = \u03B1
#beta = \u03B2
#sigma = \u03C3

# Plotting the 3D surface using persp() for skewness
windows()
par(mfrow=c(1,1))
# Open the Cairo PDF device
cairo_pdf("mean_2.pdf")

# Create a color palette
library(RColorBrewer)
library(scatterplot3d)

library(viridis)

red_darkorange_palette <- colorRampPalette(c("red", "darkorange"))
# Red to Purple hex codes

# Generate the color sequence
colors <- red_darkorange_palette(length(omega_vals))


# Create the 3D surface plot for mean values using persp()
persp(theta_vals, lambda_vals, Z_mean, theta = 30, phi = 30, 
      col = colors, xlab = expression("\u03B8"), ylab = expression("\u03BB"), zlab = "mean",
      main = "", axes = TRUE, ticktype='detailed')
# Calculate zfacet values
zfacet <- Z_mean[-1, -1]
facetcol <- cut(zfacet, length(colors))

# Add a color bar legend using image.plot() from the fields package
image.plot(legend.only = TRUE, zlim = range(zfacet, na.rm = TRUE), col = colors, horizontal = FALSE)

# Close the Cairo PDF device
dev.off()