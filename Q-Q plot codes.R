# Set working directory (optional, adjust path as needed)
setwd("C:\\Users\\1012 G2\\Documents\\R files\\ARS-W") # Uncomment and adjust if needed

# Suppress warnings for cleaner output during execution,
# but be mindful of potential issues this might hide.
options(warn = -1)

# Load necessary libraries
library(AdequacyModel) # For goodness.fit (though not used directly for QQ plot here)
library(fitdistrplus) # Useful for distribution fitting, but mle2 is used here
library(bbmle) # For Maximum Likelihood Estimation (mle2)
library(ggplot2) # For enhanced plotting, though base R graphics are used for qqplot
library(psych) # General purpose psychology tools, not strictly needed for this
library(pracma) # Practical numerical math functions, not strictly needed for this
library(wesanderson) # Color palettes, not used in this specific plot
# library(kolmim) # Kolmogorov-Smirnov test, not used here
# library(MASS) # Modern Applied Statistics with S, not used here
# library(muhaz) # Hazard estimation, not used here

# Probability Density Function (PDF) for ARSW distribution
# par: vector of parameters (theta, lambda, omega)
# x: data points
dARSW <- function(par, x) {
  theta = par[1]
  lambda = par[2]
  omega = par[3]
  
  # Return 0 for invalid parameters to guide optimization away from these regions
  if (theta <= 0 || lambda <= 0 || omega <= 0) {
    return(rep(0, length(x)))
  }
  
  # Initialize PDF values to 0
  pdf_values <- numeric(length(x))
  # Identify valid indices where x is non-negative
  valid_indices <- which(x >= 0)
  
  if (length(valid_indices) > 0) {
    x_valid <- x[valid_indices]
    
    # G(x) is the CDF of the Weibull distribution
    G_x <- 1 - exp(-lambda * x_valid^omega)
    
    # g(x) is the PDF of the Weibull distribution
    # Handle x_valid = 0 specifically for omega < 1 to avoid NaN/Inf issues
    g_x <- ifelse(x_valid == 0 & omega < 1, Inf,
                  lambda * omega * x_valid^(omega - 1) * exp(-lambda * x_valid^omega))
    
    sin_G_x <- sin(pi / 2 * G_x)
    cos_G_x <- cos(pi / 2 * G_x)
    
    numerator <- pi * theta * cos_G_x * g_x
    denominator <- (theta + (2 - theta) * sin_G_x)^2
    
    # Calculate PDF values for valid indices
    pdf_values[valid_indices] <- numerator / denominator
  }
  
  # Ensure PDF values are non-negative and handle potential NaN/Inf
  pdf_values[is.nan(pdf_values) | is.infinite(pdf_values) | pdf_values < 0] <- 1e-300 # Small positive number instead of 0 for log-likelihood stability
  
  return(pdf_values)
}

# Cumulative Distribution Function (CDF) for ARSW distribution
# par: vector of parameters (theta, lambda, omega)
# x: data points (renamed from 'q' to 'x' for consistency with dARSW and goodness.fit)
pARSW = function(par, x) {
  theta = par[1]
  lambda = par[2]
  omega = par[3]
  
  # Return 0 or a very small number for invalid parameters
  if (theta <= 0 || lambda <= 0 || omega <= 0) {
    warning("Parameters theta, lambda, and omega must be positive for pARSW. Returning 0.")
    return(rep(0, length(x)))
  }
  
  # Handle negative x values: CDF is 0 for x < 0
  if (any(x < 0)) {
    warning("x must be non-negative. Returning 0 for negative x values.")
  }
  
  # Initialize CDF values to 0
  cdf_values <- numeric(length(x))
  # Identify valid indices where x is non-negative
  valid_indices <- which(x >= 0)
  
  if (length(valid_indices) > 0) {
    x_valid <- x[valid_indices]
    
    # G(x) is the CDF of the Weibull distribution
    G_x <- 1 - exp(-lambda * x_valid^omega)
    sin_G_x <- sin(pi / 2 * G_x)
    
    numerator_frac <- theta * (1 - sin_G_x)
    denominator_frac <- theta + (2 - theta) * sin_G_x
    fraction_term <- numerator_frac / denominator_frac
    
    # Calculate CDF values for valid indices
    cdf_values[valid_indices] <- 1 - fraction_term
  }
  
  # Ensure values are within [0, 1] due to potential numerical inaccuracies
  cdf_values[cdf_values < 0] <- 0
  cdf_values[cdf_values > 1] <- 1
  
  return(cdf_values)
}

# Negative Log-Likelihood function for ARSW distribution
# This function is minimized by mle2 to find the best parameters
LL_ARSW = function(theta, lambda, omega) {
  # Return Inf if parameters are invalid (non-positive)
  if (theta <= 0 || lambda <= 0 || omega <= 0) {
    return(Inf)
  }
  
  # Calculate PDF values for the given data x and parameters
  pdf_vals <- dARSW(c(theta, lambda, omega), x)
  
  # Return Inf if any PDF value is non-positive (log(0) is -Inf, sum becomes -Inf)
  if (any(pdf_vals <= 0)) {
    return(Inf)
  }
  
  # Return the negative sum of log-likelihoods
  -sum(log(pdf_vals))
}

# Quantile Function (Inverse CDF) for ARSW distribution
# p: probabilities (between 0 and 1)
# par: vector of parameters (theta, lambda, omega)
# lower_bound, upper_bound: search interval for uniroot
qARSW <- function(p, par, lower_bound = 1e-6, upper_bound = 10000) { # Increased upper_bound significantly
  if (any(p < 0 | p > 1)) stop("p must be between 0 and 1")
  
  theta <- par[1]
  lambda <- par[2]
  omega <- par[3]
  
  sapply(p, function(prob) {
    quantile_func_for_uniroot <- function(val) pARSW(par, val) - prob
    
    # Try to find a suitable upper bound if the initial one doesn't bracket
    current_upper <- upper_bound
    f_lower <- quantile_func_for_uniroot(lower_bound)
    f_upper <- quantile_func_for_uniroot(current_upper)
    
    # Loop to expand upper bound if needed, with a safety limit
    # This helps uniroot find a root even if the initial interval is too small
    while ((is.nan(f_lower) || is.nan(f_upper) || is.infinite(f_lower) || is.infinite(f_upper) || (f_lower * f_upper > 0)) && current_upper < 1e5) { # Even higher cap for expansion
      current_upper <- current_upper * 2
      f_upper <- quantile_func_for_uniroot(current_upper)
    }
    
    # If still no bracketing after expansion, return NA
    if ((is.nan(f_lower) || is.nan(f_upper) || is.infinite(f_lower) || is.infinite(f_upper) || (f_lower * f_upper > 0))) {
      warning(paste("qARSW: No root found in search interval for p=", prob, ". Returning NA.", sep=""))
      return(NA)
    }
    
    # Find the root safely using tryCatch for error handling
    tryCatch(
      uniroot(quantile_func_for_uniroot, lower = lower_bound, upper = current_upper, tol = 1e-8)$root,
      error = function(e) {
        warning(paste("qARSW: Error in uniroot for p=", prob, ": ", e$message, ". Returning NA.", sep=""))
        return(NA)
      }
    )
  })
}

# Data-III (your provided dataset)
x = c(3.7, 2.74, 2.73, 2.5, 3.6, 3.11, 3.27, 2.87, 1.47, 3.11, 4.42, 2.41, 3.19, 3.22, 1.69, 3.28, 3.09, 1.87, 3.15, 4.9,
      3.75, 2.43, 2.95, 2.97, 3.39, 2.96, 2.53, 2.67, 2.93, 3.22, 3.39, 2.81, 4.2, 3.33, 2.55, 3.31, 3.31, 2.85, 2.56, 3.56,
      3.15, 2.35, 2.55, 2.59, 2.38, 2.81, 2.77, 2.17, 2.83, 1.92, 1.41, 3.68, 2.97, 1.36, 0.98, 2.76, 4.91, 3.68, 1.84, 1.59,
      3.19, 1.57, 0.81, 5.56, 1.73, 1.59, 2, 1.22, 1.12, 1.71, 2.17, 1.17, 5.08, 2.48, 1.18, 3.51, 2.17, 1.69, 1.25, 4.38,
      1.84, 0.39, 3.68, 2.48, 0.85, 1.61, 2.79, 4.7, 2.03, 1.8, 1.57, 1.08, 2.03, 1.61, 2.12, 1.89, 2.88, 2.82, 2.05, 3.65)

# --- Parameter Estimation ---
# Initial parameter guesses (crucial for optimization)
# These values are chosen based on typical ranges or prior knowledge from your tables
start_params <- list(theta = 0.5, lambda = 1.0, omega = 1.0)

# Perform Maximum Likelihood Estimation using mle2
# method = "L-BFGS-B" allows for lower and upper bounds on parameters,
# which is essential for ensuring positive parameters.
# The 'data' argument passes 'x' to the LL_ARSW function's environment.
fit_ARSW <- mle2(LL_ARSW, start = start_params, data = list(x = x),
                 method = "L-BFGS-B",
                 lower = c(theta = 1e-6, lambda = 1e-6, omega = 1e-6))

# Extract the estimated parameters (MLEs)
mle_params <- coef(fit_ARSW)
print("Estimated Parameters (MLEs):")
print(mle_params)

# Check for convergence issues in MLE
if (fit_ARSW@details$convergence != 0) {
  warning("MLE optimization did not converge successfully. This might affect the Q-Q plot accuracy.")
  print(fit_ARSW@details$message)
}

# Ensure MLE parameters are finite before passing to qARSW
if (any(is.na(mle_params)) || any(is.infinite(mle_params))) {
  stop("Estimated parameters (MLEs) are not finite. Cannot proceed with Q-Q plot. Review MLE convergence.")
}

# --- Q-Q Plot Generation ---
# 1. Sort the empirical data
data_sorted <- sort(x)

# 2. Generate probabilities for theoretical quantiles
# These probabilities correspond to the positions of the sorted data points
# (i - 0.5) / n is a common plotting position formula
p_values <- (1:length(data_sorted) - 0.5) / length(data_sorted)

# 3. Calculate theoretical quantiles using the estimated parameters and qARSW function
message("Calculating theoretical quantiles...")
theoretical_quantiles <- qARSW(p_values, mle_params)

# Debugging: Print theoretical_quantiles to inspect for non-finite values
message("Theoretical Quantiles before filtering (first 10):")
print(head(theoretical_quantiles, 10))
message("Theoretical Quantiles before filtering (last 10):")
print(tail(theoretical_quantiles, 10))

# Filter out any non-finite values (NA, Inf, -Inf) from theoretical quantiles
# data_sorted is guaranteed to be finite from the input `x`
valid_theoretical_indices <- is.finite(theoretical_quantiles)

data_sorted_filtered <- data_sorted[valid_theoretical_indices]
theoretical_quantiles_filtered <- theoretical_quantiles[valid_theoretical_indices]

# IMPORTANT: Check if filtered data is empty before plotting
if (length(data_sorted_filtered) == 0) {
  stop("After filtering, no finite data points remain for plotting the Q-Q plot. This means all theoretical quantiles were non-finite. Review the ARSW distribution's definition, estimated parameters, and qARSW function's behavior.")
}

# Check if the ranges are effectively zero, which can also cause 'xlim'/'ylim' errors
if (diff(range(data_sorted_filtered)) == 0) {
  stop("Empirical quantiles have zero range after filtering. Cannot create Q-Q plot. This should not happen if original data 'x' is varied.")
}
if (diff(range(theoretical_quantiles_filtered)) == 0) {
  stop("Theoretical quantiles have zero range after filtering. Cannot create Q-Q plot. This suggests the quantile function is returning a constant value for all probabilities.")
}

# 4. Create the Q-Q plot
# windows() opens a new graphics device (useful on Windows OS)
windows()

qqplot(data_sorted_filtered, theoretical_quantiles_filtered,
       pch = 17, col = "red",
       xlab = "Empirical Quantiles (Data-III)",
       ylab = "ARSW Theoretical Quantiles",
       main = "",
       cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)

# Add the 45-degree reference line (y = x)
abline(0, 1, col = "black", lwd = 2)

# Add a legend
legend("topleft", legend = c("Data Quantiles", "Reference Line"),
       col = c("red", "black"), pch = c(17, NA), lty = c(NA, 1), lwd = c(NA, 2),
       bty = "n")

# You can also add a text label for the estimated parameters on the plot if desired
# text(x = max(data_sorted_filtered) * 0.7, y = min(theoretical_quantiles_filtered) * 1.3,
#      labels = paste("MLE: \n",
#                     "theta = ", round(mle_params[1], 4), "\n",
#                     "lambda = ", round(mle_params[2], 4), "\n",
#                     "omega = ", round(mle_params[3], 4)),
#      adj = c(0, 0), cex = 0.9)
