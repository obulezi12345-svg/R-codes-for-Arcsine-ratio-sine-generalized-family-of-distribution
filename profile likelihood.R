# Set working directory (optional, adjust path as needed)
setwd("C:\\Users\\1012 G2\\Documents\\R files\\ARS-W") # Uncomment and adjust if needed

# Suppress warnings for cleaner output during execution,
# but be mindful of potential issues this might hide.
options(warn = -1)

# Load necessary libraries
library(bbmle) # For Maximum Likelihood Estimation (mle2)
library(ggplot2) # For enhanced plotting

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

# Data-I (your provided dataset)
x = c(0.2,0.3, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.7, 0.8, 0.8, 1.0,1.0, 1.0, 1.0, 1.1, 1.3, 1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.2, 2.5,2.7, 3.0, 3.0, 3.3, 3.3, 4.0, 4.0, 4.5, 4.7, 5.0, 5.4, 5.4, 7.0,7.5, 8.8, 9.0, 10.3, 22.0, 24)

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
  warning("MLE optimization did not converge successfully. This might affect the accuracy of the profile plots.")
  print(fit_ARSW@details$message)
}

# Ensure MLE parameters are finite before proceeding
if (any(is.na(mle_params)) || any(is.infinite(mle_params))) {
  stop("Estimated parameters (MLEs) are not finite. Cannot proceed with profile plots. Review MLE convergence.")
}

# --- Likelihood Profile Plots ---

# Function to plot the likelihood profile for a given parameter
plot_likelihood_profile <- function(param_name, mle_fit, data_x) {
  mle_params <- coef(mle_fit)
  
  param_mle <- mle_params[param_name]
  
  # Define a range for the parameter being profiled
  # Use a multiplier for the range, ensuring it stays positive.
  # Adjust range dynamically based on the MLE value.
  range_multiplier <- 0.5 # +/- 50% of the MLE value
  
  # Calculate initial bounds
  lower_bound_range <- param_mle * (1 - range_multiplier)
  upper_bound_range <- param_mle * (1 + range_multiplier)
  
  # Ensure lower bound is not negative or too close to zero
  lower_bound_range <- max(1e-6, lower_bound_range) 
  
  # If the calculated range is very small, or if MLE is very small,
  # adjust to a more sensible absolute range to ensure enough points for plotting.
  if (abs(upper_bound_range - lower_bound_range) < 1e-3 * param_mle) { # If range is less than 0.1% of MLE
    # Try a fixed absolute increment
    abs_increment <- max(0.01, param_mle * 0.1) # At least 0.01 or 10% of MLE
    lower_bound_range <- max(1e-6, param_mle - abs_increment)
    upper_bound_range <- param_mle + abs_increment
  }
  
  # Generate a sequence of values for the parameter being profiled
  param_values <- seq(lower_bound_range, upper_bound_range, length.out = 100)
  
  # Calculate negative log-likelihood for each value in the sequence
  nll_values <- sapply(param_values, function(val) {
    temp_theta <- mle_params["theta"]
    temp_lambda <- mle_params["lambda"]
    temp_omega <- mle_params["omega"]
    
    # Assign the current profiled value to the correct parameter
    if (param_name == "theta") {
      temp_theta <- val
    } else if (param_name == "lambda") {
      temp_lambda <- val
    } else if (param_name == "omega") {
      temp_omega <- val
    }
    
    # Call the LL_ARSW function with the current set of parameters
    LL_ARSW(theta = temp_theta, lambda = temp_lambda, omega = temp_omega)
  })
  
  # Convert negative log-likelihood to likelihood
  likelihood_values <- exp(-nll_values)
  
  # Create a data frame for ggplot
  plot_data <- data.frame(
    Parameter = param_values,
    Likelihood = likelihood_values # Changed to Likelihood
  )
  
  # Get the maximum likelihood (at the MLE)
  # Recalculate it to ensure it's based on the exact MLEs from fit_ARSW
  max_likelihood <- exp(-LL_ARSW(theta = mle_params["theta"], lambda = mle_params["lambda"], omega = mle_params["omega"]))
  
  # Determine the x-axis label based on the parameter name
  x_label_expr <- switch(param_name,
                         "theta" = expression(theta),
                         "lambda" = expression(lambda),
                         "omega" = expression(omega),
                         param_name) # Fallback if name not matched
  
  windows() # Opens a new graphics device for each plot
  # Plot using ggplot2
  p <- ggplot(plot_data, aes(x = Parameter, y = Likelihood)) + # Changed y-axis aesthetic
    geom_line(color = "red", size = 1.5) + # Increased line size
    geom_point(aes(x = param_mle, y = max_likelihood), color = "black", size = 4, shape = 8) + # Increased point size
    labs(
      title = "", # Title is intentionally left blank as per your previous request
      x = x_label_expr, # Use the symbolic expression for x-axis
      y = "Likelihood (Data-II)" # Y-axis label
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Increased title size
          axis.title = element_text(size = 16), # Increased axis label size
          axis.text = element_text(size = 14)) # Increased axis tick text size
  
  print(p)
}

message("\nGenerating Likelihood Profile Plots...")

# Plot profile for theta
plot_likelihood_profile("theta", fit_ARSW, x)

# Plot profile for lambda
plot_likelihood_profile("lambda", fit_ARSW, x)

# Plot profile for omega
plot_likelihood_profile("omega", fit_ARSW, x)
