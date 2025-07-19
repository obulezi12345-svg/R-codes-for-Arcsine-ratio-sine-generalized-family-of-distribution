# # Set working directory (optional, adjust path as needed)
setwd("C:\\Users\\1012 G2\\Documents\\R files\\ARS-W") # Uncomment and adjust if needed

options(warn = -1)

library(AdequacyModel)
library(fitdistrplus)
#library(kolmim)
library(bbmle)
library(MASS)
#library(muhaz)
library(ggplot2)
library(psych)
library(pracma)
library(wesanderson)
#library(muhaz)


#Data-I
x1 = c(0.1, 0.74, 1, 1.08, 1.16, 1.3, 1.53, 1.71, 1.97, 2.3, 2.54, 3.47,
       0.33, 0.77, 1.02, 1.08, 1.2, 1.34, 1.59, 1.72, 2.02, 2.31, 2.54, 3.61,
       0.44, 0.92, 1.05, 1.09, 1.21, 1.36, 1.6, 1.76, 2.13, 2.4, 2.78, 4.02,
       0.56, 0.93, 1.07, 1.12, 1.22, 1.39, 1.63, 1.83, 2.15, 2.45, 2.93, 4.32,
       0.59, 0.96, 1.07, 1.13, 1.22, 1.44, 1.63, 1.95, 2.16, 2.51, 3.27, 4.58,
       0.72, 1, 1.08, 1.15, 1.24, 1.46, 1.68, 1.96, 2.22, 2.53, 3.42, 5.55)



#The third dataset is on the active repair times (hours) for an airborne communication transceiver,
#studied by Gadde and Durgamamba [28] and Ameeq et al. [29]
x =c(0.2,0.3, 0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.7, 0.8, 0.8, 1.0,1.0, 1.0, 1.0, 1.1, 1.3, 1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.2, 2.5,2.7, 3.0, 3.0, 3.3, 3.3, 4.0, 4.0, 4.5, 4.7, 5.0, 5.4, 5.4, 7.0,7.5, 8.8, 9.0, 10.3, 22.0, 24)

#40 observations representing time to failure (103h) of the turbocharger from one type of engine as given by Xuet al., [29]. These data were also utilized by Afify, Al-Mofleh, and Dey
x2 =c(7, 1.6, 7.7, 6.5, 2.0, 7.1, 8.3, 2.6, 8.1, 3.0, 6.7, 8.4,7.8, 3.5, 7.9, 8.4, 8.5, 3.9, 8.7, 4.5, 8.8, 4.6, 9.0, 6.0, 4.8, 8.0, 5.0, 7.3, 5.1, 7.3, 5.3,7.3, 5.4, 6.1, 5.6, 6.0, 5.8, 6.3, 7.0, 6.5)

data=x
length(x)
# Compute Q1, Q3, and IQR
Q1 <- quantile(data, 0.25)
Q3 <- quantile(data, 0.75)
IQR_value <- Q3 - Q1
length(x)
# Define outlier thresholds
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

# Identify outliers
outliers <- data[data < lower_bound | data > upper_bound]

# Print results
cat("Q1:", Q1, "\nQ3:", Q3, "\nIQR:", IQR_value, 
    "\nLower Bound:", lower_bound, 
    "\nUpper Bound:", upper_bound, 
    "\nOutliers:", outliers, "\n")


# Load library for more statistics
library(moments)  

# Compute key statistics
mean_x   <- mean(x)       # Mean
median_x <- median(x)     # Median
sd_x     <- sd(x)         # Standard Deviation
var_x    <- var(x)        # Variance
range_x  <- range(x)      # Min & Max
iqr_x    <- IQR(x)        # Interquartile Range (Q3 - Q1)
skew_x   <- skewness(x)   # Skewness
kurt_x   <- kurtosis(x)   # Kurtosis

# Print results
cat("Mean:", mean_x, "\n")
cat("Median:", median_x, "\n")
cat("Standard Deviation:", sd_x, "\n")
cat("Variance:", var_x, "\n")
cat("Range:", range_x[1], "to", range_x[2], "\n")
cat("Interquartile Range:", iqr_x, "\n")
cat("Skewness:", skew_x, "\n")
cat("Kurtosis:", kurt_x, "\n")





###################  ARSW distribution
dARSW <- function(par, x) {
  theta = par[1]
  lambda = par[2]
  omega = par[3]
  
  # Keep the check, but optimization should ideally prevent this
  if (theta <= 0 || lambda <= 0 || omega <= 0) {
    return(rep(0, length(x)))
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

# Cumulative Distribution Function (CDF)
# IMPORTANT CHANGE: Renamed 'q' to 'x' to match 'goodness.fit' expectation
pARSW = function(par, x) {
  theta = par[1]
  lambda = par[2]
  omega = par[3]
  if (theta <= 0 || lambda <= 0 || omega <= 0) {
    # Changed stop to return 0 or a very small number for numeric stability if called outside goodness.fit
    # For goodness.fit, it's best to let the underlying PDF/LL handle the constraints.
    # If this path is reached by goodness.fit, it means something went wrong with the params.
    warning("Parameters theta, lambda, and omega must be positive for pARSW.")
    return(rep(0, length(x)))
  }
  if (any(x < 0)) { # Changed q to x here too
    warning("x must be non-negative. Returning 0 for negative x values.")
  }
  
  # Handle x < 0 separately to ensure CDF is 0
  cdf_values <- numeric(length(x))
  valid_indices <- which(x >= 0)
  
  if (length(valid_indices) > 0) {
    G_x <- 1 - exp(-lambda * x[valid_indices]^omega) # Use x[valid_indices]
    sin_G_x <- sin(pi / 2 * G_x)
    
    numerator_frac <- theta * (1 - sin_G_x)
    denominator_frac <- theta + (2 - theta) * sin_G_x
    fraction_term <- numerator_frac / denominator_frac
    cdf_values[valid_indices] <- 1 - fraction_term
  }
  
  # Ensure values are within [0, 1]
  cdf_values[cdf_values < 0] <- 0
  cdf_values[cdf_values > 1] <- 1
  
  return(cdf_values)
}

LL_ARSW = function(theta, lambda, omega) {
  # Add a check here for parameters being valid, and return Inf if not
  if (theta <= 0 || lambda <= 0 || omega <= 0) {
    return(Inf)
  }
  
  # Ensure dARSW returns positive values for log()
  pdf_vals <- dARSW(c(theta, lambda, omega), x)
  
  # Handle cases where pdf_vals might be 0, leading to -Inf in log
  if (any(pdf_vals <= 0)) {
    return(Inf)
  }
  
  -sum(log(pdf_vals))
}

fit_ARSW = mle2(minuslogl = LL_ARSW,
                start = list(theta = 1, lambda = 1, omega = 0.5),
                data = list(x = x), # Ensure 'x' is passed correctly as data
                method = "L-BFGS-B", # Use a bounded method like L-BFGS-B
                lower = c(theta = 1e-6, lambda = 1e-6, omega = 1e-6)) # Set lower bounds

summary(fit_ARSW)

## gof for Data-I
gof.test_ARSW = goodness.fit(dARSW, pARSW, starts = coef(fit_ARSW), data = x,
                             method = "N", domain = c(0, 1000))
print(gof.test_ARSW) # Use print to display the results of goodness.fit




######### FW distribution
dFW = function(par, x){
  theta=par[1]
  lambda=par[2]
  (lambda / x^2 + theta) * exp(theta * x - lambda / x) * exp(-exp(theta * x - lambda / x))
}
pFW = function(par, x){
  theta=par[1]
  lambda=par[2]
  1 - exp(-exp(theta * x - lambda / x))
}
LL_FW  <- function(theta,lambda){-sum(log(dFW(c(theta,lambda),x)))}
fit_FW   <- mle2(minuslog=LL_FW,start=list(theta=0.01,lambda=1),
                  data=list(x),method = "N")
fit_FW
fit_FW@details$convergence
gof.test_FW   <- goodness.fit(dFW,pFW,starts=c(coef(fit_FW)),data=x,
                               method="N", domain = c(0,1000))

gof.test_FW




########## Gumbel distribution

## define pdf and cdf of Gumbel distribution
dGumbel <- function(parm, x){
  theta <- parm[1]
  lambda <- parm[2]
  exp(-(x - theta) / lambda) * exp(-exp(-(x - theta) / lambda)) / lambda
}

pGumbel <- function(parm, x){
  theta <- parm[1]
  lambda <- parm[2]
  exp(-exp(-(x - theta) / lambda))
}
LL_Gumbel  <- function(theta,lambda){-sum(log(dGumbel(c(theta,lambda),x)))}
fit_Gumbel   <- mle2(minuslog=LL_Gumbel,start=list(theta=1,lambda=1),
                     data=list(x),method = "N")
fit_Gumbel
fit_Gumbel@details$convergence
gof.test_Gumbel   <- goodness.fit(dGumbel,pGumbel,starts=c(coef(fit_Gumbel)),data=x,
                                  method="N", domain = c(0,1000))

gof.test_Gumbel



##define pdf and cdf of Gamma dist.
dGamma <- function(parm, x){
  theta <- parm[1]  # shape parameter
  lambda  <- parm[2]  # rate parameter
  (lambda^theta / gamma(theta)) * x^(theta - 1) * exp(-lambda * x)
}

pGamma <- function(parm, x){
  theta <- parm[1]  # shape parameter
  lambda  <- parm[2]  # rate parameter
  pgamma(x, shape = theta, rate = lambda)  # using built-in pgamma function
}
LL_Gamma  <- function(theta,lambda){-sum(log(dGamma(c(theta,lambda),x)))}
fit_Gamma   <- mle2(minuslog=LL_Gamma,start=list(theta=1,lambda=1),
                    data=list(x),method = "N")
fit_Gamma
fit_Gamma@details$convergence
gof.test_Gamma   <- goodness.fit(dGamma,pGamma,starts=c(coef(fit_Gamma)),data=x,
                                 method="N", domain = c(0,1000))

gof.test_Gamma 




#############
#log normal distribution
pdf_lnorm = function(par,x){
  theta=par[1]
  lambda = par[2]
  dlnorm(x,meanlog=theta, sdlog=lambda)
}
cdf_lnorm = function(par,x){
  theta=par[1]
  lambda = par[2]
  plnorm(x,meanlog=theta, sdlog=lambda)
}
LL_lnorm = function(theta,lambda){-sum(log(pdf_lnorm(c(theta,lambda),x)))}
fit_lnorm = mle2(minuslog=LL_lnorm,start=list(theta=0.91,lambda=2.5),data=list(x),method = "N")
fit_lnorm
fit_lnorm@details$convergence ## must equal zero
## gof
gof.test_lnorm=goodness.fit(pdf = pdf_lnorm, cdf = cdf_lnorm, starts = c(coef(fit_lnorm)), data=x,
                            method = "N", domain = c(0,1000))
gof.test_lnorm


##define pdf and cdf of Burr XII dist.
dBurr <- function(parm,x){
  theta <- parm[1]
  lambda  <- parm[2]
  theta*lambda*x^(theta-1)*(1+x^theta)^(-(lambda+1))
}
pBurr <- function(parm,x){
  theta <- parm[1]
  lambda  <- parm[2]
  1-(1+x^theta)^(-lambda)
}
LL_Burr  <- function(theta,lambda){-sum(log(dBurr(c(theta,lambda),x)))}
fit_Burr  <- mle2(minuslog=LL_Burr,start=list(theta=1, lambda=1),data=list(x),method = "N")
fit_Burr
fit_Burr@details$convergence  ## must equal zero

## gof
gof.test_Burr <- goodness.fit(dBurr,pBurr,starts=c(coef(fit_Burr)),data=x,
                              method="N", domain = c(0,1000))
gof.test_Burr









################# Weibull distribution
dWeibull <- function(parm, x) {
  theta <- parm[1]
  lambda <- parm[2]
  (lambda / theta) * (x / theta)^(lambda - 1) * exp(-(x / theta)^lambda)
}

pWeibull <- function(parm, x) {
  theta <- parm[1]
  lambda <- parm[2]
  1 - exp(-(x / theta)^lambda)
}
LL_Weibull  <- function(theta, lambda){-sum(log(dWeibull(c(theta, lambda),x)))}
fit_Weibull   <- mle2(minuslog=LL_Weibull,start=list(theta=1,lambda=1),
                      data=list(x),method = "N")
fit_Weibull
fit_Weibull@details$convergence
gof.test_Weibull   <- goodness.fit(dWeibull,pWeibull,starts=c(coef(fit_Weibull)),data=x,
                                   method="N", domain = c(0,1000))

gof.test_Weibull 


######## LCW
dLCW =function(par, x){
  theta=par[1]
  lambda=par[2]
  pi*theta*lambda*(x^(lambda-1))*(exp(-theta*x^lambda))*((1/(sin(pi*(1-exp(-theta*x^lambda)))))^2)*(exp(cos(pi*(1-exp(-theta*x^lambda)))/sin(pi*(1-exp(-theta*x^lambda)))))*(1+exp(cos(pi*(1-exp(-theta*x^lambda)))/sin(pi*(1-exp(-theta*x^lambda)))))^(-2)
}

pLCW =function(par, x){
  theta=par[1]
  lambda=par[2]
  (1+exp(cos(pi*(1-exp(-theta*x^lambda)))/sin(pi*(1-exp(-theta*x^lambda)))))^(-1)
}
LL_LCW  <- function(theta,lambda){-sum(log(dLCW(c(theta,lambda),x)))}
fit_LCW   <- mle2(minuslog=LL_LCW,start=list(theta=1,lambda=0.1),
                   data=list(x),method = "N")
fit_LCW
fit_LCW@details$convergence
gof.test_LCW   <- goodness.fit(dLCW,pLCW,starts=c(coef(fit_LCW)),data=x,
                                method="N", domain = c(0,1000))

gof.test_LCW




# --- Plotting Code for Data-I ---
# Load necessary package (already loaded above)
# library(ggplot2)
windows()
# Convert your data into a data frame
df <- data.frame(value = x) # Corrected: used 'x' instead of 'data'

# Create the violin plot with a superimposed boxplot
ggplot(df, aes(x = "", y = value)) +  # x is an empty string for a single-group plot
  geom_violin(trim = FALSE, fill = "gray", alpha = 1.5, color = "darkgreen", linewidth = 3) +  # Violin plot with red outline
  geom_boxplot(width = 0.2, fill = "skyblue", outlier.shape = NA, color = "darkgreen", linewidth = 1.5) +  # Boxplot with red lines
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.y = element_blank(),  # Remove y-axis text labels
    axis.ticks.y = element_blank(),# Remove y-axis ticks
    axis.title.x = element_text(size = 25)
  ) +
  labs(title = "",
       x = "Data-III",
       y = NULL)  # Remove y-axis label



# Load necessary package
# Load necessary package
library(stats)

data = x
data = sort(data)
y    = seq(0, ceiling(max(x)), length.out = 500)
Cdf_ecdf = ecdf(x = data)
cdf_emp  = Cdf_ecdf(data)
Sue_ecdf = 1 - cdf_emp
p = (rank(data)) / (length(data) + 1)

# Open a window for plotting
windows()
par(mfrow = c(1, 1))

# Plot histogram without vertical axis line, label, and legend
h <- hist(data, probability = TRUE, col = "darkgray", alpha = 1.5, 
          ylab = "", ylim = c(0, 0.3), xlab = "Data-III",  
          main = "", cex.main = 2, cex.axis = 2, cex.lab = 2,
          axes = FALSE, border = NA)  # Remove default borders

# Add thick borders to histogram bars
rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density, 
     border = "black", lwd = 3)  # Adjusts the bar width

# Add only the x-axis (horizontal axis)
axis(1)  # Add x-axis with labels and ticks

# Keep the horizontal axis (y-axis) without labels or vertical line
axis(2, labels = FALSE, tick = FALSE)  # Removes the vertical axis line and labels

# Superimpose the WALD density function
lines(y, dARSW(gof.test_ARSW$mle, y), col = "red", lwd = 5)


# Removed legend





windows()
par(mfrow = c(1, 1))

## Graph of Empirical CDF
plot(ecdf(data), verticals = TRUE, ylab = "", xlab = "Data-III",  
     main = "", cex.main = 2, cex.axis = 1.5, cex.lab = 3)

# Superimpose the Wald CDF
lines(y, pARSW(gof.test_ARSW$mle, y), col = 2, lwd = 3)

# Legend removed





windows()
par(mfrow = c(1, 1))

## Graph of Empirical Survival Function
plot(stepfun(data, c(1, Sue_ecdf)), verticals = TRUE, 
     xlab = "Data-III", ylab = "", main = "", 
     cex.main = 2, cex.axis = 1.5, cex.lab = 3, lwd = 3)

# Superimpose the Wald Survival Function
lines(y, 1 - pARSW(gof.test_ARSW$mle, y), col = "red", lwd = 5)

# Legend removed


windows()
TTT(x, col = 2, lwd = 5, grid = TRUE, lty = 2) 

# Remove x-axis labels
axis(1, labels = FALSE)

# Add "Data-I" as the x-axis label
mtext("Data-III", side = 1, line = 2.5, cex = 2)

# Add the main title
title(main = "", cex.main = 2, cex.lab = 1.5, cex.axis = 1.5)




windows()
plot(pARSW(gof.test_ARSW$mle,data),p,type="l",pch=1,xlim=c(0,1),ylim=c(0,1),col="red",xlab="Data-III",
     ylab="",lwd=3, main="", cex.main=2, cex.lab=3, cex.axis=1.5)
segments(0,0,1,1, lwd=3)


























