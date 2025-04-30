# ===================================================
# LOAD LIBRARIES
# ===================================================
library(CASdatasets)
library(dplyr)
library(rjags)
library(coda)

# ===================================================
# RYTGAARD1990 EXAMPLE DATA
# ===================================================
rytgaard1990_input <- list(
  y = c(2.495, 2.120, 2.095, 1.700, 1.650, 1.985, 1.810, 1.625, 
        3.215, 2.105, 1.765, 1.715, 19.180, 1.915, 1.790, 1.755),
  n = c(5, 3, 4, 0, 4)
)

# Add metadata for JAGS model
rytgaard1990_input$N_y <- length(rytgaard1990_input$y)
rytgaard1990_input$N_n <- length(rytgaard1990_input$n)
rytgaard1990_input$min_y <- min(rytgaard1990_input$y)

# ===================================================
# LOAD AND PREPARE CASDATASETS DATA
# ===================================================
data(itamtplcost)

# Extract year and rename 'UltimateCost' for clarity
itamtplcost$year <- format(as.Date(itamtplcost$Date, format = "%d/%m/%Y"), "%Y")

itamtplcost <- itamtplcost %>%
  rename(claim_amount = UltimateCost) %>%
  select(year, claim_amount)

itamtplcost_input <- list(
  y = itamtplcost$claim_amount,
  n = itamtplcost %>% group_by(year) %>% summarise(count=n()) %>% pull(count)
)

itamtplcost_input$N_y <- length(itamtplcost_input$y)
itamtplcost_input$N_n <- length(itamtplcost_input$n)
itamtplcost_input$min_y <- min(itamtplcost_input$y)

# ===================================================
# DEFINE BAYESIAN MODEL (POISSON-PARETO)
# ===================================================
model_code <- "
model {
  for(i in 1:N_y) {
    y[i] ~ dpar(alpha, beta)
  }
  
  for(i in 1:N_n) {
    n[i] ~ dpois(theta)
  }
  
  alpha ~ dgamma(1, 0.0001)
  beta ~ dgamma(1, 0.0001)I(, min_y)
  theta ~ dgamma(1, 0.0001)
}
"

# Initial values for three MCMC chains
inits_list <- list(
  list(alpha = 0.00001, beta = 0.00001, theta = 0.00001, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 123),
  list(alpha = 100000, beta = 1, theta = 100000, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 456),
  list(alpha = 3.076, beta = 1.625, theta = 3.2, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 789)
)

# ===================================================
# RUN JAGS MODEL
# ===================================================
set.seed(123)

jags_model <- jags.model(
  textConnection(model_code),
  data = rytgaard1990_input,
  inits = inits_list,
  n.chains = 3,
  n.adapt = 20000
  )

# Sample from posterior
samples <- coda.samples(
  jags_model,
  variable.names = c("alpha", "beta", "theta"),
  n.iter = 30000
  )

# Posterior summary & diagnostics
summary(samples)
plot(samples)

# Density for each parameter
plot(density(as.matrix(samples[, "alpha"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples[, "beta"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples[, "theta"])), main="", xlab="", ylab="")

# Compute the average expected claim amount
pareto_ev <- function(a, b) ifelse(a > 1, (a * b) / (a - 1), NA)

posterior <- as.matrix(samples)
expected_values <- pareto_ev(posterior[, "alpha"], posterior[, "beta"])

mean(expected_values, na.rm=TRUE)
sd(expected_values, na.rm=TRUE)
quantile(expected_values, c(0.025, 0.5, 0.975), na.rm=TRUE)

plot(density(expected_values, na.rm=TRUE), main="", xlab="", ylab="")

# ===================================================
# PLOT 1: EMPIRICAL VS POSTERIOR PARETO CDF
# ===================================================
y <- rytgaard1990_input$y

alpha <- mean(posterior[, "alpha"])
beta <- mean(posterior[, "beta"])

# Define Pareto CDF
pareto_cdf <- function(x, alpha, beta) {
  ifelse(x < beta, 0, 1 - (beta / x)^alpha)
}

# Plot empirical CDF
plot(
  ecdf(y),
  col = "blue",
  main = "",
  xlab = "",
  ylab = "",
  lwd = 2
)

# Overlay posterior predictive Pareto CDF
x_vals <- seq(min(y), max(y), length.out = 500)
lines(x_vals, pareto_cdf(x_vals, alpha, beta),
      col = "red", lwd = 2, lty = 2)

legend("bottomright",
       legend = c("Empirical CDF", "Pareto(3.076,1.591) CDF"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# ===================================================
# PLOT 2: EMPIRICAL VS POSTERIOR POISSON CDF
# ===================================================
n <- rytgaard1990_input$n
theta <- mean(posterior[, "theta"])

# Plot empirical CDF for Poisson counts
plot(
  ecdf(n),
  col = "blue",
  main = "",
  xlab = "",
  ylab = "",
  lwd = 2,
  xlim = c(min(n), max(n) + 3)  # extend x-axis for visibility
)

# Overlay theoretical Poisson CDF
x_vals <- min(n):(max(n) + 3)
lines(x_vals, ppois(x_vals, lambda = theta),
      col = "red", lwd = 2, lty = 2, type = "s")  # 's' for step

legend("bottomright",
       legend = c("Empirical CDF", "Poisson(3.399) CDF"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# ===================================================
# DIAGNOSTICS: Trace Plots, Gelman-Rubin, and Autocorrelation
# ===================================================

# Trace plots for each parameter
traceplot(samples[, c("alpha", "beta", "theta")], main="", xlab="")

# Gelman-Rubin diagnostic
gelman.diag(samples, autoburnin = FALSE)

# Gelman-Rubin plots
params <- c("alpha", "beta", "theta")
for (param in params) {
  gelman.plot(samples[, param], xlab = "", ylab = "")
}

# Autocorrelation diagnostics
autocorr.plot(samples[, "alpha"], ask=FALSE)
autocorr.plot(samples[, "beta"], ask=FALSE)
autocorr.plot(samples[, "theta"], ask=FALSE)

autocorr.diag(samples[, "alpha"], lags = 1:10)
autocorr.diag(samples[, "beta"], lags = 1:10)
autocorr.diag(samples[, "theta"], lags = 1:10)

# Rerunning chains with updated settings
jags_model_rerun <- jags.model(
  textConnection(model_code),
  data = rytgaard1990_input,
  inits = inits_list,
  n.chains = 3,
  n.adapt = 20000
)

# Sample from posterior
samples_rerun <- coda.samples(
  jags_model_rerun,
  variable.names = c("alpha", "beta", "theta"),
  n.iter = 300000,
  thin = 10
)

# Trace plots for rerun samples
traceplot(samples_rerun[, c("alpha", "beta", "theta")], main="", xlab="")

# Gelman-Rubin diagnostic for rerun samples
gelman.diag(samples_rerun, autoburnin = FALSE)

# Autocorrelation diagnostics for rerun samples
autocorr.plot(samples_rerun[, "alpha"], ask=FALSE)
autocorr.plot(samples_rerun[, "beta"], ask=FALSE)
autocorr.plot(samples_rerun[, "theta"], ask=FALSE)

autocorr.diag(samples_rerun[, "alpha"], lags = 1:10)
autocorr.diag(samples_rerun[, "beta"], lags = 1:10)
autocorr.diag(samples_rerun[, "theta"], lags = 1:10)

# ===================================================
# PREDICTIONS: Poisson and Pareto Distribution
# ===================================================

# Poisson predictions (n values)
poisson_predictions <- numeric(15)
for (n_value in 0:14) {
  poisson_probabilities <- (posterior[, "theta"]^n_value) * exp(-posterior[, "theta"]) / factorial(n_value)
  poisson_predictions[n_value + 1] <- mean(poisson_probabilities)
}
poisson_predictions

# Pareto predictions (y values)
pareto_predictions <- numeric(1000)
uniform_samples <- runif(1000)

for (i in seq_along(uniform_samples)) {
  U <- uniform_samples[i]
  pareto_samples <- posterior[, "beta"] / (1 - U)^(1 / posterior[, "alpha"])
  pareto_predictions[i] <- mean(pareto_samples)
}

# Plot Pareto predictive distribution
plot(density(pareto_predictions, na.rm=TRUE), main="", xlab="", ylab="")

# ===================================================
# PREDICTIONS: Poisson-Pareto Distribution
# ===================================================

set.seed(123)

# Function to simulate Poisson samples
pois <- function(lambda, U) {
  p <- exp(-lambda)
  F <- p
  n <- 0
  
  while (U > F) {
    n <- n + 1
    p <- p * lambda / n
    F <- F + p
  }
  return(n)
}

# Simulation parameters
n_sim <- 1000
theta_values <- posterior[, "theta"]
alpha_values <- posterior[, "alpha"]
beta_values <- posterior[, "beta"]

# Simulate Poisson samples
uniform_poisson <- runif(n_sim)
poisson_results <- numeric(n_sim)

for (i in seq_along(uniform_poisson)) {
  U <- uniform_poisson[i]
  poisson_samples <- numeric(length(theta_values))
  
  for (j in seq_along(theta_values)) {
    poisson_samples[j] <- pois(theta_values[j], U)
  }
  
  poisson_results[i] <- round(mean(poisson_samples))
}

# Simulate aggregate claims
final_results <- numeric(length(poisson_results))

for (i in seq_along(poisson_results)) {
  n <- poisson_results[i]
  pareto_samples <- numeric(n)
  
  uniform_pareto <- runif(n)
  
  for (j in seq_along(uniform_pareto)) {
    U <- uniform_pareto[j]
    pareto_samples[j] <- mean(beta_values / (1 - U)^(1 / alpha_values))
  }
  
  final_results[i] <- sum(pareto_samples)
}

# Fit candidate distributions via moment matching
mean_S <- mean(final_results)
var_S <- var(final_results)

gamma_shape <- mean_S^2 / var_S
gamma_rate  <- mean_S / var_S

sigma2_ln <- log(1 + var_S / mean_S^2)
mu_ln     <- log(mean_S) - sigma2_ln / 2

normal_mean <- mean_S
normal_sd   <- sqrt(var_S)

# Define support for density plots
x_vals <- seq(min(final_results), max(final_results), length.out = 1000)

gamma_densities   <- dgamma(x_vals, shape = gamma_shape, rate = gamma_rate)
lognorm_densities <- dlnorm(x_vals, meanlog = mu_ln, sdlog = sqrt(sigma2_ln))
normal_densities  <- dnorm(x_vals, mean = normal_mean, sd = normal_sd)

hist_density <- hist(final_results, plot = FALSE, probability = TRUE)$density
max_y <- max(gamma_densities, lognorm_densities, normal_densities, hist_density)

# Plotting
hist(final_results, probability = TRUE, breaks = 100,
     col = "gray90", border = "white", main = "",
     xlab = "", ylab = "", ylim = c(0, max_y * 1.05))
curve(dgamma(x, shape = gamma_shape, rate = gamma_rate), col = "blue", lwd = 2, add = TRUE)
curve(dlnorm(x, meanlog = mu_ln, sdlog = sqrt(sigma2_ln)), col = "green", lwd = 2, add = TRUE)
curve(dnorm(x, mean = normal_mean, sd = normal_sd), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Gamma", "Lognormal", "Normal"),
       col = c("blue", "green", "red"), lwd = 2)

# Summary statistics
median(final_results)
quantile(final_results, probs = c(0.90, 0.95, 0.99))
max(final_results)

## New data


