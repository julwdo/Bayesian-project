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

# Express claim_amount in millions
itamtplcost$claim_amount <- round(itamtplcost$claim_amount / 1000, 3)

itamtplcost_input <- list(
  y = itamtplcost$claim_amount,
  n = as.numeric(itamtplcost %>% group_by(year) %>% summarise(count=n()) %>% pull(count))
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

# ---------------------------------------------------
# RUN JAGS MODEL
# ---------------------------------------------------
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

# ---------------------------------------------------
# PLOT 1: EMPIRICAL VS POSTERIOR PARETO CDF
# ---------------------------------------------------
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

# ---------------------------------------------------
# PLOT 2: EMPIRICAL VS POSTERIOR POISSON CDF
# ---------------------------------------------------
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

# ---------------------------------------------------
# DIAGNOSTICS: Trace Plots, Gelman-Rubin, and Autocorrelation
# ---------------------------------------------------

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

# ---------------------------------------------------
# PREDICTIONS: Poisson and Pareto Distribution
# ---------------------------------------------------
posterior_rerun <- as.matrix(samples_rerun)

# Poisson predictions (n values)
poisson_predictions <- numeric(15)
for (n_value in 0:14) {
  poisson_probabilities <- (posterior_rerun[, "theta"]^n_value) * exp(-posterior_rerun[, "theta"]) / factorial(n_value)
  poisson_predictions[n_value + 1] <- mean(poisson_probabilities)
}
poisson_predictions

# Pareto predictions (y values)
pareto_predictions <- numeric(1000)
uniform_samples <- runif(1000)

for (i in seq_along(uniform_samples)) {
  U <- uniform_samples[i]
  pareto_samples <- posterior_rerun[, "beta"] / (1 - U)^(1 / posterior_rerun[, "alpha"])
  pareto_predictions[i] <- mean(pareto_samples)
}

# Plot Pareto predictive distribution
plot(density(pareto_predictions, na.rm=TRUE), main="", xlab="", ylab="")

# ---------------------------------------------------
# PREDICTIONS: Poisson-Pareto Distribution
# ---------------------------------------------------
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
theta_values <- posterior_rerun[, "theta"]
alpha_values <- posterior_rerun[, "alpha"]
beta_values <- posterior_rerun[, "beta"]

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
gamma_rate <- mean_S / var_S

sigma2_ln <- log(1 + var_S / mean_S^2)
mu_ln <- log(mean_S) - sigma2_ln / 2

normal_mean <- mean_S
normal_sd <- sqrt(var_S)

# Define support for density plots
x_vals <- seq(min(final_results), max(final_results), length.out = 1000)

gamma_densities <- dgamma(x_vals, shape = gamma_shape, rate = gamma_rate)
lognorm_densities <- dlnorm(x_vals, meanlog = mu_ln, sdlog = sqrt(sigma2_ln))
normal_densities <- dnorm(x_vals, mean = normal_mean, sd = normal_sd)

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

# ===================================================
# NEW DATA: COMPARISON OF STATISTICS
# ===================================================

summary(rytgaard1990_input$y)
summary(itamtplcost_input$y)

summary(rytgaard1990_input$n)
summary(itamtplcost_input$n)

max(rytgaard1990_input$y)-min(rytgaard1990_input$y)
max(itamtplcost_input$y)-min(itamtplcost_input$y)

# ===================================================
# NEW DATA: POISSON-WEIBULL MODEL
# ===================================================

model_code_weib <- "
model {
  for(i in 1:N_y) {
    y[i] ~ dweib(alpha, beta)
  }
  
  for(i in 1:N_n) {
    n[i] ~ dpois(theta)
  }
  
  alpha ~ dgamma(1, 0.1)
  beta ~ dgamma(1, 0.1)
  theta ~ dgamma(1, 0.0001)
}
"

# Initial values for three MCMC chains
inits_list_weib <- list(
  list(alpha = 1.6, beta = 1.5, theta = 0.1, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 123),
  list(alpha = 2, beta = 2, theta = 10, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 456),
  list(alpha = 1.595, beta = 1.133, theta = 28.563, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 789)
)

# ---------------------------------------------------
# RUN JAGS MODEL
# ---------------------------------------------------
set.seed(123)

jags_model_weib <- jags.model(
  textConnection(model_code_weib),
  data = itamtplcost_input,
  inits = inits_list_weib,
  n.chains = 3,
  n.adapt = 10000
)

samples_weib <- coda.samples(
  jags_model_weib,
  variable.names = c("alpha", "beta", "theta"),
  n.iter = 300000,
  thin = 10
)

summary(samples_weib)

gelman.diag(samples_weib, autoburnin = FALSE)

params <- c("alpha", "beta", "theta")
for (param in params) {
  gelman.plot(samples_weib[, param], xlab = "", ylab = "")
}

autocorr.plot(samples_weib[, "alpha"], ask=FALSE)
autocorr.plot(samples_weib[, "beta"], ask=FALSE)
autocorr.plot(samples_weib[, "theta"], ask=FALSE)

autocorr.diag(samples_weib[, "alpha"], lags = 1:10)
autocorr.diag(samples_weib[, "beta"], lags = 1:10)
autocorr.diag(samples_weib[, "theta"], lags = 1:10)

traceplot(samples_weib[, c("alpha", "beta", "theta")], main="", xlab="")

plot(density(as.matrix(samples_weib[, "alpha"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples_weib[, "beta"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples_weib[, "theta"])), main="", xlab="", ylab="")

posterior_weib <- as.matrix(samples_weib)

# ---------------------------------------------------
# PREDICTIONS: Poisson-Weibull Distribution
# ---------------------------------------------------

# Simulation parameters
n_sim <- 1000
alpha_values_weib <- posterior_weib[, "alpha"]
beta_values_weib <- posterior_weib[, "beta"]
theta_values_weib <- posterior_weib[, "theta"]

message("Starting Poisson simulation...")

# Simulate Poisson samples
uniform_poisson <- runif(n_sim)
poisson_results <- numeric(n_sim)

for (i in seq_along(uniform_poisson)) {
  if (i %% 100 == 0) message("Simulating Poisson sample ", i, "/", n_sim)
  
  U <- uniform_poisson[i]
  poisson_samples <- numeric(length(theta_values_weib))
  
  for (j in seq_along(theta_values_weib)) {
    poisson_samples[j] <- pois(theta_values_weib[j], U)
  }
  
  poisson_results[i] <- round(mean(poisson_samples))
}

message("Poisson simulation completed.")
message("Starting aggregate claim simulation...")

# Simulate aggregate claims
aggregate_claims_weib <- numeric(length(poisson_results))

for (i in seq_along(poisson_results)) {
  if (i %% 100 == 0) message("Simulating aggregate claim ", i, "/", length(poisson_results))
  
  n <- poisson_results[i]
  weib_samples <- numeric(n)
  uniform_weib <- runif(n, min = 1e-10, max = 1 - 1e-10)  # avoid 0 and 1 exactly
  
  for (j in seq_along(uniform_weib)) {
    U <- uniform_weib[j]
    
    weib_sample <- 1000*beta_values_weib * (-log(1 - U))^(1 / alpha_values_weib)
    weib_samples[j] <- mean(weib_sample)
  }
  
  aggregate_claims_weib[i] <- sum(weib_samples)
}

message("Aggregate simulation completed.")

# Summary statistics
median(aggregate_claims_weib)
quantile(aggregate_claims_weib, probs = c(0.90, 0.95, 0.99))
max(aggregate_claims_weib)

# Sample statistics
S_observed <- itamtplcost %>% group_by(year) %>% summarise(sum(claim_amount)) %>% pull('sum(claim_amount)')
plot(S_observed)

# ---------------------------------------------------
# Fit & Plot Distributions to Aggregate Claims
# ---------------------------------------------------

# Estimate distribution parameters via moment matching
mean_S <- mean(aggregate_claims_weib)
var_S <- var(aggregate_claims_weib)

gamma_shape <- mean_S^2 / var_S
gamma_rate <- mean_S / var_S

sigma2_ln <- log(1 + var_S / mean_S^2)
mu_ln <- log(mean_S) - sigma2_ln / 2

normal_mean <- mean_S
normal_sd <- sqrt(var_S)

# Define support for density plots
x_vals <- seq(min(aggregate_claims_weib), max(aggregate_claims_weib), length.out = 1000)

gamma_densities <- dgamma(x_vals, shape = gamma_shape, rate = gamma_rate)
lognorm_densities <- dlnorm(x_vals, meanlog = mu_ln, sdlog = sqrt(sigma2_ln))
normal_densities <- dnorm(x_vals, mean = normal_mean, sd = normal_sd)

hist_density <- hist(aggregate_claims_weib, plot = FALSE, probability = TRUE)$density
max_y <- max(gamma_densities, lognorm_densities, normal_densities, hist_density)

# Plotting
hist(aggregate_claims_weib, probability = TRUE, breaks = 100,
     col = "gray90", border = "white", main = "",
     xlab = "", ylab = "", ylim = c(0, max_y * 1.05))
curve(dgamma(x, shape = gamma_shape, rate = gamma_rate), col = "blue", lwd = 2, add = TRUE)
curve(dlnorm(x, meanlog = mu_ln, sdlog = sqrt(sigma2_ln)), col = "green", lwd = 2, add = TRUE)
curve(dnorm(x, mean = normal_mean, sd = normal_sd), col = "red", lwd = 2, add = TRUE)
legend("right", legend = c("Gamma", "Lognormal", "Normal"),
       col = c("blue", "green", "red"), lwd = 2)

# ===================================================
# NEW DATA: POISSON-LOGNORMAL MODEL
# ===================================================

model_code_pois_lognorm <- "
model {
  for(i in 1:N_y) {
    y[i] ~ dlnorm(mu_pois_lognorm, tau_pois_lognorm)
  }
  
  for(i in 1:N_n) {
    n[i] ~ dpois(theta_pois_lognorm)
  }
  
  mu_pois_lognorm ~ dnorm(0, 1.0E-8)
  sigma_pois_lognorm ~ dgamma(1, 0.0001)
  tau_pois_lognorm <- 1 / pow(sigma_pois_lognorm, 2)
  theta_pois_lognorm ~ dgamma(1, 0.0001)
}
"

# MLE estimates for initialization
log_y <- log(itamtplcost_input$y)
mu_hat_pois_lognorm <- mean(log_y)
sigma_hat_pois_lognorm <- sd(log_y)

tau_hat_pois_lognorm <- 1/sigma_hat_pois_lognorm^2

theta_hat_pois_lognorm <- mean(itamtplcost_input$n)

# Initial values for the model
inits_pois_lognorm <- list(
  list(mu_pois_lognorm = 0.00001,
       sigma_pois_lognorm = 0.00001,
       theta_pois_lognorm = 0.00001,
       .RNG.name = "base::Wichmann-Hill",
       .RNG.seed = 123),
  list(mu_pois_lognorm = 100000,
       sigma_pois_lognorm = 100000,
       theta_pois_lognorm = 100000,
       .RNG.name = "base::Wichmann-Hill",
       .RNG.seed = 456),
  list(mu_pois_lognorm = mu_hat_pois_lognorm,
       sigma_pois_lognorm = sigma_hat_pois_lognorm,
       theta_pois_lognorm = theta_hat_pois_lognorm,
       .RNG.name = "base::Wichmann-Hill",
       .RNG.seed = 789)
)

hist(itamtplcost_input$y)

# ---------------------------------------------------
# RUN JAGS MODEL
# ---------------------------------------------------
jags_model_pois_lognorm <- jags.model(
  textConnection(model_code_pois_lognorm),
  data = itamtplcost_input,
  inits = inits_pois_lognorm,
  n.chains = 3,
  n.adapt = 10000
)

samples_pois_lognorm <- coda.samples(
  jags_model_pois_lognorm,
  variable.names = c("mu_pois_lognorm",
                     "tau_pois_lognorm",
                     "theta_pois_lognorm"),
  n.iter = 300000,
  thin = 10
)

# Diagnostics & summaries
params <- c("mu_pois_lognorm", "tau_pois_lognorm", "theta_pois_lognorm")

# Gelman-Rubin diagnostics
for (param in params) {
  gelman.plot(samples_pois_lognorm[, param], xlab = "", ylab = "")
}
gelman.diag(samples_pois_lognorm, autoburnin = FALSE)

# Autocorrelation diagnostics
autocorr.plot(samples_pois_lognorm[, "mu_pois_lognorm"], ask=FALSE)
autocorr.plot(samples_pois_lognorm[, "tau_pois_lognorm"], ask=FALSE)
autocorr.plot(samples_pois_lognorm[, "theta_pois_lognorm"], ask=FALSE)

autocorr.diag(samples_pois_lognorm[, "mu_pois_lognorm"], lags = 1:10)
autocorr.diag(samples_pois_lognorm[, "tau_pois_lognorm"], lags = 1:10)
autocorr.diag(samples_pois_lognorm[, "theta_pois_lognorm"], lags = 1:10)

# Traceplots
traceplot(samples_pois_lognorm[, params], main="", xlab="")

# Posterior summary
summary(samples_pois_lognorm)

mean(rytgaard1990_input$y)

# Density plots
plot(density(as.matrix(samples_pois_lognorm[, "mu_pois_lognorm"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples_pois_lognorm[, "tau_pois_lognorm"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples_pois_lognorm[, "theta_pois_lognorm"])), main="", xlab="", ylab="")

# Lognormal expected value calculation
lognormal_ev <- function(mu, tau) {
  sigma2 <- 1 / tau
  exp(mu + sigma2 / 2)
}

posterior_pois_lognorm <- as.matrix(samples_pois_lognorm)

expected_y_pois_lognorm <- lognormal_ev(posterior_pois_lognorm[, "mu_pois_lognorm"],
                                        posterior_pois_lognorm[, "tau_pois_lognorm"])
mean(expected_y_pois_lognorm)
sd(expected_y_pois_lognorm)
quantile(expected_y_pois_lognorm, c(0.025, 0.5, 0.975), na.rm = TRUE)

plot(density(expected_y_pois_lognorm, na.rm=TRUE), main="", xlab="", ylab="")

# ---------------------------------------------------
# PLOT 1: Empirical vs Posterior Lognormal CDF
# ---------------------------------------------------
y_obs <- itamtplcost_input$y
mu_pois_lognorm_post_mean <- mean(posterior_pois_lognorm[, "mu_pois_lognorm"])
sigma_pois_lognorm_post_mean <- mean(1 / sqrt(posterior_pois_lognorm[, "tau_pois_lognorm"]))

tau_pois_lognorm_post_mean <- mean(posterior_pois_lognorm[, "tau_pois_lognorm"])

lognormal_cdf <- function(x, mu, sigma) {
  pnorm(log(x), mean = mu, sd = sigma)
}

plot(
  ecdf(y_obs),
  col = "blue",
  main = "",
  xlab = "",
  ylab = "",
  lwd = 2
)

x_vals <- seq(min(y_obs), max(y_obs), length.out = 500)
lines(x_vals,
      lognormal_cdf(x_vals, mu_pois_lognorm_post_mean, sigma_pois_lognorm_post_mean),
      col = "red", lwd = 2, lty = 2)

legend("bottomright",
       legend = c("Empirical CDF",
                  sprintf("Lognormal(%.3f, %.3f) CDF",
                          mu_pois_lognorm_post_mean,
                          tau_pois_lognorm_post_mean)),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# ---------------------------------------------------
# PLOT 2: EMPIRICAL VS POSTERIOR POISSON CDF
# ---------------------------------------------------
n <- itamtplcost_input$n
theta_pois_lognorm_post_mean <- mean(posterior_pois_lognorm[, "theta_pois_lognorm"])

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
lines(x_vals, ppois(x_vals, lambda = theta_pois_lognorm_post_mean),
      col = "red", lwd = 2, lty = 2, type = "s")  # 's' for step

legend("bottomright",
       legend = c("Empirical CDF", "Poisson(28.628) CDF"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# ---------------------------------------------------
# Posterior Predictive Sampling for y
# ---------------------------------------------------
lognorm_predictions <- numeric(1000)
uniform_samples <- runif(1000)

mu_pois_lognorm_posterior <- posterior_pois_lognorm[, "mu_pois_lognorm"]
sigma_pois_lognorm_posterior <- 1 / sqrt(posterior_pois_lognorm[, "tau_pois_lognorm"])

for (i in seq_along(uniform_samples)) {
  U <- uniform_samples[i]
  
  normal_samples <- qnorm(U, mean = mu_pois_lognorm_posterior, sd = sigma_pois_lognorm_posterior)
  lognorm_samples <- exp(normal_samples)
  
  lognorm_predictions[i] <- mean(lognorm_samples)
}

plot(density(lognorm_predictions, na.rm = TRUE), main = "", xlab = "", ylab = "")

# ---------------------------------------------------
# PREDICTIONS: Poisson-Lognormal Distribution
# ---------------------------------------------------
set.seed(123)

# Simulation parameters
n_sim <- 1000

theta_pois_lognorm_posterior <- posterior_pois_lognorm[, "theta_pois_lognorm"]

# Simulate number of claims
uniform_poisson <- runif(n_sim)
poisson_results <- numeric(n_sim)

for (i in seq_along(uniform_poisson)) {
  U <- uniform_poisson[i]
  
  cat("Simulation", i, "\n")
  
  poisson_samples <- numeric(length(theta_pois_lognorm_posterior))
  
  for (j in seq_along(theta_pois_lognorm_posterior)) {
    poisson_samples[j] <- pois(theta_pois_lognorm_posterior[j], U)
  }
  
  poisson_results[i] <- round(mean(poisson_samples))
}

# Compare empirical vs simulated claim counts
summary(itamtplcost_input$n) # Empirical
summary(poisson_results) # Fitted

# Simulate aggregate claims
aggregate_claims_pois_lognorm <- numeric(length(poisson_results))

for (i in seq_along(poisson_results)) {
  n <- poisson_results[i]
  
  lognorm_samples <- numeric(n)
  uniform_lognorm <- runif(n)
  
  for (j in seq_along(uniform_lognorm)) {
    U <- uniform_lognorm[j]
    
    normal_samples <- qnorm(U, mean = mu_pois_lognorm_posterior, sd = sigma_pois_lognorm_posterior)
    lognorm_samples[j] <- mean(exp(normal_samples))
  }
  
  aggregate_claims_pois_lognorm[i] <- sum(lognorm_samples)
}

# Compare to empirical aggregate claim distribution
help <- itamtplcost %>%
  group_by(year) %>%
  summarise(sum = sum(claim_amount), .groups = "drop")

summary(help$sum) # Empirical
summary(aggregate_claims_pois_lognorm) # Fitted

# ---------------------------------------------------
# Fit & Plot Distributions to Aggregate Claims
# ---------------------------------------------------

# Estimate distribution parameters via moment matching
mean_S <- mean(aggregate_claims_pois_lognorm)
var_S <- var(aggregate_claims_pois_lognorm)

gamma_shape <- mean_S^2 / var_S
gamma_rate <- mean_S / var_S

sigma2_ln <- log(1 + var_S / mean_S^2)
mu_ln <- log(mean_S) - sigma2_ln / 2

normal_mean <- mean_S
normal_sd <- sqrt(var_S)

# Define support for density plots
x_vals <- seq(min(aggregate_claims_pois_lognorm), max(aggregate_claims_pois_lognorm), length.out = 1000)

gamma_densities <- dgamma(x_vals, shape = gamma_shape, rate = gamma_rate)
lognorm_densities <- dlnorm(x_vals, meanlog = mu_ln, sdlog = sqrt(sigma2_ln))
normal_densities <- dnorm(x_vals, mean = normal_mean, sd = normal_sd)

hist_density <- hist(aggregate_claims_pois_lognorm, plot = FALSE, probability = TRUE)$density
max_y <- max(gamma_densities, lognorm_densities, normal_densities, hist_density)

# Plotting
hist(aggregate_claims_pois_lognorm, probability = TRUE, breaks = 100,
     col = "gray90", border = "white", main = "",
     xlab = "", ylab = "", ylim = c(0, max_y * 1.05))
curve(dgamma(x, shape = gamma_shape, rate = gamma_rate), col = "blue", lwd = 2, add = TRUE)
curve(dlnorm(x, meanlog = mu_ln, sdlog = sqrt(sigma2_ln)), col = "green", lwd = 2, add = TRUE)
curve(dnorm(x, mean = normal_mean, sd = normal_sd), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Gamma", "Lognormal", "Normal"),
       col = c("blue", "green", "red"), lwd = 2)

# Summary statistics
median(aggregate_claims_pois_lognorm)
quantile(aggregate_claims_pois_lognorm, probs = c(0.90, 0.95, 0.99))
max(aggregate_claims_pois_lognorm)

# ===================================================
# NEW DATA: NEGATIVE BINOMIAL-LOGNORMAL MODEL
# ===================================================

model_code_nb_lognorm <- "
model {
  for (i in 1:N_y) {
    y[i] ~ dlnorm(mu_lognorm, tau_lognorm)
  }

  for (i in 1:N_n) {
    n[i] ~ dnegbin(p_nb, r_nb)
  }

  # Priors for Lognormal part
  mu_lognorm ~ dnorm(0, 1.0E-8)
  sigma_lognorm ~ dgamma(1, 0.0001)
  tau_lognorm <- 1 / pow(sigma_lognorm, 2)

  # Priors for Negative Binomial part
  p_nb ~ dbeta(1, 1)
  r_nb ~ dgamma(0.01, 0.01)
}
"

# MLE estimates for initialization
log_y <- log(itamtplcost_input$y)
mu_hat_lognorm <- mean(log_y)
sigma_hat_lognorm <- sd(log_y)

x_bar_n <- mean(itamtplcost_input$n)
var_n <- var(itamtplcost_input$n)
r_hat_nb <- x_bar_n^2 / (var_n - x_bar_n)
p_hat_nb <- r_hat_nb / (r_hat_nb + x_bar_n)

# Initial values for the model
inits_nb_lognorm <- list(
  list(mu_lognorm = 0.01, sigma_lognorm = 0.5, p_nb = 0.2, r_nb = 5, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 111),
  list(mu_lognorm = 2, sigma_lognorm = 1.5, p_nb = 0.5, r_nb = 10, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 222),
  list(mu_lognorm = mu_hat_lognorm, sigma_lognorm = sigma_hat_lognorm, p_nb = p_hat_nb, r_nb = r_hat_nb, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 789)
)

# ---------------------------------------------------
# RUN JAGS MODEL
# ---------------------------------------------------
jags_model_nb_lognorm <- jags.model(
  textConnection(model_code_nb_lognorm),
  data = itamtplcost_input,
  inits = inits_nb_lognorm,
  n.chains = 3,
  n.adapt = 20000
)

samples_nb_lognorm <- coda.samples(
  jags_model_nb_lognorm,
  variable.names = c("mu_lognorm", "tau_lognorm", "r_nb", "p_nb"),
  n.iter = 3000000,
  thin = 100
)

# Diagnostics & summaries
params <- c("mu_lognorm", "tau_lognorm", "r_nb", "p_nb")

# Gelman-Rubin diagnostics
for (param in params) {
  gelman.plot(samples_nb_lognorm[, param], xlab = "", ylab = "")
}
gelman.diag(samples_nb_lognorm, autoburnin = FALSE)

# Autocorrelation diagnostics
autocorr.diag(samples_nb_lognorm[, "mu_lognorm"], lags = 1:10)
autocorr.diag(samples_nb_lognorm[, "tau_lognorm"], lags = 1:10)
autocorr.diag(samples_nb_lognorm[, "r_nb"], lags = 1:10)
autocorr.diag(samples_nb_lognorm[, "p_nb"], lags = 1:10)

# Traceplots
traceplot(samples_nb_lognorm[, params], main="", xlab="")

# Posterior summary
summary(samples_nb_lognorm)

# Density plots
plot(density(as.matrix(samples_nb_lognorm[, "mu_lognorm"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples_nb_lognorm[, "tau_lognorm"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples_nb_lognorm[, "r_nb"])), main="", xlab="", ylab="")
plot(density(as.matrix(samples_nb_lognorm[, "p_nb"])), main="", xlab="", ylab="")

# Lognormal expected value calculation
posterior_nb_lognorm <- as.matrix(samples_nb_lognorm)

expected_y_lognorm <- lognormal_ev(posterior_nb_lognorm[, "mu_lognorm"], posterior_nb_lognorm[, "tau_lognorm"])
mean(expected_y_lognorm)
sd(expected_y_lognorm)
quantile(expected_y_lognorm, c(0.025, 0.5, 0.975), na.rm = TRUE)

plot(density(expected_y_lognorm, na.rm=TRUE), main="", xlab="", ylab="")

# ---------------------------------------------------
# PLOT: Empirical vs Posterior Lognormal CDF
# ---------------------------------------------------
y_obs <- itamtplcost_input$y
mu_post_mean <- mean(posterior_nb_lognorm[, "mu_lognorm"])
sigma_post_mean <- mean(1 / sqrt(posterior_nb_lognorm[, "tau_lognorm"]))

tau_post_mean <- mean(posterior_pois_lognorm[, "tau_pois_lognorm"])

plot(
  ecdf(y_obs),
  col = "blue",
  main = "",
  xlab = "",
  ylab = "",
  lwd = 2
)

x_vals <- seq(min(y_obs), max(y_obs), length.out = 500)
lines(x_vals, lognormal_cdf(x_vals, mu_post_mean, sigma_post_mean), col = "red", lwd = 2, lty = 2)

legend("bottomright",
       legend = c("Empirical CDF", sprintf("Lognormal(%.3f, %.3f) CDF", mu_post_mean, tau_post_mean)),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# ---------------------------------------------------
# PLOT 2: EMPIRICAL VS POSTERIOR NEGATIVE BINOMIAL CDF
# ---------------------------------------------------
n <- itamtplcost_input$n
r_post_mean <- mean(posterior_nb_lognorm[, "r_nb"])
p_post_mean <- mean(posterior_nb_lognorm[, "p_nb"])

# Plot empirical CDF for Negative Binomial counts
plot(
  ecdf(n),
  col = "blue",
  main = "",
  xlab = "",
  ylab = "",
  lwd = 2,
  xlim = c(min(n), max(n) + 40)  # extend x-axis for visibility
)

# Overlay theoretical Negative Binomial CDF
x_vals <- min(n):(max(n) + 40)
lines(x_vals, pnbinom(x_vals, r_post_mean, p_post_mean),
      col = "red", lwd = 2, lty = 2, type = "s")  # 's' for step

legend("bottomright",
       legend = c("Empirical CDF", sprintf("NegBinomial(%.3f, %.3f) CDF", r_post_mean, p_post_mean)),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# ---------------------------------------------------
# Posterior Predictive Sampling for y
# ---------------------------------------------------
lognorm_predictions <- numeric(1000)
uniform_samples <- runif(1000)

mu_posterior <- posterior_nb_lognorm[, "mu_lognorm"]
sigma_posterior <- 1 / sqrt(posterior_nb_lognorm[, "tau_lognorm"])

for (i in seq_along(uniform_samples)) {
  U <- uniform_samples[i]
  
  normal_samples <- qnorm(U, mean = mu_posterior, sd = sigma_posterior)
  lognorm_samples <- exp(normal_samples)
  
  lognorm_predictions[i] <- mean(lognorm_samples)
}

plot(density(lognorm_predictions, na.rm = TRUE), main = "", xlab = "", ylab = "")

# ---------------------------------------------------
# PREDICTIONS: Negative Binomial-Lognormal Distribution
# ---------------------------------------------------
set.seed(123)

# Function to simulate Negative Binomial samples
negbin <- function(r, p, U) {
  k <- 0
  prob <- dnbinom(k, size = r, prob = p)
  F <- prob
  
  while (U > F) {
    k <- k + 1
    prob <- dnbinom(k, size = r, prob = p)
    F <- F + prob
  }
  
  return(k)
}

# Simulation parameters
n_sim <- 1000

r_posterior <- posterior_nb_lognorm[, "r_nb"]
p_posterior <- posterior_nb_lognorm[, "p_nb"]

# Simulate number of claims
uniform_nb <- runif(n_sim)
nb_results <- numeric(n_sim)

for (i in seq_along(uniform_nb)) {
  U <- uniform_nb[i]
  
  cat("Simulation", i, "\n")
  
  nb_samples <- numeric(length(r_posterior))
  
  for (j in seq_along(r_posterior)) {
    nb_samples[j] <- negbin(r = r_posterior[j], p = p_posterior[j], U = U)
  }
  
  nb_results[i] <- round(mean(nb_samples))
}

# Compare empirical vs simulated claim counts
summary(itamtplcost_input$n) # Empirical
summary(nb_results) # Fitted

# Simulate aggregate claims
aggregate_claims_nb_lognorm <- numeric(length(nb_results))

for (i in seq_along(nb_results)) {
  n <- nb_results[i]
  
  cat("Simulation", i, "\n")
  
  lognorm_samples <- numeric(n)
  uniform_lognorm <- runif(n)
  
  for (j in seq_along(uniform_lognorm)) {
    U <- uniform_lognorm[j]
    
    normal_samples <- qnorm(U, mean = mu_posterior, sd = sigma_posterior)
    lognorm_samples[j] <- mean(exp(normal_samples))
  }
  
  aggregate_claims_nb_lognorm[i] <- sum(lognorm_samples)
}

# Compare to empirical aggregate claim distribution
summary(help$sum) # Empirical
summary(aggregate_claims_nb_lognorm) # Fitted

# ---------------------------------------------------
# Fit & Plot Distributions to Aggregate Claims
# ---------------------------------------------------

# Estimate distribution parameters via moment matching
mean_S <- mean(aggregate_claims_nb_lognorm)
var_S <- var(aggregate_claims_nb_lognorm)

gamma_shape <- mean_S^2 / var_S
gamma_rate <- mean_S / var_S

sigma2_ln <- log(1 + var_S / mean_S^2)
mu_ln <- log(mean_S) - sigma2_ln / 2

normal_mean <- mean_S
normal_sd <- sqrt(var_S)

# Define support for density plots
x_vals <- seq(min(aggregate_claims_nb_lognorm), max(aggregate_claims_nb_lognorm), length.out = 1000)

gamma_densities <- dgamma(x_vals, shape = gamma_shape, rate = gamma_rate)
lognorm_densities <- dlnorm(x_vals, meanlog = mu_ln, sdlog = sqrt(sigma2_ln))
normal_densities <- dnorm(x_vals, mean = normal_mean, sd = normal_sd)

hist_density <- hist(aggregate_claims_nb_lognorm, plot = FALSE, probability = TRUE)$density
max_y <- max(gamma_densities, lognorm_densities, normal_densities, hist_density)

# Plotting
hist(aggregate_claims_nb_lognorm, probability = TRUE, breaks = 100,
     col = "gray90", border = "white", main = "",
     xlab = "", ylab = "", ylim = c(0, max_y * 1.05))
curve(dgamma(x, shape = gamma_shape, rate = gamma_rate), col = "blue", lwd = 2, add = TRUE)
curve(dlnorm(x, meanlog = mu_ln, sdlog = sqrt(sigma2_ln)), col = "green", lwd = 2, add = TRUE)
curve(dnorm(x, mean = normal_mean, sd = normal_sd), col = "red", lwd = 2, add = TRUE)
legend("topright", legend = c("Gamma", "Lognormal", "Normal"),
       col = c("blue", "green", "red"), lwd = 2)

# Summary statistics
median(aggregate_claims_nb_lognorm)
quantile(aggregate_claims_nb_lognorm, probs = c(0.90, 0.95, 0.99))
max(aggregate_claims_nb_lognorm)
