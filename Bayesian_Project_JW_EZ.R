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
  list(alpha = 0.5, beta = 1, theta = 5, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 123),
  list(alpha = 0.6, beta = 1.2, theta = 6, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 456),
  list(alpha = 0.4, beta = 0.9, theta = 4, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 789)
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
  xlab = "y",
  ylab = "",
  lwd = 2
)

# Overlay posterior predictive Pareto CDF
x_vals <- seq(min(y), max(y), length.out = 500)
lines(x_vals, pareto_cdf(x_vals, alpha, beta),
      col = "red", lwd = 2, lty = 2)

legend("bottomright",
       legend = c("Empirical CDF", "Pareto(3.079,1.592) CDF"),
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
  xlab = "n",
  ylab = "",
  lwd = 2,
  xlim = c(min(n), max(n) + 3)  # extend x-axis for visibility
)

# Overlay theoretical Poisson CDF
x_vals <- min(n):(max(n) + 3)
lines(x_vals, ppois(x_vals, lambda = theta),
      col = "red", lwd = 2, lty = 2, type = "s")  # 's' for step

legend("bottomright",
       legend = c("Empirical CDF", "Poisson(3.396) CDF"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# NEW
# Diagnostics

# Trace plots for each parameter
traceplot(samples[, "alpha"])
traceplot(samples[, "beta"])
traceplot(samples[, "theta"])

gelman.diag(samples, autoburnin=FALSE)

gelman.plot(samples[, "alpha"], xlab = "", ylab = "")
title(xlab = "Iterations", ylab = "PSRF")
gelman.plot(samples[, "beta"], xlab = "", ylab = "")
title(xlab = "Iterations", ylab = "PSRF")
gelman.plot(samples[, "theta"], xlab = "", ylab = "")
title(xlab = "Iterations", ylab = "PSRF")

# The average autocorrelation across chains
autocorr.diag(samples[, "alpha"], lags = 1:10)
autocorr.diag(samples[, "beta"], lags = 1:10)
autocorr.diag(samples[, "theta"], lags = 1:10)

autocorr.plot(samples[, "alpha"])
autocorr.plot(samples[, "beta"])
autocorr.plot(samples[, "theta"])

# Rerunning chains
jags_model_1 <- jags.model(
  textConnection(model_code),
  data = rytgaard1990_input,
  inits = inits_list,
  n.chains = 3,
  n.adapt = 20000
)

# Sample from posterior
samples_1 <- coda.samples(
  jags_model_1,
  variable.names = c("alpha", "beta", "theta"),
  n.iter = 300000,
  thin = 10
)

traceplot(samples_1[, "alpha"])
traceplot(samples_1[, "beta"])
traceplot(samples_1[, "theta"])

gelman.diag(samples_1, autoburnin=FALSE)

autocorr.diag(samples_1[, "alpha"], lags = 1:10)
autocorr.diag(samples_1[, "beta"], lags = 1:10)
autocorr.diag(samples_1[, "theta"], lags = 1:10)

autocorr.plot(samples_1[, "alpha"])
autocorr.plot(samples_1[, "beta"])
autocorr.plot(samples_1[, "theta"])

# Predictions
m <- nrow(posterior)

results <- numeric(15)
for (n in 0:14) {
  predictive_values <- (posterior[, "theta"]^n) * exp(-posterior[, "theta"]) / factorial(n)
  results[n + 1] <- mean(predictive_values)
}

