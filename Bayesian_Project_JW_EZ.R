rytgaard1990_input <- list(
  y = c(2.495, 2.120, 2.095, 1.700, 1.650, 1.985, 1.810, 1.625, 3.215, 2.105, 1.765, 1.715, 19.180, 1.915, 1.790, 1.755),
  n = c(5, 3, 4, 0, 4)
)

rytgaard1990_input$N_y <- length(rytgaard1990_input$y)
rytgaard1990_input$N_n <- length(rytgaard1990_input$n)
rytgaard1990_input$min_y <- min(rytgaard1990_input$y)

# Load a dataset from CASdatasets (CASdatasets requires zoo, xts, and Rtools)
#install.packages("C:/Users/Julia/Downloads/CASdatasets_1.2-0.tar.gz", repos = NULL, type = "source")
library(CASdatasets)
library(dplyr)

data(itamtplcost)
# https://dutangc.github.io/CASdatasets/reference/itamtplcost.html#ref-usage
# This dataset contains large losses (in excess of 500 Keuro) of an Italian Motor-TPL company since 1997.

itamtplcost$year <- format(as.Date(itamtplcost$Date, format = "%d/%m/%Y"), "%Y")
itamtplcost <- itamtplcost %>%
  rename(claim_amount = UltimateCost) %>%
  select(year, claim_amount)

library(rjags)
library(coda)

itamtplcost_input <- list(
  y = itamtplcost$claim_amount,
  n = itamtplcost %>% group_by(year) %>% summarise(count=n()) %>% pull(count)
)

itamtplcost_input$N_y <- length(itamtplcost_input$y)
itamtplcost_input$N_n <- length(itamtplcost_input$n)
itamtplcost_input$min_y <- min(itamtplcost_input$y)

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

inits_list <- list(
  list(alpha = 0.00001, beta = 0.00001, theta = 0.00001),   # Chain 1
  list(alpha = 100000, beta = 1, theta = 100000), # Chain 2
  list(alpha = 3.076, beta = 1.625, theta = 3.2)  # Chain 3
)

jags_model <- jags.model(
  textConnection(model_code),
  data = rytgaard1990_input,
  inits = inits_list,
  n.chains = 3,
  n.adapt = 20000
  )

#update(jags_model, 20000)

samples <- coda.samples(
  jags_model,
  variable.names = c("alpha", "beta", "theta"),
  n.iter = 50000
  )
summary(samples)
#plot(samples)
