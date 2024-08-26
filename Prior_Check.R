# Prior Predictive Check 

############ Generalized Poisson Simple, Hurdle

# Load necessary libraries
library(MASS)  
library(HHMpa)

n <- 2000

lambda_samples <- numeric(0)

while (length(lambda_samples) < n) {
  new_samples <- rgamma(n, shape = 1, rate = 1)
  new_samples <- new_samples[new_samples <= 0.9]
  lambda_samples <- c(lambda_samples, new_samples)
}

lambda_samples <- lambda_samples[1:n]

# Check the result
length(lambda_samples)  # Should be 2000
summary(lambda_samples)  

theta_samples <- rgamma(n, shape = 2, rate = 0.5)
mu_samples <- rgamma(n, shape = 1, rate = 1)
phi_samples <- rgamma(n, shape = 1, rate = 1)

neg_binom_data <- c()
gen_pois_data <- c()
for (i in 1:2000) {
  print(i)
  neg_binom <- rnbinom(100, size = phi_samples[i], prob = phi_samples[i] / (phi_samples[i] + mu_samples[i]))
  gen_pois <- rgenpois(100, theta_samples[i], lambda_samples[i])
  neg_binom_data <- c(neg_binom_data, neg_binom)
  gen_pois_data <- c(gen_pois_data, gen_pois)
}

filtered_vector_binom <- neg_binom_data[neg_binom_data < 20]
filtered_vector_pois <- gen_pois_data[gen_pois_data < 20]

# Visualize the simulated data
hist(neg_binom_data, main = "Prior Predictive Check: Negative Binomial", xlab = "Simulated Data", breaks = 20)
hist(gen_pois_data, main = "Prior Predictive Check: Generalized Poisson", xlab = "Simulated Data", breaks = 20)


############ Mixture Model










