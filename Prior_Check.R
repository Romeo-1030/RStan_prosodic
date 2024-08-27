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


# Number of prior predictive samples
n_prior_samples <- 1000

# Simulate theta1, theta2, mu1, phi1, psi from Gamma distributions as before
theta1_samples <- rgamma(n_prior_samples, shape = 16, rate = 4)
theta2_samples <- rgamma(n_prior_samples, shape = 16, rate = 4)
mu1_samples <- rgamma(n_prior_samples, shape = 8, rate = 8)
phi1_samples <- rgamma(n_prior_samples, shape = 8, rate = 8)
psi_samples <- runif(n_prior_samples, 0, 1)

# Simulate lambda1, lambda2 within [0, 1] using Beta distribution
lambda1_samples <- rbeta(n_prior_samples, shape1 = 8, shape2 = 8)
lambda2_samples <- rbeta(n_prior_samples, shape1 = 8, shape2 = 8)

# Define the number of observations
N <- 100

# Placeholder for simulated data
simulated_place <- matrix(NA, n_prior_samples, N)
simulated_back_place <- matrix(NA, n_prior_samples, N)
simulated_unit_length <- matrix(NA, n_prior_samples, N)

# Loop over prior samples to generate simulated data
for (i in 1:n_prior_samples) {
  for (n in 1:N) {
    # Simulate unit length
    simulated_unit_length[i, n] <- rnbinom(1, size = phi1_samples[i], 
                                           prob = phi1_samples[i]/(phi1_samples[i] + mu1_samples[i]))
    
    # Simulate place and back_place based on mixture
    if (runif(1) < psi_samples[i]) {
      simulated_place[i, n] <- rgenpois(1, theta1_samples[i], lambda1_samples[i])
    } else {
      simulated_back_place[i, n] <- rgenpois(1, theta2_samples[i], lambda2_samples[i])
    }
  }
}

# Plot histograms of the simulated data
hist(simulated_unit_length, breaks = 30, main = "Prior Predictive Check - Unit Length", xlab = "Unit Length")
hist(simulated_place, breaks = 30, main = "Prior Predictive Check - Place", xlab = "Place")
hist(simulated_back_place, breaks = 30, main = "Prior Predictive Check - Back Place", xlab = "Back Place")









