functions {
  real genpoiss_truncated_lpmf(int y, real theta, real lambda, int truncation) {
    
    // Check if y is greater than truncation (impossible case)
    if (y > truncation) {
      return negative_infinity();  // Return -inf for impossible values
    }

    // Compute the unnormalized probability for the given y
    real prob = theta * pow(theta + lambda * y, y - 1) * exp(-theta - lambda * y) / tgamma(y + 1);

    // Check if the probability is zero to avoid log(0)
    if (prob == 0) {
      return -1000;  // Return a very low log-probability for zero probability
    } else {
      real log_prob = log(prob);

      // Compute the normalization constant Z by summing over all valid y values (from 0 to truncation)
      real Z = 0;
      for (i in 0:truncation) {
        Z += theta * pow(theta + lambda * i, i - 1) * exp(-theta - lambda * i) / tgamma(i + 1);
      }

      // Return the normalized log-probability
      return log_prob - log(Z);  // Normalize by subtracting log(Z)
    }
  }
}

data {
  int<lower=1> N;  // total number of observations
  int place[N];  
  int unit_length[N];
}

parameters {
  real<lower = 0> theta;
  real<lower = 0, upper = 1> lambda;
  real<lower = 0> mu; 
  real<lower = 0> phi;
  real psi_intercept;  // Intercept of the psi linear function
  real psi_slope;      // Slope of the psi linear function
}

transformed parameters {
  real lprior = 0;
  lprior += gamma_lpdf(theta | 2, 0.5); // Prior for theta
  lprior += gamma_lpdf(lambda | 1, 1);  // Prior for lambda
  lprior += gamma_lpdf(mu | 1, 1);
  lprior += gamma_lpdf(phi | 1, 1);
}

model {
  target += lprior;
  for (i in 1:N) {
    real psi = inv_logit(psi_intercept + psi_slope * unit_length[i]); // Calculate psi for each observation

    // Hurdle at the last element: place[i] == unit_length[i]
    if (place[i] == unit_length[i]) {
      target += log(psi);  // Probability of place = unit_length[i]
    } 
    else {
      target += log(1 - psi);  // Probability of place < unit_length[i]
      target += genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i] - 1);  // Truncated Generalized Poisson for non-last element
    }
    target += neg_binomial_2_lpmf(unit_length[i] | mu, phi);  // Negative binomial for unit_length
  }
}

generated quantities {
  real log_lik[N];

  for (i in 1:N) {
    real psi = inv_logit(psi_intercept + psi_slope * unit_length[i]);

    // Log-likelihood calculation for the new hurdle condition
    if (place[i] == unit_length[i]) {
      log_lik[i] = log(psi);  // Log-probability of place = unit_length[i]
    } else {
      log_lik[i] = log(1 - psi) 
                + genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i] - 1);  // Truncated Generalized Poisson for non-last element
    }
    log_lik[i] += neg_binomial_2_lpmf(unit_length[i] | mu, phi);  // Negative binomial for unit_length
  }
}
