
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

// summation of psi 1 and psi 2 <= 1
parameters {
  real<lower = 0> theta;
  real<lower = 0, upper = 1> lambda;
  real<lower = 0> mu; 
  real<lower = 0> phi;
  real psi_intercept;  // Intercept of the psi linear function
  real psi_slope;      // Slope of the psi linear function
  real <lower = 0, upper = 1> alpha;
}

transformed parameters {
  real lprior = 0;
  lprior += gamma_lpdf(theta | 2, 0.5); // Specify your prior distribution for theta
  lprior += gamma_lpdf(lambda | 1, 1);//  Specify your prior distribution for lambda
  lprior += gamma_lpdf(mu | 1, 1);
  lprior += gamma_lpdf(phi | 1, 1);
}

model {
  target += lprior;
  for (i in 1:N) {
    real psi = inv_logit(psi_intercept + psi_slope * unit_length[i]); 
    
    if (place[i] == 0) {
      target += log(alpha);
    } 
    else if (place[i] == 1) {
      target += log(psi);
    }
    else {
      real lpos_0 = (theta * pow(theta, -1) * exp(-theta)) / tgamma(1);
      real lpos_1 = (theta * exp(-theta - lambda)) / tgamma(2);
      target += log(1 - psi- alpha) 
                + genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i]) 
                - log(1 - lpos_0 - lpos_1);
    }
    target += neg_binomial_2_lpmf(unit_length[i] | mu, phi);
  }
}


generated quantities {
  real log_lik[N];

  for (i in 1:N){
    real psi = inv_logit(psi_intercept + psi_slope * unit_length[i]);
    if (place[i] == 0) {
      log_lik[i] = log(alpha);
    }
    else if (place[i] == 1) {
      log_lik[i] = log(psi);
    }
    else {
      real lpos_0 = (theta * pow(theta, -1) * exp(-theta)) / tgamma(1);
      real lpos_1 = (theta * exp(-theta - lambda)) / tgamma(2);
      log_lik[i] = log(1 - psi- alpha) 
                + genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i]) 
                - log(1 - lpos_0 - lpos_1);
    }
    log_lik[i] += neg_binomial_2_lpmf(unit_length[i] | mu, phi);
  }
}


