// functions {
//   real genpoiss_truncated_lpmf(int y, real theta, real lambda, int truncation) {
// 
//     if (y > truncation) {
//       return -1000; // Probability is zero outside the truncation
//     }
//     if ((theta * pow(theta + lambda * y, y-1) * exp(-theta - lambda * y)) / tgamma(y + 1) == 0) {
//       return -1000;
//     }
//     return log((theta * pow(theta + lambda * y, y-1) * exp(-theta - lambda * y)) / tgamma(y + 1));
//   }
// }
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
}

transformed parameters {
  real lprior = 0;
  lprior += gamma_lpdf(theta | 2, 0.5); 
  lprior += gamma_lpdf(lambda | 1, 1);
  lprior += gamma_lpdf(mu | 8, 1);
  lprior += gamma_lpdf(phi | 2, 0.25);
}

model {
  target += lprior;
  
  // Loop through all observations
  for (i in 1:N) {
    target += neg_binomial_2_lpmf(unit_length[i] | mu, phi);
    target += genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i]); // Modeling unit_length variable
  }
}

generated quantities {
  real log_lik[N];
  
  // Loop to calculate log likelihood for each observation
  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_lpmf(unit_length[i] | mu, phi);
    log_lik[i] += genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i]);
  }
}


