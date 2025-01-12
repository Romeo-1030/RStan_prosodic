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

transformed data {
  int back_place[N];
  for (i in 1:N) {
    back_place[i] = unit_length[i] - place[i];
  }
}

parameters {
  real<lower = 0> theta1;
  real<lower = 0, upper = 1> lambda1;
  real<lower = 0> theta2;
  real<lower = 0, upper = 1> lambda2;
  real<lower = 0> mu; 
  real<lower = 0> phi;  
  real<lower = 0, upper = 1> psi;
}

transformed parameters {
  real lprior = 0;
  lprior += gamma_lpdf(theta1 | 2, 0.5); 
  lprior += beta_lpdf(lambda1 | 1, 1); // try something like (lambda + 1)/2 follows (support of 0-1) beta distribution
  lprior += gamma_lpdf(theta2 | 2, 0.5); 
  lprior += beta_lpdf(lambda2 | 1, 1); 
  lprior += gamma_lpdf(mu | 8, 1);
  lprior += gamma_lpdf(phi | 2, 0.25);
}

model {
  target += lprior;
  
  for (i in 1:N) {
  target += log_sum_exp(log(psi) + genpoiss_truncated_lpmf(place[i] | theta1, lambda1, unit_length[i]) 
                  + neg_binomial_2_lpmf(unit_length[i] | mu, phi),
                  log1m(psi) + genpoiss_truncated_lpmf(back_place[i] | theta2, lambda2, unit_length[i]) 
                  + neg_binomial_2_lpmf(unit_length[i] | mu, phi)); //change to the same negative binomial 
  }
}


generated quantities {
  real log_lik[N];
  for (i in 1:N){
    log_lik[i] = psi*(neg_binomial_2_lpmf(unit_length[i] | mu, phi)
    + genpoiss_truncated_lpmf(place[i] | theta1, lambda1, unit_length[i]));
    log_lik[i] += (1-psi)*(neg_binomial_2_lpmf(unit_length[i] | mu, phi)
    + genpoiss_truncated_lpmf(back_place[i] | theta2, lambda2, unit_length[i]));
  }
}
// 
// generated quantities {
//   real log_lik[N];
//   for (i in 1:N){
//     log_lik[i] = log_sum_exp(
//       log(psi) + genpoiss_truncated_lpmf(place[i] | theta1, lambda1, unit_length[i]) 
//                  + neg_binomial_2_lpmf(unit_length[i] | mu1, phi1),
//       log1m(psi) + genpoiss_truncated_lpmf(back_place[i] | theta2, lambda2, unit_length[i]) 
//                  + neg_binomial_2_lpmf(unit_length[i] | mu1, phi1)
//     );
//   }
// }



