functions {
  real genpoiss_truncated_lpmf(int y, real theta, real lambda, int truncation) {
    
    if (y > truncation) {
      return -1000; // Probability is zero outside the truncation
    }
    if ((theta * pow(theta + lambda * y, y-1) * exp(-theta - lambda * y)) / tgamma(y + 1) == 0) {
      return -1000;
    }
    return log((theta * pow(theta + lambda * y, y-1) * exp(-theta - lambda * y)) / tgamma(y + 1));
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
  real<lower = -1, upper = 1> lambda1;
  real<lower = 0> theta2;
  real<lower = -1, upper = 1> lambda2;
  real<lower = 0> mu1; 
  real<lower = 0> phi1;  
  real<lower = 0, upper = 1> psi;
}

transformed parameters {
  real lprior = 0;
  lprior += gamma_lpdf(theta1 | 16, 4); 
  lprior += gamma_lpdf(lambda1 | 8, 8); // try something like (lambda + 1)/2 follows (support of 0-1) beta distribution
  lprior += gamma_lpdf(theta2 | 16, 4); 
  lprior += gamma_lpdf(lambda2 | 8, 8); 
  lprior += gamma_lpdf(mu1 | 8, 8);
  lprior += gamma_lpdf(phi1 | 8, 8);
}

model {
  target += lprior;
  
  for (i in 1:N) {
  target += log_sum_exp(log(psi) + genpoiss_truncated_lpmf(place[i] | theta1, lambda1, unit_length[i]) 
                  + neg_binomial_2_lpmf(unit_length[i] | mu1, phi1),
                  log1m(psi) + genpoiss_truncated_lpmf(back_place[i] | theta2, lambda2, unit_length[i]) 
                  + neg_binomial_2_lpmf(unit_length[i] | mu1, phi1)); //change to the same negative binomial 
  }
}


generated quantities {
  real log_lik[N];
  for (i in 1:N){
    log_lik[i] = psi*(neg_binomial_2_lpmf(unit_length[i] | mu1, phi1)
    + genpoiss_truncated_lpmf(place[i] | theta1, lambda1, unit_length[i]));
    log_lik[i] += (1-psi)*(neg_binomial_2_lpmf(unit_length[i] | mu1, phi1)
    + genpoiss_truncated_lpmf(back_place[i] | theta2, lambda2, unit_length[i]));
  }
}




