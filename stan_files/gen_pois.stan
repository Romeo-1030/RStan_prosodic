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
  lprior += gamma_lpdf(mu | 1, 1);
  lprior += gamma_lpdf(phi | 1, 1);
}
model {
  target += lprior;
  for (i in 1:N) {
    target += neg_binomial_2_lpmf(unit_length[i] | mu, phi);
    target += genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i]); // Modeling unit_length variable
  }
  
}

generated quantities {
  real log_lik[N];
  for (i in 1:N){
    log_lik[i] = neg_binomial_2_lpmf(unit_length[i] | mu, phi);
    log_lik[i] += genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i]);
  }
}


