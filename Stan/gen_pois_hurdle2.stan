functions {
  real genpoiss_truncated_lpmf(int y, real theta, real lambda, int truncation) {
    if (y > truncation) {
      return negative_infinity();
    }

    real prob = theta * pow(theta + lambda * y, y - 1) * exp(-theta - lambda * y) / tgamma(y + 1);

    if (prob == 0) {
      return -1000;
    } else {
      real log_prob = log(prob);
      real Z = 0;
      for (i in 0:truncation) {
        Z += theta * pow(theta + lambda * i, i - 1) * exp(-theta - lambda * i) / tgamma(i + 1);
      }
      return log_prob - log(Z);
    }
  }
}
data {
  int<lower=1> N;  // Total number of observations
  int place[N];  
  int unit_length[N];
}

parameters {
  real<lower = 0> theta;
  real<lower = 0, upper = 1> lambda;
  real<lower = 0> mu; 
  real<lower = 0> phi;
  real psi_intercept;    // Intercept for psi logits
  real psi_slope;        // Slope for psi logits
  real alpha_intercept;  // Intercept for alpha logits
  real alpha_slope;      // Slope for alpha logits
}

transformed parameters {
  real<lower=0, upper=1> psi[N];
  real<lower=0, upper=1> alpha[N];
  
  real lprior = 0;
  lprior += gamma_lpdf(theta | 2, 0.5);
  lprior += gamma_lpdf(lambda | 1, 1);
  lprior += gamma_lpdf(mu | 8, 1);
  lprior += gamma_lpdf(phi | 2, 0.25);
  lprior += normal_lpdf(alpha_intercept | 0, 5);  // Prior for alpha_intercept
  lprior += normal_lpdf(alpha_slope | 0, 1);      // Prior for alpha_slope
  lprior += normal_lpdf(psi_intercept | 0, 5);    // Prior for psi_intercept
  lprior += normal_lpdf(psi_slope | 0, 1);        // Prior for psi_slope

  for (i in 1:N) {
    vector[3] logits;
    
    // Compute logits as a function of unit_length[i]
    logits[1] = psi_intercept + psi_slope * unit_length[i];   // Logit for psi
    logits[2] = alpha_intercept + alpha_slope * unit_length[i]; // Logit for alpha
    logits[3] = 0;  // Logit for the remaining probability

    // Apply softmax transformation
    vector[3] probs = softmax(logits);
    psi[i] = probs[1];
    alpha[i] = probs[2];
  }
}

model {
  // Priors
  target += lprior;

  // Likelihood
  for (i in 1:N) {
    if (place[i] == 0) {
      target += log(alpha[i]);  // Probability of place = unit_length[i]
    } 
    else if (place[i] == 1) {
      target += log(psi[i]);  // Probability of place = unit_length[i] - 1
    } 
    else {
      real lpos_0 = (theta * pow(theta, -1) * exp(-theta)) / tgamma(1);
      real lpos_1 = (theta * exp(-theta - lambda)) / tgamma(2);
      target += log(1 - psi[i] - alpha[i])
                + genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i])
                - log(1 - lpos_0 - lpos_1);
    }

    target += neg_binomial_2_lpmf(unit_length[i] | mu, phi);  // Negative binomial for unit_length
  }
}

generated quantities {
  real log_lik[N];

  for (i in 1:N) {
    if (place[i] == 0) {
      log_lik[i] = log(alpha[i]);
    } 
    else if (place[i] == 1) {
      log_lik[i] = log(psi[i]);
    } 
    else {
      real lpos_0 = (theta * pow(theta, -1) * exp(-theta)) / tgamma(1);
      real lpos_1 = (theta * exp(-theta - lambda)) / tgamma(2);
      log_lik[i] = log(1 - psi[i] - alpha[i])
                  + genpoiss_truncated_lpmf(place[i] | theta, lambda, unit_length[i])
                  - log(1 - lpos_0 - lpos_1);
    }

    log_lik[i] += neg_binomial_2_lpmf(unit_length[i] | mu, phi);
  }
}
