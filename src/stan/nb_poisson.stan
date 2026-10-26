// Negative binomial IU length and Poisson Place truncated at length_minus_1.
data {
  int<lower=1> N;
  int<lower=0> place[N];
  int<lower=0> unit_length[N];
}

parameters {
  real<lower=0> theta;
  real<lower=0> mu;
  real<lower=0> phi;
}

model {
  theta ~ gamma(2, 0.5);
  mu ~ gamma(8, 1);
  phi ~ gamma(2, 0.25);

  for (i in 1:N) {
    unit_length[i] ~ neg_binomial_2(mu, phi);
    if (place[i] > unit_length[i]) {
      target += negative_infinity();
    } else {
      target += poisson_lpmf(place[i] | theta)
                - poisson_lcdf(unit_length[i] | theta);
    }
  }
}

generated quantities {
  real log_lik[N];

  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_lpmf(unit_length[i] | mu, phi)
                 + poisson_lpmf(place[i] | theta)
                 - poisson_lcdf(unit_length[i] | theta);
  }
}
