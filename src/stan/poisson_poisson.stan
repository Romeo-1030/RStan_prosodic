// Poisson IU length and Poisson Place truncated at length_minus_1.
data {
  int<lower=1> N;
  int<lower=0> place[N];
  int<lower=0> unit_length[N];
}

parameters {
  real<lower=0> theta;
  real<lower=0> mu;
}

model {
  theta ~ gamma(2, 0.5);
  mu ~ gamma(8, 1);

  for (i in 1:N) {
    unit_length[i] ~ poisson(mu);
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
    log_lik[i] = poisson_lpmf(unit_length[i] | mu)
                 + poisson_lpmf(place[i] | theta)
                 - poisson_lcdf(unit_length[i] | theta);
  }
}
