data {
  int<lower=1> N;
  int<lower=1> K;
  array[N] int<lower=0> Y;
  matrix[N, K] X;
  vector[N] log_offset;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> phi;
}
transformed parameters {
  vector[N] eta = log_offset + alpha + X * beta;
}
model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 2.5);
  phi ~ exponential(1);
  Y ~ neg_binomial_2_log(eta, phi);
}
