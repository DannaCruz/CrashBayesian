data {
  int<lower=1> N;
  int<lower=1> K;
  array[N] int<lower=0> Y;
  matrix[N, K] X;
  vector[N] log_offset;
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  int<lower=1> N_components;
  array[N_components + 1] int<lower=1, upper=N + 1> component_start;
  array[N] int<lower=1, upper=N> component_nodes;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> phi;
  vector[N] u;
  real<lower=0> sigma_u;
}
transformed parameters {
  vector[N] eta = log_offset + alpha + X * beta + sigma_u * u;
}
model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 2.5);
  phi ~ exponential(1);
  sigma_u ~ normal(0, 1);

  // Proper graph-CAR prior with precision I + D - W. The independent
  // N(0,1) term makes the precision nonsingular for disconnected graphs and
  // singleton components; the edge penalty smooths adjacent network units.
  u ~ normal(0, 1);
  target += -0.5 * dot_self(u[node1] - u[node2]);

  Y ~ neg_binomial_2_log(eta, phi);
}
