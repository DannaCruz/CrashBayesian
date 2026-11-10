data {
 int<lower=1> N;
 int<lower=1> K;
 array[N] int<lower=0> Y;
 matrix[N,K] X;
 vector[N] log_offset;
 int<lower=0> N_edges;
 array[N_edges] int<lower=1,upper=N> node1;
 array[N_edges] int<lower=1,upper=N> node2;
 vector[3] global_location;
 matrix[3,3] global_L;
 vector[N] spatial_location;
 matrix[N,3] spatial_B;
 vector<lower=0>[N] spatial_sd;
 vector<lower=0>[N] q_diag;
 vector<lower=0>[N] mu_ref;
 vector<lower=0>[N] reference_precision;
}
parameters {
 vector[3] global_raw;
 vector[K] beta;
 vector[N] spatial_raw;
}
transformed parameters {
 vector[3] global_theta = global_location + global_L * global_raw;
 real alpha = global_theta[1];
 real<lower=0> phi = exp(global_theta[2]);
 real<lower=0> sigma_u = exp(global_theta[3]);
 vector[N] spatial_effect;
 real log_spatial_jacobian;
 vector[N] eta;
 {
  vector[N] local_scale = spatial_sd .* sqrt(reference_precision ./
    (q_diag / square(sigma_u) + (mu_ref * phi) ./ (mu_ref + phi)));
  spatial_effect = spatial_location + spatial_B * global_raw + local_scale .* spatial_raw;
  log_spatial_jacobian = sum(log(local_scale));
 }
 eta = log_offset + alpha + X * beta + spatial_effect;
}
model {
 alpha ~ normal(0,5);
 beta ~ normal(0,2.5);
 phi ~ exponential(1);
 sigma_u ~ normal(0,1);
 // Block triangular Jacobian: exp transforms and parameter-dependent
 // diagonal local scales. Only det(global_L) is a fixed omitted constant.
 target += global_theta[2] + global_theta[3] + log_spatial_jacobian;
 target += -N * log(sigma_u) - 0.5 / square(sigma_u) *
   (dot_self(spatial_effect) + dot_self(spatial_effect[node1]-spatial_effect[node2]));
 Y ~ neg_binomial_2_log(eta,phi);
}
