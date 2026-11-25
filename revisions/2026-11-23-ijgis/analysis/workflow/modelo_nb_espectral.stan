data {
  int<lower=1> N;
  int<lower=0> P;
  matrix[N,P] X;
  array[N] int<lower=0> y;
  vector[N] log_exposure;
  int<lower=0> M;
  matrix[N,M] B;
  vector<lower=-1,upper=1>[M] lambda;
  int<lower=0,upper=3> model_type; // 0 NB; 1 dense spectral; 2 selection; 3 spectral BYM2
  int<lower=0,upper=1> mixture;
  int<lower=0,upper=1> iid;
  real<lower=0,upper=1> gamma;
  real<lower=0> spike_sd;
  real<lower=0> field_normalizer;
}
transformed data {
  vector[M] weight;
  for(j in 1:M) weight[j] = inv_sqrt(1-gamma*lambda[j]);
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0> phi;
  vector[M] b;
  vector<lower=0>[model_type > 0 ? 1 : 0] tau;
  vector<lower=0,upper=1>[mixture] pi_mix;
  vector<lower=0,upper=1>[iid] rho;
  vector[iid * N] v_raw;
}
model {
  vector[N] field = rep_vector(0,N);
  if(model_type == 1) field = B * (tau[1] * weight .* b);
  if(model_type == 2) field = B * b;
  if(model_type == 3)
    field = tau[1] * (sqrt(rho[1]) * (B * (weight .* b)) / field_normalizer
                      + sqrt(1-rho[1]) * v_raw);
  alpha ~ normal(0,5);
  beta ~ normal(0,0.7);
  phi ~ exponential(0.5);
  if(model_type > 0) tau ~ normal(0,0.5);
  if(model_type == 1) b ~ std_normal();
  if(model_type == 2) {
    pi_mix ~ beta(1,4);
    for(j in 1:M) {
      // Strictly ordered component scales, unlike the original manuscript prior.
      real slab_sd = sqrt(square(spike_sd) + square(tau[1] * weight[j]));
      target += log_mix(pi_mix[1], normal_lpdf(b[j] | 0,slab_sd),
                                 normal_lpdf(b[j] | 0,spike_sd));
    }
  }
  if(model_type == 3) {
    b ~ std_normal();
    v_raw ~ std_normal();
    rho ~ beta(0.5,0.5);
  }
  y ~ neg_binomial_2_log_glm(X,log_exposure + alpha + field,beta,phi);
}
generated quantities {
  vector[M * mixture] q;
  real expected_K = 0;
  if(mixture == 1) {
    for(j in 1:M) {
      real slab_sd = sqrt(square(spike_sd) + square(tau[1] * weight[j]));
      real l1 = log(pi_mix[1]) + normal_lpdf(b[j] | 0,slab_sd);
      real l0 = log1m(pi_mix[1]) + normal_lpdf(b[j] | 0,spike_sd);
      q[j] = inv_logit(l1-l0);
    }
    expected_K = sum(q);
  }
}
