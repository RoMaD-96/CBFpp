data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  int<lower=0, upper=1> y[N];
  // int<lower=0> N;
  // matrix[N, P] X;
  // int<lower=0, upper=1> y[N];
  real<lower=0> eta;
  real<lower=0> nu;
}
parameters {
  real alpha;
  vector[P] beta;
  real<lower=0,upper=1> delta;
}
model {
  /* Power prior */
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta | 0, 1);
  target += beta_lpdf(delta | eta, nu);
  target += delta * bernoulli_logit_lpmf(y | alpha + X*beta);
  
  /* Likelihood */
  // target += bernoulli_logit_lpmf(y | alpha + X*beta);
}
