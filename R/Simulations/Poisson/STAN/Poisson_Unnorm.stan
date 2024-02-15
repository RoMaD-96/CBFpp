data{
  int<lower=0> N0;
  int y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  int<lower=0> N;
  int y[N];
  real<lower=0> eta;
  real<lower=0> nu;
}
parameters{
  real<lower=0> lambda;
  real<lower=0, upper=1> delta;
}
model{
  /* Power prior */
  target += delta * poisson_lpmf(y0 | lambda) ;
  target += gamma_lpdf(lambda | alpha0, beta0);
  target += beta_lpdf(delta | eta, nu);
  /* Likelihood */
  target += poisson_lpmf(y | lambda);
}
