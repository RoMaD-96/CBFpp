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
transformed data{
  real s = sum(y0);
  real logPprime = 0.0;
  for(i in 1:N0) logPprime += lgamma(y0[i] + 1);
}
parameters{
  real<lower=0> lambda;
  real<lower=0, upper=1> delta;
}
model{
  /* Power prior */
   target += -( -delta * logPprime + lgamma(delta * s + alpha0) - (delta * s + alpha0)* log(delta * N0 + beta0) );
  target += delta * poisson_lpmf(y0 | lambda) ;
  target += gamma_lpdf(lambda | alpha0, beta0);
  target += beta_lpdf(delta | eta, nu);
  /* Likelihood */
  target += poisson_lpmf(y | lambda);
}
generated quantities{
  int y_rep [N];
  for (n in 1:N)
    y_rep[n] = poisson_rng(lambda);
}
