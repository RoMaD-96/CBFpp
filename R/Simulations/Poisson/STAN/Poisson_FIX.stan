data{
  int<lower=0> N0;
  int y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  int<lower=0> N;
  real<lower=0, upper=1> delta;
}

parameters{
  real<lower=0> lambda;
}

model{
  /* Power prior */
  target += delta * poisson_lpmf(y0 | lambda) ;
  target += gamma_lpdf(lambda | alpha0, beta0);
}

generated quantities{
  int y_curr [N];
  for (n in 1:N)
    y_curr[n] = poisson_rng(lambda);
}
