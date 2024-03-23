data {
  int<lower=0> N0;
  int<lower=0> P;
  matrix[N0, P] X0;
  vector[N0] sex_0;
  vector[N0] treat_0;
  int<lower=0, upper=1> y0[N0];
  real<lower=0> delta;
}
parameters {
  real alpha;
  vector[P] beta;
  real beta_inter;
}
transformed parameters{
  real logL = bernoulli_logit_lpmf(y0 | alpha + X0*beta +sex_0.*treat_0*beta_inter );
  real logL_sq = square(logL);
}
model {
  /*power prior*/
  target += normal_lpdf(alpha | 0, 1);
  target += normal_lpdf(beta | 0, 1);
  target += normal_lpdf(beta_inter | 0, 1);
  target += delta * logL;
}
