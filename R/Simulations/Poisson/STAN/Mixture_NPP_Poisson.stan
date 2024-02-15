data{
  int<lower=0> N0;
  int y0[N0];
  real<lower=0> alpha0;
  real<lower=0> beta0;
  int<lower=0> N;
  int y[N];
  int n_mixt_comp;
  real<lower=0> eta[n_mixt_comp];
  real<lower=0> nu[n_mixt_comp];
}
transformed data{
  real s = sum(y0);
  real logPprime = 0.0;
  for(i in 1:N0) logPprime += lgamma(y0[i] + 1);
}

parameters{
  simplex[n_mixt_comp] theta;
  real<lower=0> lambda;
  real<lower=0, upper=1> delta;
}


model{
  
  /* Mixture */
  theta ~ dirichlet(rep_vector(2.0, n_mixt_comp));
  vector[n_mixt_comp] contributions_delta;
    for(i in 1:N) {
    for(k in 1:n_mixt_comp) {
      contributions_delta[k] = log(theta[k]) + beta_lpdf(delta | eta[k], nu[k]);
    }
    target += log_sum_exp(contributions_delta);
  }
  
  
  
  // target += beta_lpdf(delta | eta, nu);
  /* Power prior */
  target += -( -contributions_delta * logPprime + lgamma(contributions_delta * s + alpha0) - (contributions_delta * s + alpha0)* log(delta * N0 + beta0) );
  target += contributions_delta * poisson_lpmf(y0 | lambda) ;
  target += gamma_lpdf(lambda | alpha0, beta0);
  
  /* Likelihood */
  target += poisson_lpmf(y | lambda);
}
