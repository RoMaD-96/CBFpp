data {
  int<lower=0> n;                // Number of trials in current data
  int<lower=0, upper=n> y;       // Number of successes in current data
  real<lower=0> alpha_h;         // Alpha for historical Beta prior (from historical data)
  real<lower=0> beta_h;          // Beta for historical Beta prior (from historical data)
  real<lower=0> alpha_v;         // Alpha for vague Beta prior
  real<lower=0> beta_v;          // Beta for vague Beta prior
  real<lower=0, upper=1> mix_weight;      // Mixture weight between historical and vague prior
}

parameters {
  real<lower=0, upper=1> theta;              // Probability of success
  // real<lower=0, upper=1> mix_weight;      // Mixture weight between historical and vague prior
}

model {
  // Mixture prior for theta
  target += log_mix(mix_weight,
                    beta_lpdf(theta | alpha_h, beta_h),  
                    beta_lpdf(theta | alpha_v, beta_v)); 
  // Likelihood for observed data
  target += binomial_lpmf(y | n, theta);
}
