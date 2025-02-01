#   ____________________________________________________________________________
#   Libraries                                                               ####

library(rstan)
library(RBesT)
library(SAMprior)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

#   ____________________________________________________________________________
#   Data                                                                    ####

# Historical Data Parameters
N_0 <- 100
y_0 <- 20
p_0 <- y_0 / N_0

# Current Data Parameters
N <- 100

# Mixture priors parameters
alpha_h <- y_0 + 1
beta_h <- N_0 + 1 - y_0
alpha_v <- 1
beta_v <- 1
mix_weight <- 0.5

# Stan model
stan_mix_prior <- stan_model("R/Simulations/Bernoulli/Stan/RMAP.stan")


#   ____________________________________________________________________________
#   Simulation Setting                                                      ####

data_SAM <- data.frame()
par_seq <- seq(20, 45, by = 1)

for (M in 1:length(par_seq)) {
  cat(
    "Doing theta =",
    par_seq[M],
    "\n"
  )

  ##  ............................................................................
  ##  SAM prior                                                               ####

  set.seed(4231)

  y <- par_seq[M]

  prior_historical <- mixbeta(inf = c(1, alpha_h, beta_h))

  # Calculate the mixture weight of the SAM prior
  wSAM <- SAM_weight(
    if.prior = prior_historical,
    delta = 0.1, # Clinically significant difference
    n = 100,
    r = par_seq[M]
  )

  print(wSAM)

  ##  ............................................................................
  ##  SAM models                                                              ####

  rmix_model_list <- list(
    n = N,
    y = par_seq[M],
    alpha_h = alpha_h,
    beta_h = beta_h,
    alpha_v = alpha_v,
    beta_v = beta_v,
    mix_weight = wSAM
  )

  mix_model <- stan(
    file = "R/Simulations/Bernoulli/Stan/RMAP.stan",
    data = rmix_model_list,
    iter = 2000,
    control = list(adapt_delta = 0.95, max_treedepth = 15),
    refresh = 0, seed = 433
  )

  mix_model_par <- rstan::extract(mix_model)

  mean_theta <- mean(mix_model_par[["theta"]])
  sd_theta <- sd(mix_model_par[["theta"]])

  data_iter <- t(c(mean_theta, sd_theta))
  data_iter <- as.data.frame(data_iter)
  colnames(data_iter) <- c("mean_theta_SAM", "sd_theta_SAM")

  data_SAM <- rbind.data.frame(data_SAM, data_iter)
}

# Create the file name
save(data_SAM, file = "R/Simulations/Bernoulli/RData/RMAP_SAM_prior/SAM_prior.RData")

data_SAM