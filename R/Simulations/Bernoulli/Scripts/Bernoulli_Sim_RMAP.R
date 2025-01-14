#   ____________________________________________________________________________
#   Libraries                                                               ####

library(rstan)



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

data_mix <- data.frame()
par_seq <-  seq(20,45, by = 1)

for (M in 1:length(par_seq)) {
  
  cat("Doing theta =",
      par_seq[M],
      "\n")
  
  
##  ............................................................................
##  RMAP models                                                             ####
  
  set.seed(4231)
  
  y <- par_seq[M]
  
  rmix_model_list <- list(
      n = N,
      y = par_seq[M],
      alpha_h = alpha_h,
      beta_h = beta_h,
      alpha_v = alpha_v,
      beta_v = beta_v,
      mix_weight = mix_weight)
    
  mix_model <- stan(file = "R/Simulations/Bernoulli/Stan/RMAP.stan",
                  data = rmix_model_list,
                  iter = 2000,
                  control = list(adapt_delta = 0.95, max_treedepth = 15),
                  refresh = 0, seed = 433)
    
  mix_model_par <- rstan::extract(mix_model)
    
  mean_theta <- mean(mix_model_par[["theta"]])
  sd_theta <- sd(mix_model_par[["theta"]])
    
  data_iter <- t(c(mean_theta, sd_theta))
  data_iter <- as.data.frame(data_iter)
  colnames(data_iter) <- c("mean_theta_mix", "sd_theta_mix")
    
  data_mix <- rbind.data.frame(data_mix,data_iter)
}

# Create the file name
save(data_mix, file = "R/Simulations/Bernoulli/RData/RMAP_SAM_prior/RMAP_prior.RData")

data_mix