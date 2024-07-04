#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "rstan",
  "npowerPrioR",
  "progressr",
  "ggplot2",
  "readxl",
  "dplyr",
  "rstanarm",
  "bayesplot",
  "mgcv")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

source("R/CBF_Functions.R")

#   ____________________________________________________________________________
#   Data                                                                    ####


library(hdbayes)

data("E2696")
# data("E1694")

historical_data <- E2696




#   ____________________________________________________________________________
#   STAN Parameter Configuration                                            ####

N_0 <- nrow(historical_data)
X_0 <- cbind(historical_data$age,  
            historical_data$treatment, 
            historical_data$sex, 
            historical_data$perform 
            
)
#   ____________________________________________________________________________
#   Power Prior Sampling Estimates                                          ####


historical_data_norm_pp <- list(
  N0 = N_0,
  P = ncol(X_0),
  y0 = historical_data$failind,
  X0 = X_0,
  a_0 = NULL
)

##  ............................................................................
##  Grid and Step size                                                      ####

J <- 20
maxA <- 1
epsilon <- 0.05


##  ............................................................................
##  Estimates                                                               ####

prior <- stan_model("STAN/Normalizing_Constant_PP_Logit.stan")

time_eval <- system.time(
  delta_estimates <- build_grid(
    compiled.model.prior = prior,
    eps = epsilon,
    M = maxA,
    J = J,
    v1 = 10,
    v2 = 10,
    stan.list = historical_data_norm_pp,
    pars = c("alpha", "beta")
  )
)


write.csv(
  delta_estimates$result,
  row.names = FALSE,
  file = paste0(
    "Data/delta_estimates_logit.csv"
  )
)
