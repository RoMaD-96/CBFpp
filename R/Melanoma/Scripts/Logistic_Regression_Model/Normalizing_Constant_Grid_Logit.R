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

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

source("R/CBF_Functions.R")

#   ____________________________________________________________________________
#   Data                                                                    ####


data <- read_excel("Data/trial_e1684_e1690_Merged.xlsx", 
                   col_types = c("numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric"))

historical_data <- filter(data, study == 1684)[,-2]
current_data <- filter(data, study == 1690)[,-2]

## Standardization ##

log_age_hist <- log(historical_data$age)
log_age_current <- log(current_data$age)


#   ____________________________________________________________________________
#   STAN Parameter Configuration                                            ####

N_0 <- length(log_age_hist)
X_0 <- cbind(log_age_hist, # x1 - log age 
             historical_data$sex, # x2 - gender
             historical_data$trt # x3 treatment
)
Y_0_cens <- historical_data$survtime
Cens_0 <- historical_data$scens


#   ____________________________________________________________________________
#   Power Prior Sampling Estimates                                          ####


historical_data_norm_pp <- list(
  N0 = N_0,
  P = ncol(X_0),
  y0 = Cens_0,
  X0 = X_0,
  sex_0 = historical_data$sex,
  treat_0 = historical_data$trt,
  delta = 0.05
)

##  ............................................................................
##  Grid and Step size                                                      ####

J <- 20
maxA <- 1
epsilon <- 0.05


##  ............................................................................
##  Estimates                                                               ####

prior <- stan_model("R/Melanoma/STAN/Normalizing_Constant_PP_Logit.stan")

time_eval <- system.time(
  delta_estimates <- build_grid(
    compiled.model.prior = prior,
    eps = epsilon,
    M = maxA,
    J = J,
    v1 = 10,
    v2 = 10,
    stan.list = historical_data_norm_pp,
    pars = c("alpha", "beta"),
    strict = FALSE
  )
)


write.csv(
  delta_estimates$result,
  row.names = FALSE,
  file = paste0(
    "Data/delta_estimates_logit.csv"
  )
)
