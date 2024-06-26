#   ____________________________________________________________________________
#   STAN Data Block Configuration                                           ####



N_0 <- length(log_age_hist)
X_0 <- cbind(log_age_hist, # x1 - log age 
             historical_data$sex, # x2 - gender
             historical_data$trt # x3 treatment
)
Y_0_cens <- historical_data$survtime
Cens_0 <- historical_data$scens


N <- N_0
X <- X_0
Y_cens <- historical_data$survtime
Cens<- historical_data$scens

# 
# N <- length(log_age_current)
# X <- cbind(rep(1, N), # x0 - dummy for intercept
#            log_age_current, # x1 - log age 
#            current_data$sex, # x2 - gender 
#            current_data$trt # x3 treatment
# )
# Y_cens <- current_data$survtime
# Cens<- current_data$scens

random_data_obs_bf_num <- list(
  N0 = N_0,
  P = ncol(X_0),
  X0 = X_0,
  y0 = Cens_0,
  sex_0 = historical_data$sex,
  perform_0 = historical_data$perform,
  N = N,
  X = X,
  y = Cens,
  sex = historical_data$sex,
  perform = historical_data$perform,
    eta = 20,
    nu = 1
  )


random_data_obs_bf_den <- list(
    N0 = N_0,
    P = ncol(X_0),
    X0 = X_0,
    y0 = Cens_0,
    sex_0 = historical_data$sex,
    perform_0 = historical_data$perform,
    N = N,
    X = X,
    y = Cens,
    sex = historical_data$sex,
    perform = historical_data$perform,
    eta = 1, 
    nu = 1    
  )










p()
### Approximation Stuff
Ks <- c(10000)

J <- 20

### GAM Approximation 
constant_data <- read.csv("Data/delta_estimates_logit.csv")
fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

### Observed BF

models_num <- norm_post_pp(Ks, 
                           list(random_data_obs_bf_num), 
                           stan_model)
models_den <- norm_post_pp(Ks,
                           list(random_data_obs_bf_den),
                           stan_model)
### Bridge sampling and Bayes Factor
bridge_num <- bridge_sampler(models_num[[1]], silent = TRUE, method = "warp3", use_neff = FALSE)
bridge_den <- bridge_sampler(models_den[[1]], silent = TRUE, method = "warp3", use_neff = FALSE)

bf_val <- bf(bridge_num,bridge_den, log = TRUE)








ncores = max(1, parallel::detectCores() - 1)
warmup  = 1000          ## warmup for MCMC sampling
total.samples = 10000   ## number of samples post warmup
samples = ceiling(warmup + total.samples / ncores)  ## outputs approx total.samples samples




library(parallel)
a0     = seq(0, 1, length.out = ncores * 5)

## wrapper to obtain log normalizing constant in parallel package
logncfun = function(a0, ...)
  hdbayes::glm.npp.lognc(
    formula = formula, family = family, histdata = histdata, a0 = a0, ...
  )

cl = makeCluster(ncores)
clusterSetRNGStream(cl, 123)
clusterExport(cl, varlist = c('formula', 'family', 'histdata'))

## call created function
a0.lognc = parLapply(
  cl = cl, X = a0, fun = logncfun, iter = 5000, warmup = warmup, refresh = 0
)
stopCluster(cl)

a0.lognc = data.frame( do.call(rbind, a0.lognc) )
head(a0.lognc)










