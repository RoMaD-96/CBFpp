#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "doFuture",
  "parallel",
  "doParallel",
  "foreach",
  "rstan",
  "bridgesampling",
  "future",
  "dplyr",
  "progressr",
  "ggplot2",
  "ggridges",
  "tidyr"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


source("R/CBF_Functions.R")


#   ____________________________________________________________________________
#   Grid Search                                                             ####

eta <- round(seq(0.5, 6, length.out = 12),2)
nu <- round(seq(0.5, 6, length.out = 12),2)
grid <- as.data.frame(expand.grid(x = eta, y = nu))
model_id <- seq(1,nrow(grid))
grid <- cbind(grid,model_id)
model_id <- model_id[-14]
grid <- grid[-14,]
combinations <- expand.grid(model_id,14)


#   ____________________________________________________________________________
#   Simulated "True" Datasets                                               ####


##  ............................................................................
##  Historical Data                                                         ####

set.seed(442)

true.lambda_0 <- 2
N_0 <- 200
y_0 <- rpois(N_0, lambda = true.lambda_0)


#   ____________________________________________________________________________
#   Current Data                                                            ####

set.seed(441)
true.lambda <- 2.2
N <- 200
y <- rpois(N, lambda = true.lambda)



#   ____________________________________________________________________________
#   Observed Bayes Factor                                                   ####


##  ............................................................................
##  Parallel Processing                                                     ####

stan_model <- list("R/Simulations/Poisson/STAN/Poisson_Norm.stan")


# Hyperparameters

alpha_0 <- 1
beta_0 <- 1


# Create lists
random_data_obs_bf_num <- random_data_obs_bf_den <-  list()
for (i in 1:nrow(combinations)) {
  random_data_obs_bf_num[[i]] <- list(
    N0 = N_0,
    y0 = y_0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    eta = grid$x[i],
    nu = grid$y[i],
    N = length(y),
    y = as.array(y)
  )
  
  random_data_obs_bf_den[[i]] <- list(
    N0 = N_0,
    y0 = y_0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    eta = 1,
    nu = 1,
    N = length(y),
    y = as.array(y)
  )
}

# Define number of iterations
n_it_obs_bf <- nrow(combinations)


# Set up parallel backend with number of cores to use
cores <- detectCores()

# Register the parallel backend
plan(multicore, workers = 7)

# Run the loop in parallel
opts <- list(packages = c("bridgesampling","rstan"),
             seed = TRUE)


# Observed Bayes Factor
obs_bf_list <- list()
with_progress({
  p <- progressor(along = 1:n_it_obs_bf)
  system.time( obs_bf_list <- foreach(i = 1:n_it_obs_bf,.options.future = opts) %dofuture% {
    
    p()
    models_num <- unnormalized_posterior_sampl(random_data_obs_bf_num[i],
                                               stan_model,
                                               iterations = 2000,
                                               ref_rate = 0,
                                               info_model = FALSE)
    models_den <- unnormalized_posterior_sampl(random_data_obs_bf_den[i],
                                               stan_model,
                                               iterations = 2000,
                                               ref_rate = 0,
                                               info_model = FALSE)
    ### Bridge sampling and Bayes Factor
    bridge_num <- bridge_sampler(models_num[[1]], silent = TRUE)
    bridge_den <- bridge_sampler(models_den[[1]], silent = TRUE)
    
    bf(bridge_num,bridge_den, log = TRUE)
  })
})


obs_bf <- c()
for (j in 1:nrow(combinations)) {
  obs_bf[j] <- obs_bf_list[[j]]$bf
}


df_obs_bf <- cbind(obs_bf,combinations)

plan(sequential)

#   ____________________________________________________________________________
#   Simulations                                                             ####

print("Generating Posterior Samples")

##  ............................................................................
##  Generating Quantities                                                   ####

stan_model_gen <- list("R/Simulations/Poisson/STAN/Poisson_Norm_Gen.stan")

random_data_obs_bf <- list()

for (i in 1:nrow(combinations)) {
  random_data_obs_bf[[i]] <- list(
    N0 = N_0,
    y0 = y_0,
    alpha0 = alpha_0,
    beta0 = beta_0,
    eta = grid$x[i],
    nu = grid$y[i],
    N = length(y),
    y = as.array(y)
  )
}

models_gen <- unnormalized_posterior_sampl(random_data_obs_bf,
                                           stan_model_gen,
                                           iterations = 2000)





# Extract data from each model 
models_gen_par <- lapply(models_gen, function(model) rstan::extract(model))

# Get y_rep values from the extracted data
models_gen_rep <- lapply(models_gen_par, function(model) model$y_rep[3000:4000,])



rm(models_gen, models_gen_par)

##  ............................................................................
##  Parallel Processing                                                     ####

print("Doing Parallel Simulations")

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Simulations Parallel M_1                                                ####

stan_model_rep <- list("R/Simulations/Poisson/STAN/Poisson_Norm_Rep.stan")

# Create lists
bf_list <- list()
# Define number of iterations
n_iterations <- nrow(combinations)

# Set up parallel backend with number of cores to use
cores <- detectCores()

# Register the parallel backend
plan(multicore, workers = 7)

# Run the loop in parallel
opts <- list(progress = progress,
             packages = c("bridgesampling","rstan"),
             seed = TRUE)
options(future.globals.maxSize= 891289600)
# Parallel Computing
num_rep <- 1:200

with_progress({
  p <- progressor(along = 1:n_iterations)
  system.time( bf_list <- foreach(j = 1:n_iterations,.options.future = opts) %dofuture% {
    
    p()
    
    y_rep <- lapply(num_rep, function(i) models_gen_rep[[j]][i,])
    
    random_data_num <- lapply(y_rep, function(y_val) {
      list(
        N0 = N_0,
        y0 = y_0,
        alpha0 = alpha_0,
        beta0 = beta_0,
        N = length(y_val),
        y_rep = as.array(y_val),
        eta = grid$x[j],
        nu = grid$y[j]
      )
    })
    
    random_data_den <- lapply(y_rep, function(y_val) {
      list(
        N0 = N_0,
        y0 = y_0,
        alpha0 = alpha_0,
        beta0 = beta_0,
        N = length(y_val),
        y_rep = as.array(y_val),
        eta = 1,
        nu = 1
      )
    })
    
    models_num <- unnormalized_posterior_sampl(random_data_num,
                                               stan_model_rep,
                                               iterations = 2000,
                                               ref_rate = 0,
                                               info_model = FALSE)
    
    models_den <- unnormalized_posterior_sampl(random_data_den,
                                               stan_model_rep,
                                               iterations = 2000,
                                               ref_rate = 0,
                                               info_model = FALSE)
    
    bridge_num <- lapply(models_num, 
                         function(model) bridge_sampler(model, silent = TRUE))
    bridge_den <- lapply(models_den,
                         function(model) bridge_sampler(model, silent = TRUE))
    
    bayesF <- list()
    
    for (k in 1:length(num_rep)) {
      bayesF[[k]] <- bf(bridge_num[[k]], bridge_den[[k]], log = TRUE)
    }
    
    bayesF
    
  })
})


plan(sequential)

bf_data <- data.frame()
for (i in 1:nrow(combinations)) {
  for (j in 1:200) {
    bf_data[i,j] <- bf_list[[i]][[j]]$bf
  }
}



save.image(file = "R/Simulations/Poisson/RData/Poisson_Sim_2_2_20.RData")
