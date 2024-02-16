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
  "patchwork",
  "ggplot2",
  "ggridges",
  "tidyr",
  "readxl"
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
#   Datasets                                                                ####

data <- drop_na(read_excel("Data/trial_e1684_e1690_Merged.xlsx", 
                           col_types = c("numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric")))

historical_data <- filter(data, study == 1684)[,-2]
current_data <- filter(data, study == 1690)[,-2]

## Standardization ##

log_age_hist <- log(historical_data$age)
log_age_current <- log(current_data$age)


#   ____________________________________________________________________________
#   Approximately Normalized Prior                                          ####

### Grid Size ###
J <- 20

### GAM Approximation ###
constant_data <- read.csv("Data/delta_estimates_logit.csv")
fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

#   ____________________________________________________________________________
#   STAN Data Block Configuration                                           ####



N_0 <- length(log_age_hist)
X_0 <- cbind(log_age_hist, # x1 - log age 
             historical_data$sex, # x2 - gender
             historical_data$trt # x3 treatment
)
Y_0_cens <- historical_data$survtime
Cens_0 <- historical_data$scens


N<- length(log_age_current)
X <- cbind(rep(1, N), # x0 - dummy for intercept
           log_age_current, # x1 - log age 
           current_data$sex, # x2 - gender 
           current_data$trt # x3 treatment
)
Y_cens <- current_data$survtime
Cens<- current_data$scens

##  ............................................................................
##  Parallel Processing                                                     ####


stan_model <- list("R/Melanoma/STAN/Logit_Sim_Norm_PP.stan")

# Hyperparameters

alpha_0 <- 1
beta_0 <- 1


#   ____________________________________________________________________________
#   Observed Bayes Factor                                                  ####


random_data_obs_bf_num <- lapply(1:nrow(combinations), function(i) {
  list(
    N0 = N_0,
    P = ncol(X_0),
    X0 = X_0,
    y0 = Cens_0,
    N = N,
    X = X,
    y = Cens,
    eta = grid$x[i],
    nu = grid$y[i]
  )
})

random_data_obs_bf_den <- lapply(1:nrow(combinations), function(i) {
  list(
    N0 = N_0,
    P = ncol(X_0),
    X0 = X_0,
    y0 = Cens_0,
    N = N,
    X = X,
    y = Cens,
    eta = 1, 
    nu = 1    
  )
})


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
    ### Approximation Stuff
    Ks <- c(10000)

    J <- 20

    ### GAM Approximation 
    constant_data <- read.csv("Data/delta_estimates_logit.csv")
    fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

    ### Observed BF

    models_num <- norm_post_pp(Ks, 
                               random_data_obs_bf_num[i], 
                               stan_model)
    models_den <- norm_post_pp(Ks,
                               random_data_obs_bf_den[i],
                               stan_model)
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

stan_model_gen <- list("R/Melanoma/STAN/Logit_Sim_Norm_PP_Gen.stan")

random_data_obs_bf <- lapply(1:nrow(combinations), function(i) {
  list(
    N0 = N_0,
    P = ncol(X_0),
    X0 = X_0,
    y0 = Cens_0,
    N = N,
    X = X,
    y = Cens,
    eta = grid$x[i],
    nu = grid$y[i]
  )
})

Ks <- c(10000)
J <- 20

models_gen <- norm_post_pp(Ks,
                           random_data_obs_bf,
                           stan_model_gen)

# Extract data from each model 
models_gen_par <- lapply(models_gen, function(model) rstan::extract(model))

# Get y_rep values from the extracted data
models_gen_rep <- lapply(models_gen_par, function(model) model$y_rep[3000:4000,])

rm(models_gen, models_gen_par)




##  ............................................................................
##  BF Replications                                                         ####

stan_model_rep <- list("R/Melanoma/STAN/Logit_Sim_Norm_PP_Rep.stan")

# Function for parallel computing
code_block <- function(portion_rep, iteration) {
  
  print("Doing Parallel CBF Simulations")
  
  # Create lists
  bf_list <- list()
  # Define number of iterations
  n_iterations <- nrow(combinations)
  my_seq <- portion_rep[[iteration]]
  
  # Set up parallel backend with number of cores to use
  cores <- detectCores()
  
  # Register the parallel backend
  plan(multicore, workers = 7)
  
  # Run the loop in parallel
  opts <- list(progress = progress,
               packages = c("bridgesampling","rstan","parallel"),
               seed = TRUE)
  options(future.globals.maxSize= 891289600)
  
  with_progress({
    p <- progressor(along = 1:n_iterations)
    system.time( bf_list <- foreach(j = 1:n_iterations,.options.future = opts) %dofuture% {
      
      p()
      y_rep <- lapply(my_seq, function(i) models_gen_rep[[j]][i,])
      
      random_data_num <- lapply(y_rep, function(y_val) {
        list(
          N0 = N_0,
          P = ncol(X_0),
          X0 = X_0,
          y0 = Cens_0,
          N = length(y_val),
          X = X,
          y_rep = as.array(y_val),
          eta = grid$x[j],
          nu = grid$y[j]
        )
      })
      
      random_data_den <- lapply(y_rep, function(y_val) {
        list(
          N0 = N_0,
          P = ncol(X_0),
          X0 = X_0,
          y0 = Cens_0,
          N = length(y_val),
          X = X,
          y_rep = as.array(y_val),
          eta = 1,
          nu = 1
        )
      })
      
      ### Approximation Stuff
      Ks <- c(10000)
      
      J <- 20
      
      ### GAM Approximation 
      constant_data <- read.csv("delta_estimates_logit.csv")
      fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)
      
      ### BF Replications
      models_num <- norm_post_pp(Ks,
                                 random_data_num,
                                 stan_model_rep)
      
      models_den <-  norm_post_pp(Ks,
                                  random_data_den,
                                  stan_model_rep)
      
      bridge_num <- lapply(models_num, 
                           function(model) bridge_sampler(model, silent = TRUE))
      bridge_den <- lapply(models_den,
                           function(model) bridge_sampler(model, silent = TRUE))
      
      bayesF <- list()
      
      for (k in 1:length(my_seq)) {
        bayesF[[k]] <- bf(bridge_num[[k]], bridge_den[[k]], log = TRUE)
      }
      rm(models_num, models_den, bridge_num, bridge_den, constant_data, fit_gam)
      gc()
      
      bayesF
      
    })
  })
  
  plan(sequential)

  filename <- sprintf("bf_list_%d.RData", iteration)
  save(bf_list, file = filename)
  
  gc()
  return(bf_list)
}


iteration_by_block <- lapply(seq(1, 100, by = 2), function(x) c(x:(x + 1)))


# Parallel Simulations by Block

BF_rep_dist <- lapply(1:length(iteration_by_block), function(iteration)  code_block(iteration_by_block, iteration))
