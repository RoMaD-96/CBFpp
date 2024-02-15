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
constant_data <- read.csv("delta_estimates_logit.csv")
fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

# mgcv::plot.gam(fit_gam)

#   ____________________________________________________________________________
#   STAN Data Block Configuration                                           ####


N_0 <- length(log_age_hist)
X_0 <- cbind(intercept = rep(1, N_0), # x0 - dummy for intercept
             log_age = log_age_hist, # x1 - (log) age 
             sex = historical_data$sex, # x2 - gender (sex)
             perf_status = historical_data$perform # x3 performance status
)
Y_0_cens <- historical_data$survtime
Cens_0 <- historical_data$scens


N<- length(log_age_current)
X <- cbind(rep(1, N), # x0 - dummy for intercept
           log_age_current, # x1 - (log) age 
           current_data$sex, # x2 - gender (sex)
           current_data$perform # x3 performance status
)
Y_cens <- current_data$survtime
Cens<- current_data$scens

##  ............................................................................
##  Parallel Processing                                                     ####


# stan_model <- list("C:/Users/romac/Desktop/Dottorato/Research_Project_PhD/R_Projects/Power_Priors/Bayes_Factor/STAN/Poisson_Models/Poisson_Norm.stan")
stan_model <- list("STAN/Logit_Sim_Norm_PP.stan")

# Hyperparameters

alpha_0 <- 1
beta_0 <- 1


#   ____________________________________________________________________________
#   Observed Bayes Factor                                                  ####



# Create lists
random_data_obs_bf_num <- list()
random_data_obs_bf_den <- list()
for (i in 1:nrow(combinations)) {
  random_data_obs_bf_num[[i]] <- list(
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

random_data_obs_bf_den[[i]] <- list(
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
    ### Approximation Stuff
    Ks <- c(5000)

    J <- 20

    ### GAM Approximation ###
    constant_data <- read.csv("delta_estimates_logit.csv")
    fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

    ###

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

stan_model_gen <- list("STAN/Logit_Sim_Norm_PP_Gen.stan")

random_data_obs_bf <- list()

for (i in 1:nrow(combinations)) {
  random_data_obs_bf[[i]] <- list(
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
}

Ks <- c(5000)
J <- 20


#options(mc.cores = 4)

models_gen <- norm_post_pp(Ks,
                           random_data_obs_bf,
                           stan_model_gen)

#options(mc.cores = 1)

#save.image(file = "Partial_Logit.RData")

# Extract data from each model 
models_gen_par <- lapply(models_gen, function(model) rstan::extract(model))

# Get y_rep values from the extracted data
models_gen_rep <- lapply(models_gen_par, function(model) model$y_rep[3000:4000,])



rm(models_gen, models_gen_par)

save.image(file = "Logit_VM.RData")
# 
# ##  ............................................................................
# ##  Parallel Processing                                                     ####
# 
# print("Doing Parallel Simulations")
# 
# ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
# ### Simulations Parallel M_1                                                ####
# 
# combinations <- rbind.data.frame(combinations, combinations)
# combinations <- combinations %>% filter(Var1 != Var2)
# 
# 
# stan_model_rep <- list("STAN/Logit_Sim_Norm_PP_Rep.stan")
# 
# # Create lists
# bf_list <- list()
# # Define number of iterations
# n_iterations <- nrow(combinations)
# 
# # Set up parallel backend with number of cores to use
# cores <- detectCores()
# 
# # Register the parallel backend
# plan(multicore, workers = 7)
# 
# # Run the loop in parallel
# opts <- list(progress = progress,
#              packages = c("bridgesampling","rstan"),
#              seed = TRUE)
# options(future.globals.maxSize= 891289600)
# with_progress({
#   p <- progressor(along = 1:n_iterations)
#   system.time( bf_list <- foreach(j = 1:n_iterations,.options.future = opts) %dofuture% {
#     
#     p()
#     y_rep <- list()
#     random_data_num <- list()
#     random_data_den <- list()
#     if(j<=(n_iterations/2)) {
#       for (i in 1:100) {
#         y_rep[[i]] <-models_gen_rep[[combinations$Var1[j]]][i,]
#         ### Models Grid
#         
#         random_data_num[[i]] <- list(
#           N0 = N_0,
#           P = ncol(X_0),
#           X0 = X_0,
#           y0 = Cens_0,
#           N = length(y_rep[[i]]),
#           X = X,
#           y_rep = as.array(y_rep[[i]]),
#           eta = grid$x[combinations$Var1[j]],
#           nu = grid$y[combinations$Var1[j]]
#         )
#         
#         random_data_den[[i]] <- list(
#           N0 = N_0,
#           P = ncol(X_0),
#           X0 = X_0,
#           y0 = Cens_0,
#           N = length(y_rep[[i]]),
#           X = X,
#           y_rep = as.array(y_rep[[i]]),
#           eta = grid$x[combinations$Var2[j]],
#           nu = grid$y[combinations$Var2[j]]
#         )
#       }
#     } else {
#       for (i in 1:100) {
#         y_rep[[i]] <- models_gen_rep[[combinations$Var2[j]]][i,]
#         ### Models Grid
#         
#         random_data_num[[i]] <- list(
#           N0 = N_0,
#           P = ncol(X_0),
#           X0 = X_0,
#           y0 = Cens_0,
#           N = length(y_rep[[i]]),
#           X = X,
#           y_rep = as.array(y_rep[[i]]),
#           eta = grid$x[combinations$Var1[j]],
#           nu = grid$y[combinations$Var1[j]]
#         )
#         
#         random_data_den[[i]] <- list(
#           N0 = N_0,
#           P = ncol(X_0),
#           X0 = X_0,
#           y0 = Cens_0,
#           N = length(y_rep[[i]]),
#           X = X,
#           y_rep = as.array(y_rep[[i]]),
#           eta = grid$x[combinations$Var2[j]],
#           nu = grid$y[combinations$Var2[j]]
#         )
#       }
#     }
#     
#     
#     ### Approximation Stuff
#     Ks <- c(10000)
#     
#     J <- 20
#     
#     ### GAM Approximation ###
#     constant_data <- read.csv("delta_estimates_logit.csv")
#     fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)
#     
#     models_num <- norm_post_pp(Ks, 
#                                random_data_num, 
#                                stan_model_rep)
#     
#     models_den <-  norm_post_pp(Ks,
#                                 random_data_den, 
#                                 stan_model_rep)
#     
#     bridge_num <- list()
#     bridge_den <- list()
#     for (k in 1:100) {
#       bridge_num[[k]] <- bridge_sampler(models_num[[k]], silent = TRUE)
#       bridge_den[[k]] <- bridge_sampler(models_den[[k]], silent = TRUE)
#     }
# 
#     bayesF <- list()
# 
#     for (k in 1:100) {
#       bayesF[[k]] <- bf(bridge_num[[k]], bridge_den[[k]], log = TRUE)
#     }
# 
#     
#     rm(data_list, models_num, models_den, bridge_num, bridge_den)
#     gc()
#     gc()
#     
#     bayesF
#     
#   })
# })
# 
# plan(sequential)
# 
# bf_data <- data.frame()
# for (i in 1:nrow(combinations)) {
#   for (j in 1:100) {
#     bf_data[i,j] <- bf_list[[i]][[j]]$bf
#   }
# }
# 
# 
# save.image(file = "Logit.RData")