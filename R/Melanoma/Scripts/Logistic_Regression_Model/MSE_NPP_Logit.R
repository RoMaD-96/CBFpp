#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "doFuture",
  "cmdstanr",
  "parallel",
  "doParallel",
  "hdbayes",
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

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


#   ____________________________________________________________________________
#   Data                                                                    ####

data("E2696")
data("E1694")
historical_data <- E2696
current_data <- E1694
all_data <- list(current_data, historical_data)



formula <- failind ~  age + treatment + sex + perform
p <- length(attr(terms(formula), "term.labels"))
family <- binomial("logit")



mean(historical_data$treatment)
mean(current_data$treatment)



a0 <- seq(0, 1, length.out = 21)

if (requireNamespace("parallel")) {
  if (instantiate::stan_cmdstan_exists()) {
    # wrapper to obtain log normalizing constant in parallel package
    logncfun <- function(a0, ...) {
      hdbayes::glm.npp.lognc(
        formula = formula, family = family, histdata = historical_data, a0 = a0, ...
      )
    }
    cl <- parallel::makeCluster(10)
    parallel::clusterSetRNGStream(cl, 123)
    parallel::clusterExport(cl, varlist = c("formula", "family", "historical_data"))
    # call created function
    time.npp.1 <- system.time(
      a0.lognc <- parLapply(
        cl = cl, X = a0, fun = logncfun, iter_warmup = 1000,
        iter_sampling = 1000, chains = 4, refresh = 0
      )
    )
    parallel::stopCluster(cl)
  }
  a0.lognc <- data.frame(do.call(rbind, a0.lognc))
}

#   ____________________________________________________________________________
#   Grid Search                                                             ####

eta <- round(seq(0.5, 6, length.out = 12), 2)
nu <- round(seq(0.5, 6, length.out = 12), 2)
grid <- as.data.frame(expand.grid(x = eta, y = nu))
model_id <- seq(1, nrow(grid))
grid <- cbind(grid, model_id)
model_id <- model_id[-14]
grid <- grid[-14, ]
combinations <- expand.grid(model_id, 14)


#   ____________________________________________________________________________
#   Models                                                                  ####


# model_list <- list()
# for (i in 1:nrow(grid)) {
#   model_list[[i]] <- glm.npp(
#     formula = formula, family = family, data.list = all_data,
#     a0.lognc = a0.lognc$a0,
#     lognc = matrix(a0.lognc$lognc, ncol = 1),
#     a0.shape1 = grid$x[i], a0.shape2 = grid$y[i],
#     iter_warmup = 1000, iter_sampling = 1000,
#     chains = 4, parallel_chains = 4,
#     refresh = 100, seed = 433)
# }
#
# save(model_list, file = "R/Melanoma/RData/model_list_MSE.RData")



#   ____________________________________________________________________________
#   MSE method Shen et al (2023)                                            ####


#   ____________________________________________________________________________
#   Derive simulated datasets                                               ####


fit_bhm <- glm.bhm(
  formula, family, list(historical_data),
  meta.mean.mean = 0, meta.mean.sd = 10,
  meta.sd.mean = 0, meta.sd.sd = 0.5,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 12345
)

beta_inter <- round(mean(fit_bhm$`(Intercept)`), 4)
beta_age_0 <- round(mean(fit_bhm$age), 4)
beta_treat_0 <- round(mean(fit_bhm$treatment), 4)
beta_sex_0 <- round(mean(fit_bhm$sex), 4)
beta_perform_0 <- round(mean(fit_bhm$perform), 4)


betas <- c(beta_inter, beta_treat_0, beta_sex_0, beta_age_0, beta_perform_0)
design_dataset <- cbind.data.frame(rep(1, nrow(current_data)), current_data[, 3:6])
design_matrix <- as.matrix(design_dataset)

eta <- design_matrix %*% betas
eta <- as.vector(eta)
p <- 1 / (1 + exp(-eta))

curr_sim_list <- list()

for (i in 1:50) {
  y_curr_sim <- rbinom(n = nrow(current_data), size = 1, prob = p)
  current_data_sim <- cbind.data.frame(failind = y_curr_sim, current_data[, 3:6])
  curr_sim_list[[i]] <- current_data_sim
}



model_list <- list()
# Define number of iterations
n_iterations <- nrow(grid)

# Set up parallel backend with number of cores to use
cores <- detectCores()

# Register the parallel backend
plan(multicore, workers = 7)

# Run the loop in parallel

opts <- list(
  packages = c("hdbayes"),
  seed = TRUE
)
# options(future.globals.maxSize= 891289600)

with_progress({
  p <- progressor(along = 1:n_iterations)
  system.time(model_list <- foreach(model = 1:n_iterations, .options.future = opts) %dofuture% {
    models <- list()

    for (i in 1:length(curr_sim_list)) {
      all_data_sim <- list(curr_sim_list[[i]], historical_data)

      model_sim <- glm.npp(
        formula = formula, family = family, data.list = all_data_sim,
        a0.lognc = a0.lognc$a0,
        lognc = matrix(a0.lognc$lognc, ncol = 1),
        a0.shape1 = grid$x[model], a0.shape2 = grid$y[model],
        iter_warmup = 1000, iter_sampling = 1000,
        chains = 4, parallel_chains = 1,
        refresh = 0, seed = 12345
      )


      models[[i]] <- model_sim
    }

    p()

    models
  })
})

plan(sequential)



#   ____________________________________________________________________________
#   When mu_star = 0 (treat+d_mtd)                                          ####

betas_alt <- c(beta_inter, 0, beta_sex_0, beta_age_0, beta_perform_0)
design_dataset_alt <- cbind.data.frame(rep(1, nrow(current_data)), current_data[, 3:6])
design_matrix_alt <- as.matrix(design_dataset_alt)

eta_alt <- design_matrix_alt %*% betas_alt
eta_alt <- as.vector(eta_alt)
p_alt <- 1 / (1 + exp(-eta_alt))

curr_sim_list_alt <- list()

for (i in 1:50) {
  y_curr_sim_alt <- rbinom(n = nrow(current_data), size = 1, prob = p_alt)
  current_data_sim_alt <- cbind.data.frame(failind = y_curr_sim_alt, current_data[, 3:6])
  curr_sim_list_alt[[i]] <- current_data_sim_alt
}



model_list_alt <- list()
# Define number of iterations
n_iterations <- nrow(grid)

# Set up parallel backend with number of cores to use
cores <- detectCores()

# Register the parallel backend
plan(multicore, workers = 7)

# Run the loop in parallel

opts <- list(
  packages = c("hdbayes"),
  seed = TRUE
)
# options(future.globals.maxSize= 891289600)

with_progress({
  p <- progressor(along = 1:n_iterations)
  system.time(model_list_alt <- foreach(model = 1:n_iterations, .options.future = opts) %dofuture% {
    models <- list()

    for (i in 1:length(curr_sim_list_alt)) {
      all_data_sim <- list(curr_sim_list_alt[[i]], historical_data)

      model_sim <- glm.npp(
        formula = formula, family = family, data.list = all_data_sim,
        a0.lognc = a0.lognc$a0,
        lognc = matrix(a0.lognc$lognc, ncol = 1),
        a0.shape1 = grid$x[model], a0.shape2 = grid$y[model],
        iter_warmup = 1000, iter_sampling = 1000,
        chains = 4, parallel_chains = 1,
        refresh = 0, seed = 12345
      )


      models[[i]] <- model_sim
    }

    p()

    models
  })
})

plan(sequential)


save(model_list_alt, file = "R/Melanoma/RData/model_list_MSE_alt.RData")





#   ____________________________________________________________________________
#   Computing MSE                                                           ####





compute_weighted_mse <- function(models, models_alt, y_0_bar, d_MTD, w) {
  # Initialize a vector to store weighted MSEs
  weighted_MSEs <- numeric(length(models))

  # Loop over each model
  for (i in 1:143) {
    posterior_treat_mean <- numeric(50)
    posterior_treat_mean_alt <- numeric(50)

    for (j in 1:50) {
      # Extract posterior mean for beta (treatment)
      posterior_treat_mean[j] <- mean(models[[i]][[j]]$treatment)
      posterior_treat_mean_alt[j] <- mean(models_alt[[i]][[j]]$treatment)
    }

    # Scenario 1: beta^* = y_0_bar
    beta_star1 <- y_0_bar
    MSE1 <- sum((posterior_treat_mean - beta_star1)^2) / 50

    # Scenario 2: beta^* = y_0_bar + d_MTD
    beta_star2 <- 0
    MSE2 <- sum((posterior_treat_mean_alt - beta_star2)^2) / 50
    # Compute weighted MSE
    weighted_MSEs[i] <- w * MSE1 + (1 - w) * MSE2
  }

  return(weighted_MSEs)
}

d_MTD <- -beta_treat_0
w <- 0.5


weighted_MSE_vec <- compute_weighted_mse(model_list, model_list_alt, beta_treat_0, d_MTD, w)

grid[which.min(weighted_MSE_vec), ]
