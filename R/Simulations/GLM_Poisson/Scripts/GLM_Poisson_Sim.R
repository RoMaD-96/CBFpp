#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "doFuture",
  "cmdstanr",
  "parallel",
  "doParallel",
  "hdbayes",
  "foreach",
  "future",
  "dplyr",
  "progressr",
  "rethinking",
  "tidyr",
  "readxl",
  "RBesT",
  "SAMprior"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

source("R/CBF_Functions.R", echo=TRUE)

#   ____________________________________________________________________________
#   Grid Search                                                             ####

eta <- round(seq(0.5, 6, by = 0.5),2)
nu <- round(seq(0.5, 6, by =0.5),2)
grid <- as.data.frame(expand.grid(x = eta, y = nu))
model_id <- seq(1,nrow(grid))
grid <- cbind(grid,model_id)
model_id <- model_id[-14]
grid <- grid[-14,]
combinations <- expand.grid(model_id,14)


# Simulation Scenarios

n_hist <- 500  # Historical data sample size
n_curr <- 500  # Current data sample size

# Historical data parameters
beta0_hist <- 0.3   # Intercept for historical data
beta1_hist <- 0.5   
beta2_hist <- 0.2   
beta3_hist <- -0.3  

# Simulate predictor variables for historical data
set.seed(442)
x1_hist <- rnorm(n_hist)
set.seed(443)
x2_hist <- rnorm(n_hist)
set.seed(444)
x3_hist <- rnorm(n_hist)

# Historical mean response
lambda_hist <- exp(beta0_hist + beta1_hist * x1_hist + beta2_hist * x2_hist + beta3_hist * x3_hist)

set.seed(433)
y_hist <- rpois(n_hist, lambda = lambda_hist)

historical_data <- data.frame(y = y_hist, 
                              x1 = x1_hist, 
                              x2 = x2_hist, 
                              x3 = x3_hist, 
                              lambda = lambda_hist)

# Current data 
# We'll vary beta1_curr while keeping beta0_curr = beta0_hist, beta2_curr = beta2_hist, beta3_curr = beta3_hist.

beta1_curr_values <- c(0.50, 0.53, 0.54, 0.55, 0.57, 0.58, 0.59,
                       0.60, 0.61, 0.62, 0.63, 0.64, 0.65)

beta2_curr <- beta2_hist
beta3_curr <- beta3_hist

# Simulate predictors for current data
set.seed(4231)
x1_curr <- rnorm(n_curr)
set.seed(4232)
x2_curr <- rnorm(n_curr)
set.seed(4233)
x3_curr <- rnorm(n_curr)

current_datasets <- list()

for (i in seq_along(beta1_curr_values)) {
  if (i == 1) {
    # First dataset matches historical data exactly
    curr_data <- historical_data
  } else {
    beta1_curr <- beta1_curr_values[i]
    
    # Mean response for current data scenario
    lambda_curr <- exp(beta0_hist + beta1_curr * x1_curr + beta2_curr * x2_curr + beta3_curr * x3_curr)  
    
    # Response variable for current data
    set.seed(433)
    y_curr <- rpois(n_curr, lambda = lambda_curr)
    
    curr_data <- data.frame(y = y_curr, 
                            x1 = x1_curr, 
                            x2 = x2_curr, 
                            x3 = x3_curr, 
                            lambda = lambda_curr)
  }
  
  current_datasets[[i]] <- curr_data
}

names(current_datasets) <- paste0("Disagreement_Level_", seq_along(beta1_curr_values))

# Mean of the historical data
mean_hist_y <- mean(historical_data$y)

# Mean of the current data for each scenario
mean_curr_y <- sapply(current_datasets, function(df) mean(df$y))
print(mean_curr_y)

#   ____________________________________________________________________________
#   Datasets                                                                ####

formula <- y ~ x1 + x2 + x3
p <- length(attr(terms(formula), "term.labels"))
family <- poisson(link = "log")

#   ____________________________________________________________________________
#   Approximate normalizing constant                                        ####

a0 <- seq(0, 1, length.out = 21)

if( requireNamespace("parallel") ){
  if (instantiate::stan_cmdstan_exists()) {
    # wrapper to obtain log normalizing constant in parallel package
    logncfun <- function(a0, ...){
      hdbayes::glm.npp.lognc(
        formula = formula, family = family, histdata = historical_data, a0 = a0, ...
      )
    }
    cl <- parallel::makeCluster(10)
    parallel::clusterSetRNGStream(cl, 123)
    parallel::clusterExport(cl, varlist = c('formula', 'family', 'historical_data'))
    # call created function
    time.npp.1 <- system.time(
      a0.lognc <- parLapply(
        cl = cl, X = a0, fun = logncfun, iter_warmup = 1000,
        iter_sampling = 1000, chains = 4, refresh = 0
      )
    )
    parallel::stopCluster(cl)
  }
  a0.lognc <- data.frame( do.call(rbind, a0.lognc) )
}

#   ____________________________________________________________________________
#   Simulation Studies                                                      ####

results_list <- list()

for (m in seq_along(current_datasets)) {
  cat("Processing current_datasets model:", m, "\n")
  
  current_data <- current_datasets[[m]]
  all_data <- list(current_data, historical_data)
  
##  ............................................................................
##  Parallel Processing                                                     ####
  

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Observed Bayes Factor                                                   ####
  
  model_unif_H0 <- glm.npp(
      formula = formula, family = family, data.list = all_data,
      a0.lognc = a0.lognc$a0,
      lognc = matrix(a0.lognc$lognc, ncol = 1),
      a0.shape1 = 1, a0.shape2 = 1, 
      iter_warmup = 1000, iter_sampling = 1000,
      chains = 4, parallel_chains = 4,
      refresh = 100, seed = 433
  )
  
  set.seed(4321)
  log_ML_H0 <- glm.logml.npp(model_unif_H0)
  
  # Define number of iterations
  n_it_obs_bf <- nrow(combinations)
  
  # Register the parallel backend
  cores <- detectCores() - 1
  plan(multisession, workers = cores)
  
  # Run the loop in parallel
  opts <- list(packages = c("hdbayes"),
               seed = TRUE)
  
  
  # Observed Bayes Factor
  
  obs_bf_list <- list()
  y_rep_list <- list()
  with_progress({
    p <- progressor(along = 1:n_it_obs_bf)
    system.time( obs_bf_list <- foreach(i = 1:n_it_obs_bf,.options.future = opts) %dofuture% {
      
      # Observed BF
      model_H1 <- glm.npp(
        formula = formula, family = family, data.list = all_data,
        a0.lognc = a0.lognc$a0,
        lognc = matrix(a0.lognc$lognc, ncol = 1),
        a0.shape1 = grid$x[i], a0.shape2 = grid$y[i], 
        iter_warmup = 1000, iter_sampling = 1000,
        chains = 4, parallel_chains = 1,
        refresh = 100, seed = 433
      )
      
      set.seed(4321)
      log_ML_H1 <- glm.logml.npp(model_H1)
      
      p()
      
      log_ML_H1$logml - log_ML_H0$logml
      
    })
  })
  
  
  obs_bf <- c()
  for (j in 1:nrow(combinations)) {
    obs_bf[j] <- obs_bf_list[[j]]
  }
  
  df_obs_bf <- cbind.data.frame(obs_bf,grid$x,grid$y)
  colnames(df_obs_bf) <- c("obs_bf", "eta_H1", "nu_H1")
  
  plan(sequential)
  
  
  print("Generating Posterior Samples")
  
##  ............................................................................
##  Generating Quantities                                                   ####
  
  df_obs_bf_pos <- filter(df_obs_bf, obs_bf>=0)
  
  n_iter <- nrow(df_obs_bf_pos)
  
  # Register the parallel backend
  plan(multisession, workers = cores)
  
  # Run the loop in parallel
  opts <- list(packages = c("hdbayes", "cmdstanr"),
               seed = TRUE)
  
  y_rep_list <- list()
  
  with_progress({
    p <- progressor(along = 1:n_iter)
    system.time( y_rep_list <- foreach(i = 1:n_iter,.options.future = opts) %dofuture% {
      
      # Observed BF
      source("R/Melanoma/Scripts/Logistic_Regression_Model/Fun_Post_Pred.R", echo=TRUE)
      
      model_H1 <- glm.npp_post_pred(
        formula = formula, family = family, data.list = all_data,
        a0.lognc = a0.lognc$a0,
        lognc = matrix(a0.lognc$lognc, ncol = 1),
        a0.shape1 = df_obs_bf_pos$eta_H1[i], a0.shape2 = df_obs_bf_pos$nu_H1[i], 
        iter_warmup = 1000, iter_sampling = 1000,
        chains = 4, parallel_chains = 1,
        refresh = 100, seed = 433
      )
      
      num_columns <- nrow(historical_data)
      
      columns_list <- lapply(1:num_columns, function(i) model_H1[[paste0("y_rep[", i, "]")]])
      
      y_rep_data <- do.call(cbind.data.frame, columns_list)
      colnames(y_rep_data) <- paste0("yRep_", 1:num_columns)
      p()
      y_rep_data
      
      
    })
  })
  
  plan(sequential)
  
  print("Doing Parallel CBF Simulations")
  
  current_data_y_rep <- list()
  for (model in 1:nrow(df_obs_bf_pos)) {
    # Initialize current_data_y_rep[[model]] as a list so that it can hold elements for i
    current_data_y_rep[[model]] <- list()
    
    for (post_pred_data in 1:100) {
      row_data <- y_rep_list[[model]][post_pred_data + 3000,]
      column_data <- t(row_data)
      
      # Convert to a data frame
      column_y_rep <- as.data.frame(column_data)
      
      # Current dataset with y_rep
      current_data_y_rep[[model]][[post_pred_data]] <- cbind.data.frame(column_y_rep[,1], current_data[,c(-1)])
      colnames(current_data_y_rep[[model]][[post_pred_data]]) <- colnames(historical_data)
    }
  }
  
  
  # Formula considering the posterior predictive samples
  bf_list <- list()
  
  # Define number of iterations
  n_iterations <- nrow(df_obs_bf_pos)
  
  # Register the parallel backend
  plan(multisession, workers = cores)
  
  # Run the loop in parallel
  opts <- list(packages = c("hdbayes"),
               seed = TRUE)
  #options(future.globals.maxSize= 891289600)
  
  with_progress({
    p <- progressor(along = 1:n_iterations)
    system.time( bf_list <- foreach(model = 1:n_iterations,.options.future = opts) %dofuture% {
      
      
      bayesF <- list()
      
      for (post_pred_data in 1:100) {
        
        all_data_rep <- list(current_data_y_rep[[model]][[post_pred_data]], historical_data)
        
        model_H1_rep <- glm.npp(
          formula = formula, family = family, data.list = all_data_rep,
          a0.lognc = a0.lognc$a0,
          lognc = matrix(a0.lognc$lognc, ncol = 1),
          a0.shape1 = df_obs_bf_pos$eta_H1[model], a0.shape2 = df_obs_bf_pos$nu_H1[model], 
          iter_warmup = 1000, iter_sampling = 1000,
          chains = 4, parallel_chains = 1,
          refresh = 0, seed = 433
        )
        
        log_ML_H1_rep <- glm.logml.npp(model_H1_rep)
        
        
        
        bf_rep <- log_ML_H1_rep$logml - log_ML_H0$logml
        
        bayesF[[post_pred_data]] <- bf_rep
      }
      
      p()
      
      bayesF
      
    })
  })
  
  plan(sequential)
  
  
  
  
  bf_data_rep <- as.data.frame(do.call(rbind, lapply(bf_list, unlist)))
  
  colnames(bf_data_rep) <- paste0("BF_rep", 1:100)
  
  print(bf_data_rep)
  
  bf_data_long <- bf_data_rep %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")
  
  df_obs_bf <- filter(df_obs_bf, obs_bf>=0)
  
  # Optimal CBF prior
  optimal_model_index <- calculate_maximized_value(bf_data_long$Value, df_obs_bf$obs_bf, 0.75)
  
  
  
  sim_results <- list(
    obs_bf_data = df_obs_bf,
    opt_index = optimal_model_index,
    opt_bf =  df_obs_bf[optimal_model_index,])
  
  results_list[[m]] <- sim_results
  
}

#save(results_list, file = "R/Simulations/GLM_Poisson/RData/GLM_Poisson_Sim.RData")



#   ____________________________________________________________________________
#   Running NPP CBF optimal models                                          ####

cbf_models <- list()

n_iterations <- length(results_list)

# Register the parallel backend
plan(multisession, workers = cores)

# Run the loop in parallel

opts <- list(packages = c("hdbayes"),
             seed = TRUE)
#options(future.globals.maxSize= 891289600)

with_progress({
  p <- progressor(along = 1:n_iterations)
  system.time( cbf_models <- foreach(model = 1:n_iterations,.options.future = opts) %dofuture% {
    
    
    current_data <- current_datasets[[model]]
    all_data <- list(current_data, historical_data)
    
      model_cbf <- glm.npp(
        formula = formula, family = family, data.list = all_data,
        a0.lognc = a0.lognc$a0,
        lognc = matrix(a0.lognc$lognc, ncol = 1),
        a0.shape1 = results_list[[model]]$opt_bf$eta_H1,
        a0.shape2 = results_list[[model]]$opt_bf$nu_H1, 
        iter_warmup = 1000, iter_sampling = 1000,
        chains = 4, parallel_chains = 1,
        refresh = 0, seed = 433
      )
      
    
    p()
    
    model_cbf
    
  })
})

plan(sequential)


#   ____________________________________________________________________________
#   Running NPP uniform models                                              ####

unif_models <- list()

n_iterations <- length(results_list)

# Register the parallel backend
plan(multisession, workers = cores)

# Run the loop in parallel

opts <- list(packages = c("hdbayes"),
             seed = TRUE)
#options(future.globals.maxSize= 891289600)

with_progress({
  p <- progressor(along = 1:n_iterations)
  system.time( unif_models <- foreach(model = 1:n_iterations,.options.future = opts) %dofuture% {
    
    
    current_data <- current_datasets[[model]]
    all_data <- list(current_data, historical_data)
    
    model_unif <- glm.npp(
      formula = formula, family = family, data.list = all_data,
      a0.lognc = a0.lognc$a0,
      lognc = matrix(a0.lognc$lognc, ncol = 1),
      a0.shape1 = 1,
      a0.shape2 = 1, 
      iter_warmup = 1000, iter_sampling = 1000,
      chains = 4, parallel_chains = 1,
      refresh = 0, seed = 433
    )
    
    
    p()
    
    model_unif
    
  })
})

plan(sequential)

#   ____________________________________________________________________________
#   Running RMAP models                                                     ####

rmap_models <- list()

n_iterations <- length(results_list)

# Set up parallel backend with number of cores to use
cores <- detectCores()

# Register the parallel backend
plan(multisession, workers = 13)

# Run the loop in parallel

opts <- list(packages = c("hdbayes"),
             seed = TRUE)
#options(future.globals.maxSize= 891289600)

with_progress({
  p <- progressor(along = 1:n_iterations)
  system.time( rmap_models <- foreach(model = 1:n_iterations,.options.future = opts) %dofuture% {
    
    
    current_data <- current_datasets[[model]]
    all_data <- list(current_data, historical_data)
    
    model_rmap <- glm.rmap(
      formula = formula, family = family, data.list = all_data,
      w = 0.5,
      iter_warmup = 1000, iter_sampling = 1000,
      chains = 4, parallel_chains = 1,
      refresh = 0, seed = 433)
    
    
    p()
    
    model_rmap
    
  })
})

plan(sequential)


#   ____________________________________________________________________________
#   Running SAM prior models                                                ####

# SAM prior function for Poisson data
source("R/Simulations/GLM_Poisson/Scripts/SAM_function.R", echo=TRUE)


alpha0 <- 1    
beta0  <- 1    

alpha_post <- alpha0 + sum(historical_data$y)       
beta_post  <- beta0 + nrow(historical_data)             


prior_historical <- mixgamma(c(1, alpha_post, beta_post),
                               param = 'ab')


##  ............................................................................
##  SAM Models                                                              ####

samp_models <- list()

n_iterations <- length(results_list)

# Register the parallel backend
plan(multisession, workers = cores)

# Run the loop in parallel
opts <- list(packages = c("hdbayes"),
             seed = TRUE)
#options(future.globals.maxSize= 891289600)

with_progress({
  p <- progressor(along = 1:n_iterations)
  system.time( samp_models <- foreach(model = 1:n_iterations,.options.future = opts) %dofuture% {
    
    
    current_data <- current_datasets[[model]]
    all_data <- list(current_data, historical_data)
    
    
    wSAM <- SAM_weight_pois(if.prior = prior_historical,
                            delta = 0.2, 
                            data = current_data$y)
    
    
    model_samp <- glm.rmap(
      formula = formula, family = family, data.list = all_data,
      w = wSAM,
      iter_warmup = 1000, iter_sampling = 1000,
      chains = 4, parallel_chains = 1,
      refresh = 0, seed = 433)
    
    
    p()
    
    model_samp
    
  })
})

plan(sequential)


save(cbf_models, unif_models,
     rmap_models, samp_models, file = "R/Simulations/GLM_Poisson/RData/GLM_Poisson_Sim_Models.RData")
