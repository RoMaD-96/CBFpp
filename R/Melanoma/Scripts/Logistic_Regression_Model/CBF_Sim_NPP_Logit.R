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
  "ggplot2",
  "ggridges",
  "tidyr"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#   ____________________________________________________________________________
#   Grid Search                                                             ####

eta <- round(seq(0.5, 6, by = 0.5), 2)
nu <- round(seq(0.5, 6, by = 0.5), 2)
grid <- expand.grid(x = eta, y = nu)
grid <- as.data.frame(grid)
model_id <- seq(1, nrow(grid))
grid <- cbind(grid, model_id)
model_id <- model_id[-14]
grid <- grid[-14, ]
combinations <- expand.grid(model_id, 14)

#   ____________________________________________________________________________
#   Datasets                                                                ####

set.seed(442)

data("E2696")
data("E1694")
historical_data <- E2696
current_data <- E1694
all_data <- list(current_data, historical_data)

formula <- failind ~ age + treatment + sex + perform
p <- length(attr(terms(formula), "term.labels"))
family <- binomial("logit")

#   ____________________________________________________________________________
#   Approximate normalizing constant                                        ####

a0 <- seq(0, 1, length.out = 21)

if (requireNamespace("parallel")) {
  if (instantiate::stan_cmdstan_exists()) {
    logncfun <- function(a0, ...) {
      hdbayes::glm.npp.lognc(
        formula = formula,
        family = family,
        histdata = historical_data,
        a0 = a0,
        ...
      )
    }
    cl <- parallel::makeCluster(10)
    parallel::clusterSetRNGStream(cl, 123)
    parallel::clusterExport(cl, varlist = c("formula", "family", "historical_data"))
    time.npp.1 <- system.time(
      a0.lognc <- parLapply(
        cl = cl,
        X = a0,
        fun = logncfun,
        iter_warmup = 1000,
        iter_sampling = 1000,
        chains = 4,
        refresh = 0
      )
    )
    parallel::stopCluster(cl)
    a0.lognc <- data.frame(do.call(rbind, a0.lognc))
  }
}

#   ____________________________________________________________________________
#   Observed Bayes Factor                                                  ####

model_unif_H0 <- glm.npp(
  formula = formula,
  family = family,
  data.list = all_data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = 1,
  a0.shape2 = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  seed = 12345
)

set.seed(4321)
log_ML_H0 <- glm.logml.npp(model_unif_H0)

# Define number of iterations
n_it_obs_bf <- nrow(combinations)

# Register the parallel backend
cores <- detectCores() - 1

plan(multisession, workers = cores)

# Run the loop in parallel
opts <- list(packages = c("hdbayes"), seed = TRUE)

# Observed Bayes Factor
obs_bf_list <- list()

with_progress({
  p <- progressor(along = 1:n_it_obs_bf)
  system.time( 
    obs_bf_list <- foreach(i = 1:n_it_obs_bf,.options.future = opts) %dofuture% {

  # Observed BF

    model_H1 <- glm.npp(
      formula = formula, family = family, data.list = all_data,
      a0.lognc = a0.lognc$a0,
      lognc = matrix(a0.lognc$lognc, ncol = 1),
      a0.shape1 = grid$x[i], a0.shape2 = grid$y[i], 
      iter_warmup = 1000, iter_sampling = 1000,
      chains = 4, parallel_chains = 1,
      refresh = 100, seed = 12345
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

df_obs_bf <- cbind.data.frame(obs_bf, grid$x, grid$y)
colnames(df_obs_bf) <- c("obs_bf", "eta_H1", "nu_H1")

plan(sequential)

save(df_obs_bf, file = "R/Melanoma/RData/obs_bf_melanoma.RData")


#   ____________________________________________________________________________
#   Posterior predictive samples                                            ####

df_obs_bf_pos <- filter(df_obs_bf, obs_bf>=0)

n_iter <- nrow(df_obs_bf_pos)

# Register the parallel backend
plan(multisession, workers = cores)

# Run the loop in parallel
opts <- list(packages = c("hdbayes", "cmdstanr"), seed = TRUE)

y_rep_list <- list()

with_progress({
  p <- progressor(along = 1:n_iter)
  system.time( 
    y_rep_list <- foreach(i = 1:n_iter,.options.future = opts) %dofuture% {
    
    # Code with posterior predictive samples
    source("R/Melanoma/Scripts/Logistic_Regression_Model/Fun_Post_Pred.R", echo=TRUE)
    
    model_H1 <- glm.npp_post_pred(
      formula = formula, family = family, data.list = all_data,
      a0.lognc = a0.lognc$a0,
      lognc = matrix(a0.lognc$lognc, ncol = 1),
      a0.shape1 = df_obs_bf_pos$eta_H1[i], a0.shape2 = df_obs_bf_pos$nu_H1[i], 
      iter_warmup = 1000, iter_sampling = 1000,
      chains = 4, parallel_chains = 1,
      refresh = 100, seed = 12345
    )
    
    num_columns <- nrow(current_data)
    
    columns_list <- lapply(1:num_columns, function(i) model_H1[[paste0("y_rep[", i, "]")]])
    
    y_rep_data <- do.call(cbind.data.frame, columns_list)
    
    colnames(y_rep_data) <- paste0("yRep_", 1:num_columns)
    
    p()
    
    y_rep_data

    
  })
})

plan(sequential)

save(y_rep_list, file = "R/Melanoma/RData/obs_rep.RData")


##  ............................................................................
##  Replicated samples                                                      ####

current_data_y_rep <- list()
for (model in 1:nrow(df_obs_bf_pos)) {
  current_data_y_rep[[model]] <- list()
  
  for (post_pred_data in 1:100) {
    row_data <- y_rep_list[[model]][post_pred_data + 3000,]
    column_data <- t(row_data)
    column_y_rep <- as.data.frame(column_data)

    current_data_y_rep[[model]][[post_pred_data]] <- cbind.data.frame(column_y_rep[,1], current_data[,-2])
    current_data_y_rep[[model]][[post_pred_data]] <- current_data_y_rep[[model]][[post_pred_data]][,c(2,1,3:6)]
    colnames(current_data_y_rep[[model]][[post_pred_data]]) <- colnames(historical_data)
  }
}

##  ............................................................................
##  Replicated Bayes factor                                                 ####

bf_list <- list()

# Define number of iterations
n_iterations <- nrow(df_obs_bf_pos)

# Register the parallel backend
plan(multisession, workers = cores)
  
# Run the loop in parallel
opts <- list(packages = c("hdbayes"), seed = TRUE)

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
          refresh = 0, seed = 12345
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

save(bf_list, file = "R/Melanoma/RData/bf_list.RData")

