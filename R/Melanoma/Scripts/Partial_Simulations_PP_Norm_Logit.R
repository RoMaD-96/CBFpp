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


load("Logit_VM.RData")


stan_model_rep <- list("STAN/Logit_Sim_Norm_PP_Rep.stan")


code_block <- function(portion_rep, iteration) {
  
  print("Doing Parallel Simulations from 49")
  
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
      Ks <- c(5000)
      
      J <- 20
      
      ### GAM Approximation ###
      constant_data <- read.csv("delta_estimates_logit.csv")
      fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)
      
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
  plan(sequential)
  
  filename <- sprintf("bf_list_%d.RData", iteration)
  save(bf_list, file = filename)
  
  gc()
  return(bf_list)
}


iteration_by_block <- lapply(seq(1, 100, by = 2), function(x) c(x:(x + 1)))


print("Starting Parallel Simulations by Block")

BF_rep_dist <- lapply(49:length(iteration_by_block), function(iteration)  code_block(iteration_by_block, iteration))


#save.image(file = "Logit_2.RData")