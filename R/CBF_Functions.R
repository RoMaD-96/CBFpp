#   ____________________________________________________________________________
#   CBF Optimal Criterion                                                   ####


calculate_maximized_value <- function(distribution_values, observed_BF, hpdi_prob = 0.95) {
  results <- numeric(length(observed_BF)) # Initialize a results vector
  
  # Calculate the number of distributions
  num_distributions <- length(distribution_values) / length(observed_BF)
  
  for (i in seq_along(observed_BF)) {
    # Select values for the current distribution
    current_distribution_values <- distribution_values[((i - 1) * num_distributions + 1):(i * num_distributions)]
    
    # Calculate the CDF for the current distribution values
    cdf_values <- ecdf(current_distribution_values)
    
    # Calculate the portion of the CDF greater than zero
    portion_cdf_above_zero <-  1-cdf_values(0)
    
    # Calculate the HPDI
    hpdi <- rethinking::HPDI(current_distribution_values, prob = hpdi_prob)
    
    # Check conditions and calculate the result
    if (portion_cdf_above_zero >= 0.5 && observed_BF[i] >= hpdi[1] && observed_BF[i] <= hpdi[2]) {
      results[i] <- portion_cdf_above_zero * observed_BF[i]
    } else {
      results[i] <- 0
    }
  }
  
  # Check if all results are zero
  if (max(results) == 0) {
    print("The model under the Null Hypothesis is optimal")
  } else {
    # Return the index of the maximum value in results
    print(which.max(results))
  }
}


#   ____________________________________________________________________________
#   STAN Functions                                                          ####


##  ............................................................................
##  Fixed Power Prior Posterior                                             ####


fixed_posterior <- function(stan_setup_list,
                            stan_file_list,
                            iterations = 2000,
                            treedepth_max = 15,
                            adaptive_delta = 0.95) {
  fixed_pp_stan <- list()
  for (i in 1:length(stan_setup_list)) {
    cat("Doing delta =", stan_setup_list[[i]]$delta, "\ value")
    fixed_pp_stan[[i]] <- stan(
      file = unlist(stan_file_list),
      iter = iterations,
      data = stan_setup_list[[i]],
      control = list(adapt_delta = adaptive_delta,
                     max_treedepth = treedepth_max)
    )
    cat("\n")
  }
  # List Names
  if (is.null(names(stan_setup_list)) == FALSE) {
    names(fixed_pp_stan) <- names(stan_setup_list)
  }
  else {
    print("Please assign names in the list")
  }
  
  return(fixed_pp_stan)
}


##  ............................................................................
##  Unnormalized Power Prior Posterior                                      ####

unnormalized_posterior_sampl <- function(stan_setup_list,
                                         stan_file_list,
                                         iterations = 2000,
                                         treedepth_max = 15,
                                         adaptive_delta = 0.95,
                                         chains = 4, 
                                         ref_rate = max(iterations/10, 1),
                                         info_model = TRUE) {
  random_pp_stan <- list()
  stan_seed <- 4231
  model <- stan_model(unlist(stan_file_list))
  for (i in 1:length(stan_setup_list)) {
    if (info_model == TRUE) {
      cat("Doing Beta with (eta,nu) =",
          c(stan_setup_list[[i]]$eta, stan_setup_list[[i]]$nu),
          "\n , considering seed =", c(stan_seed), "\n")
    }
    random_pp_stan[[i]] <-
      sampling(
        model,
        data = stan_setup_list[[i]],
        iter = iterations,
        chains = chains,
        seed = stan_seed,
        refresh = ref_rate,
        control = list(adapt_delta = adaptive_delta,
                       max_treedepth = treedepth_max)
      )
    if (info_model == TRUE) {
      cat("\n")
    }
  }
  
  # List Names
  if (is.null(names(stan_setup_list)) == FALSE) {
    names(random_pp_stan) <- names(stan_setup_list)
  }
  else {
    print("Please assign names in the stan_setup_list")
  }
  
  return(random_pp_stan)
  
}

### Using stan() function

unnormalized_posterior <- function(stan_setup_list,
                                   stan_file_list,
                                   iterations = 2000,
                                   treedepth_max = 15,
                                   adaptive_delta = 0.95) {
  random_pp_stan <- list()
  for (i in 1:length(stan_setup_list)) {
    cat("Doing Beta with (eta,nu) =",
        c(stan_setup_list[[i]]$eta, stan_setup_list[[i]]$nu),
        "\n")
    random_pp_stan[[i]] <- stan(
      file = unlist(stan_file_list),
      data = stan_setup_list[[i]],
      iter = iterations,
      control = list(adapt_delta = adaptive_delta,
                     max_treedepth = treedepth_max)
    )
    cat("\n")
  }
  
  # List Names
  if (is.null(names(stan_setup_list)) == FALSE) {
    names(random_pp_stan) <- names(stan_setup_list)
  }
  else {
    print("Please assign names in the stan_setup_list")
  }
  
  return(random_pp_stan)
  
}



##  ............................................................................
##  Approximate Normalized                                                  ####


### Function with a unique K ####

norm_post_pp <- function(K,
                         stan_setup_list,
                         stan_file_list,
                         iterations = 2000,
                         treedepth_max = 15,
                         adaptive_delta = 0.95,
                         chains = 4, 
                         ref_rate = max(iterations/10, 1)) {
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit_gam, newdata = data.frame(a0 = pred_a0s)))
  norm_pp_stan <- list()
  for (i in 1:length(stan_setup_list)) {
    cat("Doing Beta with (eta,nu) =",
        c(stan_setup_list[[i]]$eta, stan_setup_list[[i]]$nu),
        "\n")
    stan_setup_list[[i]]$pred_grid_x <- a0_grid$a0
    stan_setup_list[[i]]$pred_grid_y <- a0_grid$lc_pred
    model <- stan_model(unlist(stan_file_list))
    norm_pp_stan[[i]] <-  sampling(
      model,
      data = stan_setup_list[[i]],
      iter = iterations,
      chains = chains,
      refresh = ref_rate,
      control = list(adapt_delta = adaptive_delta,
                     max_treedepth = treedepth_max))
    cat("\n")
  }
  # List Names
  if (is.null(names(stan_setup_list)) == FALSE) {
    names(norm_pp_stan) <- names(stan_setup_list)
  } 
  else {
    print("Please assign names in the list")
  }
  
  return(norm_pp_stan)
  
}



### Function with multiple K's ####

norm_posterior_K <- function( K,
                              stan_setup_list,
                              stan_file_list,
                              iterations = 2000,
                              treedepth_max = 15,
                              adaptive_delta = 0.95) {
  cat("Doing Beta with (eta,nu) =",
      c(stan_setup_list$eta, stan_setup_list$nu),
      "\n")
  cat("Doing k=", K, "\n")
  pred_a0s <- seq(0, max(constant_data$a0), length.out = K)
  a0_grid <- data.frame(a0 = pred_a0s,
                        lc_pred = predict(fit_gam, newdata = data.frame(a0 = pred_a0s)))
  stan_setup_list$pred_grid_x <- a0_grid$a0
  stan_setup_list$pred_grid_y <- a0_grid$lc_pred
  approx_norm_posterior <-
    stan(
      file = unlist(stan_file_list),
      data = stan_setup_list,
      iter = iterations,
      control = list(adapt_delta = adaptive_delta,
                     max_treedepth = treedepth_max),
      refresh = 0
    )
  return(approx_norm_posterior)
}



app_norm_post <- function(K,
                          stan_setup_list,
                          stan_file_list,
                          iterations = 2000,
                          treedepth_max = 15) {
  norm_post_value <- list()
  for (i in 1:length(stan_setup_list)) {
    norm_post_value[[i]] <- lapply(K, stan_setup_list = stan_setup_list[[i]], stan_file_list = stan_file_list, norm_posterior_K)
    names(norm_post_value[[i]]) <- paste(K)
  }
  
  # List Names
  if (is.null(names(stan_setup_list)) == FALSE) {
    names(norm_post_value) <- names(stan_setup_list)
  }
  else {
    print("Please assign names in the list")
  }
  return(norm_post_value)
}
