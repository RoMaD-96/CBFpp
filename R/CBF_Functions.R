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
    portion_cdf_above_zero <- 1 - cdf_values(0)

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

