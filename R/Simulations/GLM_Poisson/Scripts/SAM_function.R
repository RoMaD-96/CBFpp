SAM_weight_pois <- function(if.prior,
                            theta.h,
                            method.w,
                            prior.odds,
                            data,
                            delta,
                            y,     # total events, if no individual data
                            ...)
{
  
  if (!missing(data)) {
    if (is.vector(data) && is.numeric(data)) {
      # A vector of Poisson counts
      individual_counts <- data
    } else {
      stop("`data` format not recognized. Provide a numeric vector of counts.")
    }
  } else {
    # If data is missing, user must supply y directly
    if (missing(y)) {
      stop("Must provide either `data` (a vector of Poisson counts) or `y` (total events).")
    }
    individual_counts <- NA  # indicates we're using summary y
  }
  
  # Validate method.w
  if(!missing(method.w)){
    assertthat::assert_that(method.w %in% c('LRT', 'PPR'))
  }
  if (missing(method.w)) {
    message("Using the LRT (Likelihood Ratio Test) as the default method to calculate mixture weight for SAM priors.")
    method.w <- 'LRT'
  }
  
  # Validate prior.odds if method.w == 'PPR'
  if(method.w == 'PPR' & missing(prior.odds)){
    message("Missing the prior odds, set as 1.")
    prior.odds <- 1
  }
  
  # If theta.h is missing, use posterior mean from the informative prior
  if (missing(theta.h)) {
    message("Using the posterior mean from the informative prior as the estimate of lambda based on historical data.")
    theta.h <- summary(if.prior)['mean']
  }
  lambda_h_hat <- theta.h
  
  # Compute the ratio R
  
  logLikPois <- function(counts, lambda) {
    # sum of log(dpois(counts[i], lambda))
    sum(dpois(counts, lambda = lambda, log = TRUE))
  }
  
  lambda_plus  <- pmax(lambda_h_hat + delta, 1e-8)
  lambda_minus <- pmax(lambda_h_hat - delta, 1e-8)
  
  if (!all(is.na(individual_counts))) {
    # We have a vector of counts
    log_numer <- max(
      logLikPois(individual_counts, lambda_plus),
      logLikPois(individual_counts, lambda_minus)
    )
    log_denom <- logLikPois(individual_counts, lambda_h_hat)
  } else {
    # Treat y as a single Poisson observation with mean = lambda
    log_numer <- max(
      dpois(y, lambda = lambda_plus, log = TRUE),
      dpois(y, lambda = lambda_minus, log = TRUE)
    )
    log_denom <- dpois(y, lambda = lambda_h_hat, log = TRUE)
  }
  
  R <- exp(log_numer - log_denom)
  
  if (method.w == 'PPR') {
    R <- R / prior.odds
  }
  
  # Compute the SAM weight
  w <- 1 / (1 + R)
  
  return(w)
}