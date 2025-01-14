#' get Stan data for PP

get.stan.data_pp = function(
    formula,
    family,
    data.list,
    a0_vals,
    offset.list       = NULL,
    beta_mean         = NULL,
    beta_sd           = NULL,
    disp_mean         = NULL,
    disp_sd           = NULL
) {
  data_checks(formula, family, data.list, offset.list)
  
  res          = stack_data(formula = formula,
                            data.list = data.list)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  p            = ncol(X)
  N            = length(y)
  K            = length(end.index)
  fam.indx     = get.dist_link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]
  
  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }
  
  ## Default prior on regression coefficients is N(0, 10^2)
  if ( !is.null(beta_mean) ){
    if ( !( is.vector(beta_mean) & (length(beta_mean) %in% c(1, p)) ) )
      stop("beta_mean must be a scalar or a vector of length ", p, " if beta_mean is not NULL")
  }
  beta_mean = to_vector(param = beta_mean, default.value = 0, len = p)
  if ( !is.null(beta_sd) ){
    if ( !( is.vector(beta_sd) & (length(beta_sd) %in% c(1, p)) ) )
      stop("beta_sd must be a scalar or a vector of length ", p, " if beta_sd is not NULL")
  }
  beta_sd = to_vector(param = beta_sd, default.value = 10, len = p)
  
  ## check a0 values
  if ( !( is.vector(a0_vals) & (length(a0_vals) %in% c(1, K-1)) ) )
    stop("a0_vals must be a scalar or a vector of length ", K-1)
  a0_vals = to_vector(param = a0_vals, len = K-1)
  if ( any(a0_vals < 0 | a0_vals > 1 ) )
    stop("Each element of a0_vals must be a scalar between 0 and 1")
  a0_vals = c(1, a0_vals) # first element = 1 for current data
  
  ## Default half-normal prior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp_mean) ){
    if ( !( is.vector(disp_mean) & (length(disp_mean) == 1) ) )
      stop("disp_mean must be a scalar if disp_mean is not NULL")
  }
  disp_mean = to_vector(param = disp_mean, default.value = 0, len = 1)
  if ( !is.null(disp_sd) ){
    if ( !( is.vector(disp_sd) & (length(disp_sd) == 1) ) )
      stop("disp_sd must be a scalar if disp_sd is not NULL")
  }
  disp_sd = to_vector(param = disp_sd, default.value = 10, len = 1)
  
  standata = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'mean_beta'       = beta_mean,
    'sd_beta'         = beta_sd,
    'a0_vals'         = a0_vals,
    'disp_mean'       = disp_mean,
    'disp_sd'         = disp_sd,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )
  return(standata)
}



#' get Stan data for NPP

get.stan.data_npp = function(
    formula,
    family,
    data.list,
    a0.lognc,
    lognc,
    offset.list       = NULL,
    beta_mean         = NULL,
    beta_sd           = NULL,
    disp_mean         = NULL,
    disp_sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL
) {
  data_checks(formula, family, data.list, offset.list)
  
  res          = stack_data(formula = formula,
                            data.list = data.list)
  y            = res$y
  X            = res$X
  start.index  = res$start.index
  end.index    = res$end.index
  p            = ncol(X)
  N            = length(y)
  K            = length(end.index)
  fam.indx     = get.dist_link(family)
  dist         = fam.indx[1]
  link         = fam.indx[2]
  
  ## Default offset for each data set is a vector of 0s
  if ( is.null(offset.list) ){
    offset = rep(0, N)
  }else {
    offset = unlist(offset.list)
  }
  
  ## Default prior on regression coefficients is N(0, 10^2)
  if ( !is.null(beta_mean) ){
    if ( !( is.vector(beta_mean) & (length(beta_mean) %in% c(1, p)) ) )
      stop("beta_mean must be a scalar or a vector of length ", p, " if beta_mean is not NULL")
  }
  beta_mean = to_vector(param = beta_mean, default.value = 0, len = p)
  if ( !is.null(beta_sd) ){
    if ( !( is.vector(beta_sd) & (length(beta_sd) %in% c(1, p)) ) )
      stop("beta_sd must be a scalar or a vector of length ", p, " if beta_sd is not NULL")
  }
  beta_sd = to_vector(param = beta_sd, default.value = 10, len = p)
  
  ## check a0.lognc and lognc
  if ( length(a0.lognc) != nrow(lognc) )
    stop('the number of rows in lognc must be the same as the length of a0.lognc')
  if ( ncol(lognc) != (K - 1) )
    stop('the number of columns in lognc must be the same as the number of historical data sets')
  if ( any(is.na(a0.lognc) ) )
    stop('a0.lognc must not have missing values')
  if ( any(is.na(lognc)) )
    stop('lognc must not have missing values')
  if ( any(a0.lognc < 0) || any(a0.lognc > 1) )
    stop('each element of a0.lognc should be between 0 and 1')
  
  ## Default half-normal prior on dispersion parameters (if exist) is N^{+}(0, 10^2)
  if ( !is.null(disp_mean) ){
    if ( !( is.vector(disp_mean) & (length(disp_mean) == 1) ) )
      stop("disp_mean must be a scalar if disp_mean is not NULL")
  }
  disp_mean = to_vector(param = disp_mean, default.value = 0, len = 1)
  if ( !is.null(disp_sd) ){
    if ( !( is.vector(disp_sd) & (length(disp_sd) == 1) ) )
      stop("disp_sd must be a scalar if disp_sd is not NULL")
  }
  disp_sd = to_vector(param = disp_sd, default.value = 10, len = 1)
  
  ## Default lower bound for each a0 is 0; default upper bound for each a0 is 1
  if ( !is.null(a0.lower) ){
    if ( !( is.vector(a0.lower) & (length(a0.lower) %in% c(1, K - 1)) ) )
      stop("a0.lower must be a scalar or a vector of length ", K - 1, " if a0.lower is not NULL")
  }
  a0.lower = to_vector(param = a0.lower, default.value = 0, len = K - 1)
  if ( !is.null(a0.upper) ){
    if ( !( is.vector(a0.upper) & (length(a0.upper) %in% c(1, K - 1)) ) )
      stop("a0.upper must be a scalar or a vector of length ", K - 1, " if a0.upper is not NULL")
  }
  a0.upper = to_vector(param = a0.upper, default.value = 1, len = K - 1)
  
  standata = list(
    'K'               = K,
    'N'               = N,
    'start_idx'       = start.index,
    'end_idx'         = end.index,
    'p'               = p,
    'y'               = y,
    'X'               = X,
    'mean_beta'       = beta_mean,
    'sd_beta'         = beta_sd,
    'disp_mean'       = disp_mean,
    'disp_sd'         = disp_sd,
    's'               = length(a0.lognc),
    'a0_lognc'        = a0.lognc,
    'lognc'           = lognc,
    'a0_shape1'       = a0.shape1,
    'a0_shape2'       = a0.shape2,
    'a0_lower'        = a0.lower,
    'a0_upper'        = a0.upper,
    'dist'            = dist,
    'link'            = link,
    'offs'            = offset
  )
  return(standata)
}


#' get distribution and link for a GLM

get.dist_link = function(family) {
  fams  = c('binomial', 'poisson', 'gaussian', 'Gamma', 'inverse.gaussian')
  links = c(
    'identity', 'log', 'logit', 'inverse', 'probit', 'cauchit', 'cloglog'
    ,'sqrt', '1/mu^2'
  )
  fam.id  = which(family$family == fams)
  link.id = which(family$link == links)
  c(fam.id, link.id)
}

#' get Stan data for normal/half-normal prior

get.stan.data_post = function(
    formula,
    family,
    data.list,
    offset.list       = NULL,
    beta_mean         = NULL,
    beta_sd           = NULL,
    disp_mean         = NULL,
    disp_sd           = NULL
) {
  if( length(data.list) > 1 ){
    data.list   = list(data.list[[1]])
  }
  if ( length(offset.list) > 1 ) {
    offset.list = list(offset.list[[1]])
  }
  if ( length(disp_mean) > 1 ) {
    disp_mean = disp_mean[1]
  }
  if ( length(disp_sd) > 1 ) {
    disp_sd = disp_sd[1]
  }
  
  standata_pp = get.stan.data_pp(
    formula     = formula,
    family      = family,
    data.list   = data.list,
    a0_vals     = 0,
    offset.list = offset.list,
    beta_mean   = beta_mean,
    beta_sd     = beta_sd,
    disp_mean   = disp_mean,
    disp_sd     = disp_sd
  )
  
  standata = list(
    'n'               = standata_pp$N,
    'p'               = standata_pp$p,
    'y'               = standata_pp$y,
    'X'               = standata_pp$X,
    'mean_beta'       = standata_pp$mean_beta,
    'sd_beta'         = standata_pp$sd_beta,
    'disp_mean'       = standata_pp$disp_mean,
    'disp_sd'         = standata_pp$disp_sd,
    'dist'            = standata_pp$dist,
    'link'            = standata_pp$link,
    'offs'            = standata_pp$offs
  )
  return(standata)
}




#' check if the input data are in appropriate forms for all methods except for LEAP

data_checks = function(
    formula, family, data.list, offset.list
) {
  if ( !inherits(formula, "formula") )
    stop('formula must be of type "formula"')
  if ( !inherits(family, 'family') )
    stop('family must be of type "family" (e.g., cannot be a character--use binomial() instead of "binomial"). See help(family)')
  if ( !formula.tools::is.two.sided(formula) )
    stop('formula must be two-sided')
  yname = formula.tools::lhs.vars(formula)
  if ( length(yname) != 1 )
    stop('formula must contain exactly 1 lhs variable name')
  varnames = all.vars(formula)
  if ( !( is.list(data.list) ) )
    stop("data.list must be a list of data.frames")
  for( i in seq_len( length(data.list) ) ){
    if ( !( is.data.frame(data.list[[i]]) ) )
      stop("element ", i, " in data.list must be a data.frame")
    if ( any( is.na(data.list[[i]]) ) )
      stop("element ", i, " in data.list cannot contain missing values")
    if ( !( all( varnames %in% names(data.list[[i]]) ) ) )
      stop("formula contains terms not in element ", i, " in data.list")
  }
  if ( !( is.null(offset.list) ) ){
    if ( !( is.list(offset.list) ) )
      stop("offset.list must be a list of vectors if offset.list is not NULL")
    if ( length(offset.list) != length(data.list) )
      stop("offset.list and data.list must have equal lengths if offset.list is not NULL")
    for( i in seq_len( length(offset.list) ) ){
      if ( !( is.vector(offset.list[[i]]) ) )
        stop("element ", i, " in offset.list must be a vector")
      if ( any( is.na(offset.list[[i]]) ) )
        stop("element ", i, " in offset.list cannot contain missing values")
      if ( length(offset.list[[i]]) != nrow(data.list[[i]]) )
        stop("the length of element ", i, " in offset.list must be equal to the number of rows in element ", i, " in data.list if offset.list is not NULL")
    }
  }
}



#' reshape the input data list into a list which contains the response vector y (from all data sets), the covariate matrix X (from all data sets),

stack_data = function(
    formula, data.list
) {
  ## get stacked design matrix and response vector using formula
  X = lapply(data.list, function(s){
    stats::model.matrix(formula, s)
  })
  X = do.call(rbind, X)
  y = lapply(data.list, function(s){
    s[, all.vars(formula)[1]]
  })
  y = unlist(y)
  
  ## reset indices of X
  rownames(X) = NULL
  
  # get starting and ending indices for each data
  num.obs = sapply(data.list, function(s){
    nrow(s)
  })
  end.index   = cumsum(num.obs)
  start.index = c(1, end.index[-length(data.list)] + 1)
  
  return(list(X = X, y = y, start.index = start.index, end.index = end.index))
}


#' transfer a scalar/vector/NULL into a vector of given length and (default) values

to_vector = function(
    param, default.value = 0, len = 1
) {
  if ( is.null(param) ){
    param = rep(default.value, len)
  }else if ( length(param) == 1 ){
    param = rep(as.numeric(param), len)
  }else {
    param = as.numeric(param)
  }
  return(param)
}



rename_params = function(
    fit, oldnames, newnames
) {
  pars = fit$metadata()$model_params
  pars = c(pars[1], oldnames, (pars[!pars %in% oldnames])[-1])
  d    = fit$draws(format = 'draws_df', variables = pars)
  posterior::variables(d)[posterior::variables(d) %in% oldnames] = newnames
  return(d)
}

#   ____________________________________________________________________________
#   Posterior of normalized power prior                                     ####


glm.npp_post_pred = function(
    formula,
    family,
    data.list,
    a0.lognc,
    lognc,
    offset.list       = NULL,
    beta_mean         = NULL,
    beta_sd           = NULL,
    disp_mean         = NULL,
    disp_sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for NPP
  standata = get.stan.data_npp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    a0.lognc       = a0.lognc,
    lognc          = lognc,
    offset.list    = offset.list,
    beta_mean      = beta_mean,
    beta_sd        = beta_sd,
    disp_mean      = disp_mean,
    disp_sd        = disp_sd,
    a0.shape1      = a0.shape1,
    a0.shape2      = a0.shape2,
    a0.lower       = a0.lower,
    a0.upper       = a0.upper
  )
  

  stan_file <- "R/Melanoma/CmdStan/NPP_Logit_Post_Pred.stan"
  glm_npp_posterior_post_pred <- cmdstan_model(stan_file)
  
  ## fit model in cmdstanr
  fit = glm_npp_posterior_post_pred$sample(data = standata,
                                    iter_warmup = iter_warmup,
                                    iter_sampling = iter_sampling,
                                    chains = chains,...)
  
  ## rename parameters
  p        = standata$p
  X        = standata$X
  K        = standata$K
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('a0s[', 1:(K-1), ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:(K-1)))
  d        = rename_params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standata
  return(d)
}




#   ____________________________________________________________________________
#   NPP mixture prior                                                       ####


glm.npp_mix = function(
    formula,
    family,
    data.list,
    a0.lognc,
    lognc,
    offset.list       = NULL,
    beta_mean         = NULL,
    beta_sd           = NULL,
    disp_mean         = NULL,
    disp_sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for NPP
  standata = get.stan.data_npp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    a0.lognc       = a0.lognc,
    lognc          = lognc,
    offset.list    = offset.list,
    beta_mean      = beta_mean,
    beta_sd        = beta_sd,
    disp_mean      = disp_mean,
    disp_sd        = disp_sd,
    a0.shape1      = a0.shape1,
    a0.shape2      = a0.shape2,
    a0.lower       = a0.lower,
    a0.upper       = a0.upper
  )
  
  
  stan_file <- "R/Melanoma/CmdStan/NPP_Logit_Mix.stan"
  glm_npp_post_mix <- cmdstan_model(stan_file)
  
  ## fit model in cmdstanr
  fit = glm_npp_post_mix$sample(data = standata,
                                iter_warmup = iter_warmup,
                                iter_sampling = iter_sampling,
                                chains = chains,...)
  
  ## rename parameters
  p        = standata$p
  X        = standata$X
  K        = standata$K
  oldnames = paste0("beta[", 1:p, "]")
  newnames = colnames(X)
  
  if ( !family$family %in% c('binomial', 'poisson') ) {
    oldnames = c(oldnames, 'dispersion[1]')
    newnames = c(newnames, 'dispersion')
  }
  oldnames = c(oldnames, paste0('a0s[', 1:(K-1), ']'))
  newnames = c(newnames, paste0('a0_hist_', 1:(K-1)))
  d        = rename_params(fit = fit, oldnames = oldnames, newnames = newnames)
  ## add data used for the variables specified in the data block of the Stan program as an attribute
  attr(x = d, which = 'data') = standata
  return(d)
}

#   ____________________________________________________________________________
#   Log-Likelihood models                                                   ####

glm.npp_post_loo = function(
    formula,
    family,
    data.list,
    a0.lognc,
    lognc,
    offset.list       = NULL,
    beta_mean         = NULL,
    beta_sd           = NULL,
    disp_mean         = NULL,
    disp_sd           = NULL,
    a0.shape1         = 1,
    a0.shape2         = 1,
    a0.lower          = NULL,
    a0.upper          = NULL,
    iter_warmup       = 1000,
    iter_sampling     = 1000,
    chains            = 4,
    ...
) {
  ## get Stan data for NPP
  standata = get.stan.data_npp(
    formula        = formula,
    family         = family,
    data.list      = data.list,
    a0.lognc       = a0.lognc,
    lognc          = lognc,
    offset.list    = offset.list,
    beta_mean      = beta_mean,
    beta_sd        = beta_sd,
    disp_mean      = disp_mean,
    disp_sd        = disp_sd,
    a0.shape1      = a0.shape1,
    a0.shape2      = a0.shape2,
    a0.lower       = a0.lower,
    a0.upper       = a0.upper
  )
  
  
  stan_file <- "R/Melanoma/CmdStan/NPP_Logit_Log_Lik.stan"
  glm_npp_posterior_log_lik <- cmdstan_model(stan_file)
  
  ## fit model in cmdstanr
  fit = glm_npp_posterior_log_lik$sample(data = standata,
                                           iter_warmup = iter_warmup,
                                           iter_sampling = iter_sampling,
                                           chains = chains,...)
  
 
  return(fit)
}



