#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "ggthemes",
  "progressr",
  "ggplot2",
  "ggpubr",
  "colorspace",
  "modi",
  "overlapping")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#   ____________________________________________________________________________
#   Functions                                                               ####

# Marginal Likelihood Function
intFun <- function(delta) {
  dnorm(x = theta, mean = theta_0, sd = sqrt(se^2 + se_0^2/delta)) *
    dbeta(x = delta, shape1 = eta_null, shape = nu_null)
}


# Joint Posterior Function
joint_post <- function(delta_seq, theta_seq, theta, st_dev, theta_0, st_dev_0, eta, nu) {
  
  # Marginal Likelihood Function
  intFun <- function(delta) {
    dnorm(x = theta, mean = theta_0, sd = sqrt(se^2 + se_0^2/delta)) *
      dbeta(x = delta, shape1 = eta_null, shape = nu_null)
  }
  
  num_post <- dnorm(x = theta, mean = theta_seq, sd =  st_dev)*dnorm( x = theta_seq, mean = theta_0, sd = sqrt(st_dev_0^2/delta_seq))*
    dbeta(x = delta_seq, shape1 = eta, shape2 = nu)
  
  normconst <- integrate(f = intFun, lower = 0, upper = 1)$value
  joint_post <- num_post/normconst
  return(joint_post)
}


# Marginal Posterior of Theta 
theta_m_post <- function(theta_seq, theta, theta_0, st_dev, st_dev_0, eta, nu){
  
  # Marginal Likelihood Function
  intFun <- function(delta) {
    dnorm(x = theta, mean = theta_0, sd = sqrt(se^2 + se_0^2/delta)) *
      dbeta(x = delta, shape1 = eta_null, shape = nu_null)
  }
  
  normconst <- integrate(f = intFun, lower = 0, upper = 1)$value
  
  theta_post <- sapply(X = theta_seq, FUN = function(t) {
    integrate(f = function(delta) {
      (dnorm(x = theta, mean = t, sd = st_dev) *
         dnorm(x = t, mean = theta_0, sd = sqrt(st_dev_0^2/delta)) *
         dbeta(x = delta, shape1 = eta, shape2 = nu))/normconst
    }, lower = 0, upper = 1)$value
  })
  return(theta_post)
}

# Marginal Posterior of Delta
delta_m_post <- function(delta_seq, theta, theta_0, st_dev, st_dev_0, eta, nu){
  # Marginal Likelihood Function
  intFun <- function(delta) {
    dnorm(x = theta, mean = theta_0, sd = sqrt(se^2 + se_0^2/delta)) *
      dbeta(x = delta, shape1 = eta_null, shape = nu_null)
  }
  
  normconst <- integrate(f = intFun, lower = 0, upper = 1)$value
  num_delta <- dnorm(x = theta, mean = theta_0, sd = (st_dev+st_dev_0/sqrt(delta_seq)))*
    dbeta(x = delta_seq, shape1 = eta, shape2 = nu)
  delta_post <- num_delta/normconst
  
  return(delta_post)
}

# Marginal Likelihood Function
marg_lik <- function(delta_seq, eta, nu, theta, theta_0, st_dev, st_dev_0) {
  intFun_m <- function(delta_seq) {
    dnorm(x = theta, mean = theta_0, sd = sqrt(st_dev^2 + st_dev^2/delta_seq)) *
      dbeta(x = delta_seq, shape1 = eta, shape = nu)
  }  
  marg_lik_val <- integrate(f = intFun_m, lower = 0, upper = 1)$value
  return(marg_lik_val)
  
}

# Observed BF
bf_fun <- function(delta_seq, theta, theta_0, st_dev, st_dev_0, eta_null, nu_null, eta_alt, nu_alt, log_val = TRUE){
  
  # Marginal Likelihoods
  marg_lik_null <- marg_lik(delta_seq, eta_null, nu_null, theta, theta_0, st_dev, st_dev_0)
  marg_lik_alt <- marg_lik(delta_seq, eta_alt, nu_alt, theta, theta_0, st_dev, st_dev_0)
  
  # Log-value assessment 
  if(log_val == TRUE){
    bayes_factor <- log(marg_lik_alt/marg_lik_null, 10) 
  } 
  else{
    bayes_factor <- marg_lik_alt/marg_lik_null 
  }
  
  return(bayes_factor)
}


#   ____________________________________________________________________________
#   Starting Simulations                                                    ####

par_seq <-  seq(0,6, by = 0.2)
for (M in 1:length(par_seq)) {
  
  cat("Doing theta =",
      par_seq[M],
      "\n")
  
##  ............................................................................
##  Data                                                                    ####
  
set.seed(4231)
  
# Historical Data Parameters
true_theta_0 <- 0
true_se_0 <- 1
  
#y_0 <- rnorm(200, true_theta_0, true_se_0)
  
theta_0 <- true_theta_0
se_0 <- true_se_0

# Current Data Parameters
true_theta <- par_seq[M]
true_se <- 1

#y <-  rnorm(200, true_theta, true_se)

theta <- true_theta
se <- true_se
  
# overlap(list(y=rnorm(1000,true_theta,true_se), y_0=rnorm(1000,true_theta_0, true_se_0)))
  
##  ............................................................................
##  Model and Grids Specification                                           ####
  
# Beta Parameters
eta_null <- 1
nu_null <- 1
  
eta_alt <- seq(0.5, 6,by = 0.5)
nu_alt <- seq(0.5, 6,by = 0.5)
  
grid_beta <- expand.grid(eta_val = eta_alt, nu_val = nu_alt)
grid_beta <- grid_beta[-14,]
  
  
delta_seq <- seq(0, 1, length.out = 500)
theta_seq <- seq(-20, 20, by = 0.01)
par_grid <- expand.grid(delta = delta_seq, theta = theta_seq)
  
##  ............................................................................
##  Observed Bayes Factor                                                   ####
  
obs_log_bf <- c()
for (i in 1:nrow(grid_beta)) {
  obs_log_bf[i] <- round(bf_fun(delta = delta_seq, theta = theta, theta_0 = theta_0, 
                                st_dev = se, st_dev_0 = se_0, eta_null = eta_null, nu_null = nu_null,
                                eta_alt = grid_beta$eta_val[i], nu_alt = grid_beta$nu_val[i], 
                                log_val = TRUE),4)  
  grid_beta$obs_bf[i] <- obs_log_bf[i]
}
  
##  ............................................................................
##  Posterior Predictive                                                    ####

# MC Sampling
n_rep <- 200 # number of replications to generate
sampled_theta_list <- list()
y_rep <- list()
bf_list <- list()
theta_rep <- list()
sd_rep <- list()
post_theta <- list()
pb <- progress::progress_bar$new(
  format = "  Running [:bar] :percent eta: :eta",
  total = nrow(grid_beta), clear = FALSE, width= 60)
  
for (i in 1:nrow(grid_beta)) {
  pb$tick()
  post_theta[[i]] <- theta_m_post(theta_seq = theta_seq, theta = theta,
                                  theta_0 = theta_0, st_dev = se, st_dev_0 = se_0,
                                  eta = grid_beta$eta_val[i], nu = grid_beta$nu_val[i])
  y_rep[[i]] <- list()  # Initialize y_rep[[i]] as a list
  bf_list[[i]] <- list()  # Initialize bf_list[[i]] as a list
  theta_rep[[i]] <- list()
  sampled_theta_list[[i]] <- list()
  sd_rep[[i]] <- list()
  for (j in 1:n_rep) {
    sampled_theta_list[[i]][[j]] <- sample(theta_seq, size = 1, replace = TRUE, prob = (post_theta[[i]]/sum(post_theta[[i]])))
    y_rep[[i]][[j]] <- rnorm(200, mean = sampled_theta_list[[i]][[j]], sd = se) 
    theta_rep[[i]][[j]] <- round(mean(y_rep[[i]][[j]]),4)
    sd_rep[[i]][[j]] <- round(sd(y_rep[[i]][[j]]),4)
    bf_list[[i]][j] <-  round(bf_fun(delta_seq = delta_seq, theta = theta_rep[[i]][[j]], theta_0 = theta_0, 
                                     st_dev = sd_rep[[i]][[j]], st_dev_0 = se_0, eta_null = eta_null, nu_null = nu_null, 
                                     eta_alt = grid_beta$eta_val[i], nu_alt = grid_beta$nu_val[i], log_val = TRUE),4)
    
  }
}

bf_data <- data.frame()
for (i in 1:nrow(grid_beta)) {
  for (j in 1:n_rep) {
    bf_data[i,j] <- bf_list[[i]][[j]]
  }
}
  
# Define the directory path
path <- "R/Simulations/Gaussian/RData"
  
# Create the file name
file_name <- paste0("Gaussian_", formatC(M, format = "f", digits = 2), ".RData")
  
# Combine the path and file name to create a full file path
full_file_path <- file.path(path, file_name)
  
save.image(file = full_file_path)
  
}

