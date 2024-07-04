#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "ggthemes",
  "progressr",
  "ggplot2",
  "hypergeo",
  "colorspace",
  "modi")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#   ____________________________________________________________________________
#   Functions                                                               ####

# Beta-Binomial Function

dbetabinom <- function(y, N, a, b) {
  exp(lchoose(n = N, k = y) + lbeta(a + y, b + N - y) - lbeta(a, b))
}

# Marginal Likelihood Function

marg_lik <- function(delta, y, N, y_0, N_0, p, q, eta, nu) {
  intFun <- function(delta) {
    a_theta <- delta * y_0 + p
    b_theta <- delta * (N_0 - y_0) + q
    dbetabinom(y = y, N = N, a = a_theta, b = b_theta) *
      dbeta(x = delta, shape1 = eta, shape2 = nu)
  }
  normConst <- integrate(f = intFun, lower = 0, upper = 1)$value
}

# Marginal Posterior for delta

marg_post_delta <- function(delta, y, N, y_0, N_0, p, q, eta, nu, norm_const){
  numerator <- function(delta) {
    a_theta <- delta * y_0 + p
    b_theta <- delta * (N_0 - y_0) + q
    dbetabinom(y = y, N = N, a = a_theta, b = b_theta) *
      dbeta(x = delta, shape1 = eta, shape2 = nu)
  }
  dens <- numerator(delta = delta) / norm_const
}

# Marginal Posterior for theta

marg_post_theta <- function(theta, y, N, y_0, N_0, p, q, eta, nu, norm_const){
  numerator <- sapply(X = theta, FUN = function(t){
    dbinom(y, size=N, prob=t) * 
    integrate(f = function(delta){
      dbeta(x = delta, shape1 = eta, shape2 = nu) * 
       dbeta(x = t, shape1 = (delta * y_0 + p), shape2 = (delta * (N_0 - y_0) + q))
    }, lower = 0, upper = 1, subdivisions = 300L, rel.tol = .Machine$double.eps^0.55)$value
  })
  theta_post <- numerator/norm_const
}


# Observed Bayes Factor


bf_fun <- function(delta_seq, y, N, y_0, N_0, p, q, eta_null, nu_null, eta_alt, nu_alt, log_val = TRUE){
  # Marginal Likelihoods
  
  marg_lik_null <- marg_lik(delta, y, N, y_0, N_0, p, q, eta_null, nu_null)
  marg_lik_alt <- marg_lik(delta, y, N, y_0, N_0, p, q, eta_alt, nu_alt)
  
  # Log-value assessment 
  if(log_val == TRUE){
    bayes_factor <- log(marg_lik_alt/marg_lik_null, 10) 
  } 
  else{
    bayes_factor <- marg_lik_alt/marg_lik_null 
  }
  
  # Final value
  return(bayes_factor)
}



#   ____________________________________________________________________________
#   Simulation Setting                                                      ####


par_seq <-  seq(20,45, by = 1)

for (M in 1:length(par_seq)) {
  
  cat("Doing theta =",
      par_seq[M],
      "\n")


#   ____________________________________________________________________________
#   Data                                                                    ####
  
set.seed(4231)
# Historical Data Parameters
N_0 <- 100
y_0 <- 20
p_0 <- y_0 / N_0

# Current Data Parameters
N <- 100
y <- par_seq[M]
p <- y / N



df <- rbind(data.frame(binom=rbinom(N_0, y_0, p_0), binomial.p='Historical') , data.frame(binom=rbinom(N, y, p), binomial.p='Current'))


#   ____________________________________________________________________________
#   Model and Grids Specification                                           ####

# Initial Prior Parameter for Theta
p <- 1
q <- 1

# Initial Prior Parameter for Beta
eta_null <- 1
nu_null <- 1

eta_alt <- seq(0.5,6,by = 0.5)
nu_alt <- seq(0.5,6,by = 0.5)

grid_beta <- expand.grid(eta_val = eta_alt, nu_val = nu_alt)
grid_beta <- grid_beta[-14,]


delta_seq <- seq(0, 1, length.out = 1000)
theta_seq <- seq(0, 1, length.out = 1000)
par_grid <- expand.grid(delta = delta_seq, theta = theta_seq)



#   ____________________________________________________________________________
#   Observed Bayes Factor                                                   ####


obs_log_bf <- c()
for (i in 1:nrow(grid_beta)) {
  obs_log_bf[i] <- round(bf_fun(delta_seq = delta_seq, y = y, N = N, y_0 = y_0,
                                N_0 = N_0, p = p, q = q, eta_null = eta_null , nu_null = nu_null, 
                                eta_alt = grid_beta$eta_val[i], nu_alt = grid_beta$nu_val[i],
                                log_val = TRUE),4)  
  grid_beta$obs_bf[i] <- obs_log_bf[i]
}

#   ____________________________________________________________________________
#   Posterior Predictive                                                    ####

# MC Sampling

n_rep <- 100 # number of replications to generate
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
  norm_const <- marg_lik(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                         p = q, q = q, eta = grid_beta$eta_val[i], nu = grid_beta$nu_val[i])
  post_theta[[i]] <- marg_post_theta(theta = theta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                     p = q, q = q, eta = grid_beta$eta_val[i], nu = grid_beta$nu_val[i], norm_const)
  y_rep[[i]] <- list()  # Initialize y_rep[[i]] as a list
  bf_list[[i]] <- list()  # Initialize bf_list[[i]] as a list
  theta_rep[[i]] <- list()
  sampled_theta_list[[i]] <- list()
  sd_rep[[i]] <- list()
  for (j in 1:n_rep) {
      sampled_theta_list[[i]][[j]] <- sample(theta_seq, size = 1, replace = TRUE, prob = (post_theta[[i]]/sum(post_theta[[i]])))
      y_rep[[i]][[j]] <- sum(rbinom(100,1, sampled_theta_list[[i]][[j]]))
      bf_list[[i]][j] <-  round(bf_fun(delta_seq = delta_seq, y = y_rep[[i]][[j]], N = N, y_0 = y_0,
                                       N_0 = N_0, p = p, q = q, eta_null = eta_null , nu_null = nu_null, 
                                       eta_alt = grid_beta$eta_val[i], nu_alt = grid_beta$nu_val[i],
                                       log_val = TRUE),4)
      
    }
}


bf_data <- data.frame()
for (i in 1:nrow(grid_beta)) {
  for (j in 1:n_rep) {
    bf_data[i,j] <- bf_list[[i]][[j]]
  }
}

# Define the directory path
path <- "R/Simulations/Bernoulli/RData"

# Create the file name
file_name <- paste0("Bernoulli_", formatC(M, format = "f", digits = 2), ".RData")

# Combine the path and file name to create a full file path
full_file_path <- file.path(path, file_name)

save.image(file = full_file_path)

}

