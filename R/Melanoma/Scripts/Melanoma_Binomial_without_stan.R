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
  normConst <- integrate(f = intFun, lower = 0, upper = 1, 
                          rel.tol = .Machine$double.eps^0.75)$value
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
      }, lower = 0, upper = 1, subdivisions = 500L, rel.tol = .Machine$double.eps^0.75)$value
  })
  theta_post <- numerator/norm_const
}



##  ............................................................................
##  Bayes Factor Functions                                                  ####

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
#   Data                                                                    ####

data <- read_excel("~/Desktop/Research_Project_PhD_Link/R_Projects/Power_Priors/Bayes_Factor/1.Main_Work_With_Ioannis/Uniform_Reference/Samples_From_Alternative/Melanoma/Data/trial_e1684_e1690_Merged.xlsx", 
                                       col_types = c("numeric", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric", "numeric", "numeric", 
                                                     "numeric"))


historical_data <- filter(data, study == 1684)[,-2]
current_data <- filter(data, study == 1690)[,-2]

N_0 <- nrow(historical_data)
N <- nrow(current_data)
y_0 <- 174
y <- 190

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
    y_rep[[i]][[j]] <- sum(rbinom(N, 1, sampled_theta_list[[i]][[j]]))
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


#   ____________________________________________________________________________
#   Optimal Prior                                                           ####

calculate_maximized_value <- function(distribution_values, observed_BF, hpdi_prob = 0.95) {
  results <- numeric(length(observed_BF)) # Initialize a results vector
  
  # Calculate the number of points in each distributions
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







# Reshape the data into long format
bf_data_long <- bf_data %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")
bf_data_long <- as.data.frame(bf_data_long)

optimal_model_index <- calculate_maximized_value(bf_data_long$Value,obs_log_bf, 0.95)

optimal_model <- print(grid_beta[optimal_model_index,])



#   ____________________________________________________________________________
#   Plots                                                                   ####

norm_const_alt <- marg_lik(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                       p = q, q = q, eta = 0.5, nu = 6)
post_theta_alt <- marg_post_theta(theta = theta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                   p = q, q = q, eta = 0.5, nu = 6, norm_const_alt)
post_delta_alt <- marg_post_delta(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                  p = q, q = q, eta = 0.5, nu = 6,  norm_const_alt)
norm_const_unif <- marg_lik(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                           p = q, q = q, eta = 1, nu = 1)
post_theta_unif <- marg_post_theta(theta = theta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                  p = q, q = q, eta = 1, nu = 1, norm_const_unif)
post_delta_unif <- marg_post_delta(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                  p = q, q = q, eta = 1, nu = 1,  norm_const_alt)


data_plot <- data.frame(theta_seq, delta_seq, 
                        post_theta_alt, post_delta_alt,
                        post_theta_unif, post_delta_unif)

# Plot using ggplot2"#E69F00", "#56B4E9"

ggplot(data_plot, aes(x = theta_seq)) +
  geom_line(aes(y = post_theta_alt, color = "Optimal"), size = 1) +
  geom_line(aes(y = post_theta_unif, color = "Default"), size = 1) +
  labs(x = expression(theta ~ "values"),  # Changing x-axis label
       y = "Posterior Probability",
       color = expression(pi[0](delta)~":")) +  # Changing legend title to pi(delta)_0
  theme_light() +
  scale_color_manual(values = c("Optimal" = "#E69F00", "Default" = "#56B4E9")) +
  coord_cartesian(xlim = c(0.35, 0.56))+
  theme(legend.position = "top")

prior_beta_1_1 <- dbeta(delta_seq, 1, 1)
prior_beta_0_5_6 <- dbeta(delta_seq, 0.5, 6)

# Add these to your data frame
data <- data.frame(delta_seq, post_delta_alt, post_delta_unif, prior_beta_1_1, prior_beta_0_5_6)

# Plot using ggplot2
ggplot(data, aes(x = delta_seq)) +
  geom_line(aes(y = post_delta_alt, color = "Optimal"), size = 1.2) +
  geom_line(aes(y = post_delta_unif, color = "Default"), size = 1.2) +
  geom_line(aes(y = prior_beta_1_1, color = "Default"), size = 1, linetype = "twodash") +
  geom_line(aes(y = prior_beta_0_5_6, color = "Optimal"), size = 1, linetype = "twodash") +
  labs(x = expression(delta ~ "values"),
       y = "Posterior Probability",
       color = expression(pi[0](delta)~":")) +
  theme_light() +
  scale_color_manual(values = c("Optimal" = "#E69F00", "Default" = "#56B4E9")) +
  coord_cartesian(xlim = c(0, 0.40)) +
  theme(legend.position = "top")



library(ggplot2)
library(ggpubr)

# Assuming data_plot and data are already defined and have the necessary columns.

# Create the first plot without a legend name and with a bold, centered title
plot1 <- ggplot(data_plot, aes(x = theta_seq)) +
  geom_line(aes(y = post_theta_alt, color = "Optimal"), size = 1) +
  geom_line(aes(y = post_theta_unif, color = "Default"), size = 1) +
  labs(x = expression(theta ~ "values"), y = "Posterior Probability", color = NULL) +
  #ggtitle(expression("Posterior Distribution of " * theta)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16)) +  # Bold and centered title
  scale_color_manual(values = c("Optimal" = "#8A0404", "Default" = "#56B4E9")) +
  coord_cartesian(xlim = c(0.35, 0.56))

# Create the second plot without a legend name and with a bold, centered title
plot2 <- ggplot(data, aes(x = delta_seq)) +
  geom_line(aes(y = post_delta_alt, color = "Optimal"), size = 1.2) +
  geom_line(aes(y = post_delta_unif, color = "Default"), size = 1.2) +
  geom_line(aes(y = prior_beta_1_1, color = "Default"), size = 1, linetype = "twodash") +
  geom_line(aes(y = prior_beta_0_5_6, color = "Optimal"), size = 1, linetype = "twodash") +
  labs(x = expression(delta ~ "values"), y = "Posterior Probability", color = NULL) +
  #ggtitle(expression("Prior and Posterior Distributions of " * delta)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16)) +  # Bold and centered title
  scale_color_manual(values = c("Optimal" = "#8A0404", "Default" = "#56B4E9")) +
  coord_cartesian(xlim = c(0, 0.40),
                  ylim = c(0,30))

# Combine the plots using ggarrange with a shared legend
plot_comb <- ggarrange(plot1, plot2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top")


ggsave(filename = "plot_melanoma_bin.pdf",path = "1.Main_Work_With_Ioannis/Uniform_Reference/Plots", plot = plot_comb,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)


