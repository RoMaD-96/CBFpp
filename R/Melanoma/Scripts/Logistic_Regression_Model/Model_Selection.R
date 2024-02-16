#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "rethinking",
  "ggpubr",
  "rstan",
  "bridgesampling",
  "dplyr",
  "readxl",
  "progressr",
  "ggplot2",
  "ggridges",
  "tidyverse",
  "tidyr",
  "tibble"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


source("R/CBF_Functions.R")


load("R/Melanoma/RData/Logit_VM.RData")
rm(list = setdiff(ls(), c("df_obs_bf","grid")))



#   ____________________________________________________________________________
#   Optimal Model                                                           ####

data_structure <- function(path, n_files){
  
  base_path <- path
  
  # Load Data
  for (i in 1:n_files) {
    # Construct the full file path
    file_path <- paste0(base_path, "bf_list_", i, ".RData")
    load(file_path)
    assign(paste0("bf_list_", i), bf_list)
  }
  
  
  combined_list <- list()
  
  for (i in 1:n_files) {
    # Construct the name of the list
    list_name <- paste0("bf_list_", i)
    
    # Fetch the list using the 'get' function
    current_list <- get(list_name)
    
    # Add the contents of the current list to the combined list
    combined_list <- c(combined_list, current_list)
  }
  
  original_list <- list()
  original_list <- combined_list 
  
  # Number of models
  num_models <- 143
  
  # Creare una nuova lista vuota per contenere gli elementi raggruppati
  grouped_list <- vector("list", num_models)
  
  for(i in 1:num_models){
    indices <- seq(i, length(original_list), num_models)
    grouped_list[[i]] <- original_list[indices]
  }
  
  
  grouped_list <- lapply(grouped_list[1:num_models], function(x) unlist(x, recursive = FALSE))
  
  n_bf_values <- length(grouped_list[[1]])
  
  bf_data_rep <- data.frame()
  for (i in 1:num_models) {
    for (j in 1:n_bf_values) {
      bf_data_rep[i,j] <- grouped_list[[i]][[j]]$bf
    }
  }
  
  return(bf_data_rep)
}

path <- "R/Melanoma/RData/"


bf_data <- data_structure(path,48)

# Reshape the data into long format
bf_data_long <- bf_data %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")



optimal_model_index <- calculate_maximized_value(bf_data_long$Value, df_obs_bf$obs_bf, 0.95)



#   ____________________________________________________________________________
#   Datasets                                                                ####

set.seed(478)

data <- read_excel("Data/trial_e1684_e1690_Merged.xlsx", 
                   col_types = c("numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                "numeric"))



historical_data <- filter(data, study == 1684)[,-2]
current_data <- filter(data, study == 1690)[,-2]

## Standardization ##

log_age_hist <- log(historical_data$age)
log_age_current <- log(current_data$age)


#   ____________________________________________________________________________
#   Approximately Normalized Prior                                          ####

### Grid Size ###
J <- 20

### GAM Approximation ###
constant_data <- read.csv("Data/delta_estimates_logit.csv")
fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

# mgcv::plot.gam(fit_gam)

#   ____________________________________________________________________________
#   STAN Data Block Configuration                                           ####


N_0 <- length(log_age_hist)
X_0 <- cbind(log_age = log_age_hist, # x1 - log age 
             sex = historical_data$sex, # x2 - gender 
             treat_status = historical_data$trt # x3 treatment
)
Y_0_cens <- historical_data$survtime
Cens_0 <- historical_data$scens


N<- length(log_age_current)
X <- cbind(log_age_current, # x1 - log age 
           current_data$sex, # x2 - gender
           current_data$trt # x3 treatment
)
Y_cens <- current_data$survtime
Cens<- current_data$scens


#   ____________________________________________________________________________
#   Model Specification                                                     ####

stan_model <- list("R/Melanoma/STAN/Logit_Sim_Norm_PP.stan")

# Hyperparameters

alpha_0 <- 1
beta_0 <- 1

# Stan lists

# CBF Prior
opt_bf_model_list <- list(
                          N0 = N_0,
                          P = ncol(X_0),
                          X0 = X_0,
                          y0 = Cens_0,
                          N = N,
                          X = X,
                          y = Cens,
                          eta = 0.5,
                          nu = 6)
# Uniform Prior
unif_model_list <- list(
                        N0 = N_0,
                        P = ncol(X_0),
                        X0 = X_0,
                        y0 = Cens_0,
                        N = N,
                        X = X,
                        y = Cens,
                        eta = 1,
                        nu = 1) 
# KL Prior
KL_model_list <- list(
  N0 = N_0,
  P = ncol(X_0),
  X0 = X_0,
  y0 = Cens_0,
  N = N,
  X = X,
  y = Cens,
  eta = 0.7,
  nu = 1.5) 

# MSE Prior
MSE_model_list <- list(
  N0 = N_0,
  P = ncol(X_0),
  X0 = X_0,
  y0 = Cens_0,
  N = N,
  X = X,
  y = Cens,
  eta = 5.5,
  nu = 3) 

# Jeffreys Prior
jeffreys_model_list <- list(
  N0 = N_0,
  P = ncol(X_0),
  X0 = X_0,
  y0 = Cens_0,
  N = N,
  X = X,
  y = Cens,
  eta = 0.5,
  nu = 0.5) 

#   ____________________________________________________________________________
#   Running models                                                          ####


rstan_options(auto_write = TRUE)
options(mc.cores = 4)

Ks <- c(10000)

J <- 20

### GAM Approximation
constant_data <- read.csv("Data/delta_estimates_logit.csv")
fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)

### Normalizing Power Prior

model_opt_bf <- norm_post_pp(Ks, 
                           list(opt_bf_model_list), 
                           stan_model)
model_unif <- norm_post_pp(Ks,
                           list(unif_model_list),
                           stan_model)
model_KL <- norm_post_pp(Ks,
                           list(KL_model_list),
                           stan_model)
model_MSE <- norm_post_pp(Ks,
                         list(MSE_model_list),
                         stan_model)
model_jeffreys <- norm_post_pp(Ks,
                         list(jeffreys_model_list),
                         stan_model)
 
# Posterior Parameters
model_opt_bf_par <- rstan::extract(model_opt_bf[[1]])
model_unif_par <- rstan::extract(model_unif[[1]])
model_KL_par <- rstan::extract(model_KL[[1]])
model_MSE_par <- rstan::extract(model_MSE[[1]])
model_jeffreys_par <- rstan::extract(model_jeffreys[[1]])


#   ____________________________________________________________________________
#   Latex Table                                                             ####


# Dataframe Posterior Parameters
data_parameter_post <- data.frame(log_age_beta_opt = model_opt_bf_par$beta[,1],
                                  log_age_beta_unif = model_unif_par$beta[,1],
                                  log_age_beta_KL = model_KL_par$beta[,1],
                                  log_age_beta_MSE = model_MSE_par$beta[,1],
                                  log_age_beta_jeffreys = model_jeffreys_par$beta[,1],
                                  sex_beta_opt = model_opt_bf_par$beta[,2],
                                  sex_beta_unif = model_unif_par$beta[,2],
                                  sex_beta_KL = model_KL_par$beta[,2],
                                  sex_beta_MSE = model_MSE_par$beta[,2],
                                  sex_beta_jeffreys = model_jeffreys_par$beta[,2],
                                  treatment_beta_opt = model_opt_bf_par$beta[,3],
                                  treatment_beta_unif = model_unif_par$beta[,3],
                                  treatment_beta_KL = model_KL_par$beta[,3],
                                  treatment_beta_MSE = model_MSE_par$beta[,3],
                                  treatment_beta_jeffreys = model_jeffreys_par$beta[,3],
                                  delta_opt = model_opt_bf_par$delta,
                                  delta_unif = model_unif_par$delta,
                                  delta_KL = model_KL_par$delta,
                                  delta_MSE = model_MSE_par$delta,
                                  delta_opt_jeffreys = model_jeffreys_par$delta) 

# HPDI 
hpdi_age_opt <- rethinking::HPDI( data_parameter_post$log_age_beta_opt, prob = 0.95)
hpdi_age_unif <- rethinking::HPDI( data_parameter_post$log_age_beta_unif, prob = 0.95)
hpdi_age_KL <- rethinking::HPDI( data_parameter_post$log_age_beta_KL, prob = 0.95)
hpdi_age_MSE <- rethinking::HPDI( data_parameter_post$log_age_beta_MSE, prob = 0.95)
hpdi_age_jeffreys <- rethinking::HPDI( data_parameter_post$log_age_beta_jeffreys, prob = 0.95)
hpdi_sex_opt <- rethinking::HPDI( data_parameter_post$sex_beta_opt, prob = 0.95)
hpdi_sex_unif <- rethinking::HPDI( data_parameter_post$sex_beta_unif, prob = 0.95)
hpdi_sex_KL <- rethinking::HPDI( data_parameter_post$sex_beta_KL, prob = 0.95)
hpdi_sex_MSE <- rethinking::HPDI( data_parameter_post$sex_beta_MSE, prob = 0.95)
hpdi_sex_jeffreys <- rethinking::HPDI( data_parameter_post$sex_beta_jeffreys, prob = 0.95)
hpdi_treat_opt <- rethinking::HPDI( data_parameter_post$treatment_beta_opt, prob = 0.95)
hpdi_treat_unif <- rethinking::HPDI( data_parameter_post$treatment_beta_unif, prob = 0.95)
hpdi_treat_KL <- rethinking::HPDI( data_parameter_post$treatment_beta_KL, prob = 0.95)
hpdi_treat_MSE <- rethinking::HPDI( data_parameter_post$treatment_beta_MSE, prob = 0.95)
hpdi_treat_jeffreys <- rethinking::HPDI( data_parameter_post$treatment_beta_jeffreys, prob = 0.95)


# Dataframe for Latex Table
data_parameter_post_tab <- data.frame(
  prior = c("Beta(0.5,6)", "Beta(1,1)", "Beta(0.7,1.5)", "Beta(5.5,3)", "Beta(0.5,0.5)"),
  mean_age = format(round(c(mean(data_parameter_post$log_age_beta_opt),
           mean(data_parameter_post$log_age_beta_unif),
           mean(data_parameter_post$log_age_beta_KL),
           mean(data_parameter_post$log_age_beta_MSE),
           mean(data_parameter_post$log_age_beta_jeffreys)),2), nsmall = 2),
  mean_sex = format(round(c(mean(data_parameter_post$sex_beta_opt),
               mean(data_parameter_post$sex_beta_unif),
               mean(data_parameter_post$sex_beta_KL),
               mean(data_parameter_post$sex_beta_MSE),
               mean(data_parameter_post$sex_beta_jeffreys)),2), nsmall = 2),
  mean_treat = format(round(c(mean(data_parameter_post$treatment_beta_opt),
                mean(data_parameter_post$treatment_beta_unif),
                mean(data_parameter_post$treatment_beta_KL),
                mean(data_parameter_post$treatment_beta_MSE),
                mean(data_parameter_post$treatment_beta_jeffreys)
                ),2), nsmall = 2),
  sd_age = format(round(c(sd(data_parameter_post$log_age_beta_opt),
         sd(data_parameter_post$log_age_beta_unif),
         sd(data_parameter_post$log_age_beta_KL),
         sd(data_parameter_post$log_age_beta_MSE),
         sd(data_parameter_post$log_age_beta_jeffreys)),2), nsmall = 2),
  sd_sex = format(round(c(sd(data_parameter_post$sex_beta_opt),
             sd(data_parameter_post$sex_beta_unif),
             sd(data_parameter_post$sex_beta_KL),
             sd(data_parameter_post$sex_beta_MSE),
             sd(data_parameter_post$sex_beta_jeffreys)),2), nsmall = 2),
  sd_treat = format(round(c(sd(data_parameter_post$treatment_beta_opt),
              sd(data_parameter_post$treatment_beta_unif),
              sd(data_parameter_post$treatment_beta_KL),
              sd(data_parameter_post$treatment_beta_MSE),
              sd(data_parameter_post$treatment_beta_jeffreys)),2), nsmall = 2),
  hpdi_age = c(paste("(", round(hpdi_age_opt[1],2), ",", round(hpdi_age_opt[2],2), ")"),
           paste("(", round(hpdi_age_unif[1],2), ",", round(hpdi_age_unif[2],2), ")"),
           paste("(", round(hpdi_age_KL[1],2), ",", round(hpdi_age_KL[2],2), ")"),
           paste("(", round(hpdi_age_MSE[1],2), ",", round(hpdi_age_MSE[2],2), ")"),
           paste("(", round(hpdi_age_jeffreys[1],2), ",", round(hpdi_age_jeffreys[2],2), ")")),
  hpdi_sex = c(paste("(", round(hpdi_sex_opt[1],2), ",", round(hpdi_sex_opt[2],2), ")"),
           paste("(", round(hpdi_sex_unif[1],2), ",", round(hpdi_sex_unif[2],2), ")"),
           paste("(", round(hpdi_sex_KL[1],2), ",", round(hpdi_sex_KL[2],2), ")"),
           paste("(", round(hpdi_sex_MSE[1],2), ",", round(hpdi_sex_MSE[2],2), ")"),
           paste("(", round(hpdi_sex_jeffreys[1],2), ",", round(hpdi_sex_jeffreys[2],2), ")")),
  hpdi_treat = c(paste("(", round(hpdi_treat_opt[1],2), ",", round(hpdi_treat_opt[2],2), ")"),
                paste("(", round(hpdi_treat_unif[1],2), ",", round(hpdi_treat_unif[2],2), ")"),
                paste("(", round(hpdi_treat_KL[1],2), ",", round(hpdi_treat_KL[2],2), ")"),
                paste("(", round(hpdi_treat_MSE[1],2), ",", round(hpdi_treat_MSE[2],2), ")"),
                paste("(", round(hpdi_treat_jeffreys[1],2), ",", round(hpdi_treat_jeffreys[2],2), ")"))
)

# Create LaTeX Table
dfTab_beta <- data_parameter_post_tab

xtab_beta <- xtable(dfTab_beta)
colnames(xtab_beta) <- c("",
                          "Age",
                          "Sex",
                          "Treat.",
                          "Age",
                          "Sex",
                          "Treat.",
                          "Age",
                          "Sex",
                          "Treat.")

align(xtab_beta) <- rep("c", length(colnames(xtab_beta)) + 1)

## Add Multicolumns
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\toprule & & & \\multicolumn{2}{c}{Tests about the effect size $\\theta$} & \\multicolumn{2}{c}{Tests about the power parameter $\\alpha$} \\\\ \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}'
print(xtab_beta, floating = FALSE, include.rownames = FALSE, add.to.row = addtorow,
      sanitize.text.function = function(x){x}, booktabs = TRUE, hline.after = c(0, nrow(xtab_beta)))


#   ____________________________________________________________________________
#   Plots                                                                   ####


# Log-Age Plots
long_data_age <- gather(data_parameter_post, key = "variable", value = "value", 
                        log_age_beta_opt, log_age_beta_unif)

plot_age <- ggplot(long_data_age, aes(x = value, fill = variable)) +
  geom_density(data = subset(long_data_age, variable == "log_age_beta_unif"), 
               aes(x = value, fill = variable), alpha = 0.95) +
  geom_density(data = subset(long_data_age, variable == "log_age_beta_opt"), 
               aes(x = value, fill = variable), alpha = 0.8) +
  labs(title = "Log of Age", x = "", y = "Posterior Density") +
  scale_fill_manual(values = c("#D55E00", "grey20"), 
                    name = "Prior",
                    labels = c("Optimal", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_age_opt[1], xmax = hpdi_age_opt[2], y = 2.2 * 1.15, height = 0.15), 
                 color = "#D55E00", alpha = 1, size = 0.5) +
  geom_errorbarh(aes(xmin = hpdi_age_unif[1], xmax = hpdi_age_unif[2], y = 2.2 * 1.05, height = 0.15), 
                 color = "grey20", alpha = 1, size = 0.5) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.8)) +  # Adjust y-axis
  scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5), limits = c(-1.0, 1.2)) +  # Adjust x-axis
  theme_minimal(base_size = 16)+
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))

# Print the Plot
print(plot_age)



# Sex Plots
long_data_sex <- gather(data_parameter_post, key = "variable", value = "value", 
                        sex_beta_opt, sex_beta_unif)

plot_sex <- ggplot(long_data_sex, aes(x = value, fill = variable)) +
  geom_density(data = subset(long_data_sex, variable == "sex_beta_unif"), 
               aes(x = value, fill = variable), alpha = 0.95) +
  geom_density(data = subset(long_data_sex, variable == "sex_beta_opt"), 
               aes(x = value, fill = variable), alpha = 0.8) +
  labs(title = "Sex", x = "Value", y = "") +
  scale_fill_manual(values = c("#D55E00", "grey20"), 
                    name = "Prior",
                    labels = c("Optimal", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_sex_opt[1], xmax = hpdi_sex_opt[2], y = 2.2 * 1.15, height = 0.15), 
                 color = "#D55E00", alpha = 1, size = 0.5) +
  geom_errorbarh(aes(xmin = hpdi_sex_unif[1], xmax = hpdi_sex_unif[2], y = 2.2 * 1.05, height = 0.15), 
                 color = "grey20", alpha = 1, size = 0.5) +
  theme_minimal(base_size = 16) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.8)) +  # Adjust y-axis
  scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5), limits = c(-1.0, 1.2)) +  # Adjust x-axis
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))

# Print the Plot
print(plot_sex)


# Treatment Plots
long_data_treat <- gather(data_parameter_post, key = "variable", value = "value", 
                        treatment_beta_opt, treatment_beta_unif)

plot_treat <- ggplot(long_data_treat, aes(x = value, fill = variable)) +
  geom_density(data = subset(long_data_treat, variable == "treatment_beta_unif"), 
               aes(x = value, fill = variable), alpha = 0.95) +
  geom_density(data = subset(long_data_treat, variable == "treatment_beta_opt"), 
               aes(x = value, fill = variable), alpha = 0.8) +
  labs(title = "Treatment", x = "", y = "") +
  scale_fill_manual(values = c("#D55E00", "grey20"), 
                    name = "Prior",
                    labels = c("Optimal", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_treat_opt[1], xmax = hpdi_treat_opt[2], y = 2.2 * 1.15, height = 0.15), 
                 color = "#D55E00", alpha = 1, size = 0.5) +
  geom_errorbarh(aes(xmin = hpdi_treat_unif[1], xmax = hpdi_treat_unif[2], y = 2.2 * 1.05, height = 0.15), 
                 color = "grey20", alpha = 1, size = 0.5) +
  theme_minimal(base_size = 16) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.8)) +  # Adjust y-axis
  scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5), limits = c(-1.0, 1.2)) +  # Adjust x-axis
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))

# Print the plot
print(plot_treat)

plot_comb <- ggarrange(plot_age, plot_sex, plot_treat, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
print(plot_comb)

ggsave(filename = "post_param_melanoma.pdf",path = "Plots", plot = plot_comb,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)





# Delta plot
long_data_delta <- gather(data_parameter_post, key = "variable", value = "value", 
                          delta_opt, delta_unif, delta_KL, delta_MSE, delta_opt_jeffreys)
hpdi_delta_opt <- rethinking::HPDI( data_parameter_post$delta_opt, prob = 0.95)
hpdi_delta_unif <- rethinking::HPDI( data_parameter_post$delta_unif, prob = 0.95)
plot_delta <- ggplot(long_data_delta, aes(x = value, fill = variable)) +
  geom_density(data = subset(long_data_delta, variable == "delta_unif"), 
               aes(x = value, fill = variable), alpha = 0.95) +
  geom_density(data = subset(long_data_delta, variable == "delta_opt"), 
               aes(x = value, fill = variable), alpha = 0.8) +
  labs(x = expression("Values of "* delta), y = "Posterior Density") +
  scale_fill_manual(values = c("#D55E00", "grey20"), 
                    name = "Prior: ",
                    labels = c("Optimal", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_delta_opt[1], xmax = hpdi_delta_opt[2], y = 1.8 * 1.15, height = 0.15), 
                 color = "#D55E00", alpha = 1, size = 0.5) +
  geom_errorbarh(aes(xmin = hpdi_delta_unif[1], xmax = hpdi_delta_unif[2], y = 1.8 * 1.05, height = 0.15), 
                 color = "grey20", alpha = 1, size = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")

# Print the plot
print(plot_delta)


