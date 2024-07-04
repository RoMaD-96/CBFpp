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
  "xtable",
  "tibble"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


source("R/CBF_Functions.R")
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
  num_models <- 64
  
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

path <- "R/Melanoma/RData/Melanoma_Data/"


bf_data <- data_structure(path,40)

# Reshape the data into long format
bf_data_long <- bf_data %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")

load("R/Melanoma/obs_bf_melanoma.RData")
df_obs_bf <- filter(df_obs_bf, obs_bf>=0)


optimal_model_index <- calculate_maximized_value(bf_data_long$Value, df_obs_bf$obs_bf, 0.75)



#   ____________________________________________________________________________
#   Datasets                                                                ####

set.seed(4231)
set.seed(433)

library(hdbayes)
data("E2696")
data("E1694")
historical_data <- E2696
current_data <- E1694


#   ____________________________________________________________________________
#   Approximately Normalized Prior                                          ####

### Grid Size ###
J <- 20

### GAM Approximation ###
constant_data <- read.csv("Data/delta_estimates_logit.csv")
fit_gam <- mgcv::gam(lc_a0 ~ s(a0, k = J + 1), data = constant_data)
plot(fit_gam)
 #mgcv::plot.gam(fit_gam)

#   ____________________________________________________________________________
#   STAN Data Block Configuration                                           ####


N_0 <- nrow(historical_data)
X0 <- cbind(historical_data$age, 
            historical_data$treatment, 
            historical_data$sex, 
            historical_data$perform 
)


N <- nrow(current_data)
X <- cbind(current_data$age, 
           current_data$treatment,
           current_data$sex,
           current_data$perform
)
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
  P = ncol(X0),
  X0 = as.matrix(X0),
  y0 = historical_data$failind,
  N = N,
  X = as.matrix(X),
  y = current_data$failind,
  eta = 5,
  nu = 0.5
)
# Uniform Prior
unif_model_list <- list(
  N0 = N_0,
  P = ncol(X0),
  X0 = as.matrix(X0),
  y0 = historical_data$failind,
  N = N,
  X = as.matrix(X),
  y = current_data$failind,
  eta = 1,
  nu = 1
)

# Jeffreys Prior
jeffreys_model_list <- list(
  N0 = N_0,
  P = ncol(X0),
  X0 = as.matrix(X0),
  y0 = historical_data$failind,
  N = N,
  X = as.matrix(X),
  y = current_data$failind,
  eta = 0.5,
  nu = 0.5) 

# Beta(2,2)

beta_2_model_list <- list(
  N0 = N_0,
  P = ncol(X0),
  X0 = as.matrix(X0),
  y0 = historical_data$failind,
  N = N,
  X = as.matrix(X),
  y = current_data$failind,
  eta = 2,
  nu = 2) 




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

model_jeffreys <- norm_post_pp(Ks,
                         list(jeffreys_model_list),
                         stan_model)

model_beta_2 <- norm_post_pp(Ks,
                               list(beta_2_model_list),
                               stan_model) 
# Posterior Parameters
model_opt_bf_par <- rstan::extract(model_opt_bf[[1]])
model_unif_par <- rstan::extract(model_unif[[1]])
model_jeffreys_par <- rstan::extract(model_jeffreys[[1]])
model_beta_2_par <- rstan::extract(model_beta_2[[1]])

#   ____________________________________________________________________________
#   Latex Table                                                             ####

# Dataframe Posterior Parameters
data_parameter_post <- data.frame(age_beta_opt = model_opt_bf_par$beta[,1],
                                  age_beta_unif = model_unif_par$beta[,1],
                                  age_beta_jeffreys = model_jeffreys_par$beta[,1],
                                  age_beta_2 = model_beta_2_par$beta[,1],
                                  treat_beta_opt = model_opt_bf_par$beta[,2],
                                  treat_beta_unif = model_unif_par$beta[,2],
                                  treat_beta_jeffreys = model_jeffreys_par$beta[,2],
                                  treat_beta_2 = model_beta_2_par$beta[,2],
                                  sex_beta_opt = model_opt_bf_par$beta[,3],
                                  sex_beta_unif = model_unif_par$beta[,3],
                                  sex_beta_jeffreys = model_jeffreys_par$beta[,3],
                                  sex_beta_2 = model_beta_2_par$beta[,3],
                                  perform_beta_opt = model_opt_bf_par$beta[,4],
                                  perform_beta_unif = model_unif_par$beta[,4],
                                  perform_beta_jeffreys = model_jeffreys_par$beta[,4],
                                  perform_beta_2 = model_beta_2_par$beta[,4],
                                  delta_opt = model_opt_bf_par$delta,
                                  delta_unif = model_unif_par$delta,
                                  delta_opt_jeffreys = model_jeffreys_par$delta, 
                                  delta_beta_2 = model_beta_2_par$delta) 

# HPDI 
hpdi_age_opt <- rethinking::HPDI( data_parameter_post$age_beta_opt, prob = 0.95)
hpdi_age_unif <- rethinking::HPDI( data_parameter_post$age_beta_unif, prob = 0.95)
hpdi_age_jeffreys <- rethinking::HPDI( data_parameter_post$age_beta_jeffreys, prob = 0.95)
hpdi_age_beta_2 <- rethinking::HPDI( data_parameter_post$age_beta_2, prob = 0.95)
hpdi_treat_opt <- rethinking::HPDI( data_parameter_post$treat_beta_opt, prob = 0.95)
hpdi_treat_unif <- rethinking::HPDI( data_parameter_post$treat_beta_unif, prob = 0.95)
hpdi_treat_jeffreys <- rethinking::HPDI( data_parameter_post$treat_beta_jeffreys, prob = 0.95)
hpdi_treat_beta_2 <- rethinking::HPDI( data_parameter_post$treat_beta_2, prob = 0.95)
hpdi_sex_opt <- rethinking::HPDI( data_parameter_post$sex_beta_opt, prob = 0.95)
hpdi_sex_unif <- rethinking::HPDI( data_parameter_post$sex_beta_unif, prob = 0.95)
hpdi_sex_jeffreys <- rethinking::HPDI( data_parameter_post$sex_beta_jeffreys, prob = 0.95)
hpdi_sex_beta_2 <- rethinking::HPDI( data_parameter_post$sex_beta_2, prob = 0.95)
hpdi_perform_opt <- rethinking::HPDI( data_parameter_post$perform_beta_opt, prob = 0.95)
hpdi_perform_unif <- rethinking::HPDI( data_parameter_post$perform_beta_unif, prob = 0.95)
hpdi_perform_jeffreys <- rethinking::HPDI( data_parameter_post$perform_beta_jeffreys, prob = 0.95)
hpdi_perform_beta_2 <- rethinking::HPDI( data_parameter_post$perform_beta_2, prob = 0.95)





# Dataframe for Latex Table
data_parameter_post_tab <- data.frame(
  prior = c("Beta(0.5,5)", "Beta(1,1)", 
            #"Beta(0.7,1.5)", "Beta(5.5,3)",
            "Beta(0.5,0.5)", "Beta(2,2)"),
  mean_age = format(round(c(mean(data_parameter_post$age_beta_opt),
           mean(data_parameter_post$age_beta_unif),
           mean(data_parameter_post$age_beta_jeffreys),
           mean(data_parameter_post$age_beta_2)),3), nsmall = 3),
  mean_treat = format(round(c(mean(data_parameter_post$treat_beta_opt),
               mean(data_parameter_post$treat_beta_unif),
               mean(data_parameter_post$treat_beta_jeffreys),
               mean(data_parameter_post$treat_beta_2)),3), nsmall = 3),
  mean_sex = format(round(c(mean(data_parameter_post$sex_beta_opt),
                mean(data_parameter_post$sex_beta_unif),
                mean(data_parameter_post$sex_beta_jeffreys),
                mean(data_parameter_post$sex_beta_2)
                ),3), nsmall = 3),
  mean_perform = format(round(c(mean(data_parameter_post$perform_beta_opt),
                            mean(data_parameter_post$perform_beta_unif),
                            mean(data_parameter_post$perform_beta_jeffreys),
                            mean(data_parameter_post$perform_beta_2)
  ),3), nsmall = 3),
  sd_age = format(round(c(sd(data_parameter_post$age_beta_opt),
         sd(data_parameter_post$age_beta_unif),
         sd(data_parameter_post$age_beta_jeffreys),
         sd(data_parameter_post$age_beta_2)),3), nsmall = 3),
  sd_treat = format(round(c(sd(data_parameter_post$treat_beta_opt),
             sd(data_parameter_post$treat_beta_unif),
             sd(data_parameter_post$treat_beta_jeffreys),
             sd(data_parameter_post$treat_beta_2)),3), nsmall = 3),
  sd_sex = format(round(c(sd(data_parameter_post$sex_beta_opt),
              sd(data_parameter_post$sex_beta_unif),
              sd(data_parameter_post$sex_beta_jeffreys),
              sd(data_parameter_post$sex_beta_2)),3), nsmall = 3),
  sd_perform = format(round(c(sd(data_parameter_post$perform_beta_opt),
                          sd(data_parameter_post$perform_beta_unif),
                          sd(data_parameter_post$perform_beta_jeffreys),
                          sd(data_parameter_post$perform_beta_2)),3), nsmall = 3),
  
  hpdi_age = c(paste("(", round(hpdi_age_opt[1],3), ",", round(hpdi_age_opt[2],3), ")"),
           paste("(", round(hpdi_age_unif[1],3), ",", round(hpdi_age_unif[2],3), ")"),
           paste("(", round(hpdi_age_jeffreys[1],3), ",", round(hpdi_age_jeffreys[2],3), ")"),
           paste("(", round(hpdi_age_beta_2[1],3), ",", round(hpdi_age_beta_2[2],3), ")")),
  hpdi_treat = c(paste("(", round(hpdi_treat_opt[1],3), ",", round(hpdi_treat_opt[2],3), ")"),
           paste("(", round(hpdi_treat_unif[1],3), ",", round(hpdi_treat_unif[2],3), ")"),
           paste("(", round(hpdi_treat_jeffreys[1],3), ",", round(hpdi_treat_jeffreys[2],3), ")"),
           paste("(", round(hpdi_treat_beta_2[1],3), ",", round(hpdi_treat_beta_2[2],3), ")")),
  hpdi_sex = c(paste("(", round(hpdi_sex_opt[1],3), ",", round(hpdi_sex_opt[2],3), ")"),
                paste("(", round(hpdi_sex_unif[1],3), ",", round(hpdi_sex_unif[2],3), ")"),
                paste("(", round(hpdi_sex_jeffreys[1],3), ",", round(hpdi_sex_jeffreys[2],3), ")"),
               paste("(", round(hpdi_sex_beta_2[1],3), ",", round(hpdi_sex_beta_2[2],3), ")")),
  hpdi_perform = c(paste("(", round(hpdi_perform_opt[1],3), ",", round(hpdi_perform_opt[2],3), ")"),
               paste("(", round(hpdi_perform_unif[1],3), ",", round(hpdi_perform_unif[2],3), ")"),
               paste("(", round(hpdi_perform_jeffreys[1],3), ",", round(hpdi_perform_jeffreys[2],3), ")"),
               paste("(", round(hpdi_perform_beta_2[1],3), ",", round(hpdi_perform_beta_2[2],3), ")"))
)

# Create LaTeX Table
dfTab_beta <- data_parameter_post_tab

xtab_beta <- xtable(dfTab_beta)
colnames(xtab_beta) <- c(
  "",
  "Age",
  "Treat.",
  "Sex",
  "Perf.",
  "Age",
  "Treat.",
  "Sex",
  "Perf.",
  "Age",
  "Treat.",
  "Sex",
  "Perf."
)

align(xtab_beta) <- rep("c", length(colnames(xtab_beta)) + 1)

## Add Multicolumns
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\toprule & & & \\multicolumn{2}{c}{Tests about the effect size $\\theta$} & \\multicolumn{2}{c}{Tests about the power parameter $\\alpha$} \\\\ \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}'
print(xtab_beta, floating = FALSE, include.rownames = FALSE, add.to.row = addtorow,
      sanitize.text.function = function(x){x}, booktabs = TRUE, hline.after = c(0, nrow(xtab_beta)))


#   ____________________________________________________________________________
#   Plots                                                                   ####


# Age Plots
long_data_age <- gather(data_parameter_post, key = "variable", value = "value", 
                        age_beta_opt, age_beta_unif)

plot_age <- ggplot(long_data_age, aes(x = value, fill = variable)) +
  geom_density(data = subset(long_data_age, variable == "age_beta_unif"), 
               aes(x = value, fill = variable), alpha = 0.95) +
  geom_density(data = subset(long_data_age, variable == "age_beta_opt"), 
               aes(x = value, fill = variable), alpha = 0.8) +
  labs(title = "Age", x = "", y = "Posterior Density") +
  scale_fill_manual(values = c("#8A0910", "#0072B2"), 
                    name = "Prior",
                    labels = c("CBF", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_age_opt[1], xmax = hpdi_age_opt[2], y = 41* 1.15, height = 2),
                 color = "#8A0910", alpha = 1, size = 1.5) +
  geom_errorbarh(aes(xmin = hpdi_age_unif[1], xmax = hpdi_age_unif[2], y = 41 * 1.05, height = 2),
                 color = "#0072B2", alpha = 1, size = 1.5) +

  theme_light(base_size = 16)+
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
  scale_fill_manual(values = c("#8A0910", "#0072B2"), 
                    name = "Prior",
                    labels = c("CBF", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_sex_opt[1], xmax = hpdi_sex_opt[2], y = 2.0 * 1.15, height = 0.15), 
                 color = "#8A0910", alpha = 1, size = 1.5) +
  geom_errorbarh(aes(xmin = hpdi_sex_unif[1], xmax = hpdi_sex_unif[2], y = 2.0 * 1.05, height = 0.15), 
                 color = "#0072B2", alpha = 1, size = 1.5) +
  theme_light(base_size = 16) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.8)) +  # Adjust y-axis
  scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5), limits = c(-1.2, 1.2)) +  # Adjust x-axis
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))

# Print the Plot
print(plot_sex)


# Performance Plots
long_data_perform <- gather(data_parameter_post, key = "variable", value = "value", 
                        perform_beta_opt, perform_beta_unif)

plot_perform <- ggplot(long_data_perform, aes(x = value, fill = variable)) +
  geom_density(data = subset(long_data_perform, variable == "perform_beta_unif"), 
               aes(x = value, fill = variable), alpha = 0.95) +
  geom_density(data = subset(long_data_perform, variable == "perform_beta_opt"), 
               aes(x = value, fill = variable), alpha = 0.8) +
  labs(title = "Performance Status", x = "Value", y = "") +
  scale_fill_manual(values = c("#8A0910", "#0072B2"), 
                    name = "Prior",
                    labels = c("CBF", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_perform_opt[1], xmax = hpdi_perform_opt[2], y = 1.5 * 1.15, height = 0.10), 
                 color = "#8A0910", alpha = 1, size = 1.5) +
  geom_errorbarh(aes(xmin = hpdi_perform_unif[1], xmax = hpdi_perform_unif[2], y = 1.5 * 1.05, height = 0.10), 
                 color = "#0072B2", alpha = 1, size = 1.5) +
  theme_light(base_size = 16) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.0)) +  # Adjust y-axis
  scale_x_continuous(breaks = seq(-2, 2, by = 0.5), limits = c(-2, 1)) +  # Adjust x-axis
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))

# Print the Plot
print(plot_perform)


# Treatment Plots
long_data_treat <- gather(data_parameter_post, key = "variable", value = "value", 
                          treat_beta_opt, treat_beta_unif)

plot_treat <- ggplot(long_data_treat, aes(x = value, fill = variable)) +
  geom_density(data = subset(long_data_treat, variable == "treat_beta_unif"), 
               aes(x = value, fill = variable), alpha = 0.95) +
  geom_density(data = subset(long_data_treat, variable == "treat_beta_opt"), 
               aes(x = value, fill = variable), alpha = 0.8) +
  labs(title = "Treatment", x = "Value", y = "") +
  scale_fill_manual(values = c("#8A0910", "#0072B2"), 
                    name = "Prior",
                    labels = c("CBF", "Default")) +
  geom_errorbarh(aes(xmin = hpdi_treat_opt[1], xmax = hpdi_treat_opt[2], y = 1.8 * 1.15, height = 0.10), 
                 color = "#8A0910", alpha = 1, size = 1.5) +
  geom_errorbarh(aes(xmin = hpdi_treat_unif[1], xmax = hpdi_treat_unif[2], y = 1.8 * 1.05, height = 0.10), 
                 color = "#0072B2", alpha = 1, size = 1.5) +
  theme_light(base_size = 16) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.2)) +  # Adjust y-axis
  # scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5), limits = c(-1.2, 1.2)) +  # Adjust x-axis
  theme(plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))

# Print the Plot
print(plot_treat)


# 
# # Delta plot
# long_data_delta <- gather(data_parameter_post, key = "variable", value = "value", 
#                           delta_opt, delta_unif)
# hpdi_delta_opt <- rethinking::HPDI( data_parameter_post$delta_opt, prob = 0.95)
# hpdi_delta_unif <- rethinking::HPDI( data_parameter_post$delta_unif, prob = 0.95)
# plot_delta <- ggplot(long_data_delta, aes(x = value, fill = variable)) +
#   geom_density(data = subset(long_data_delta, variable == "delta_unif"), 
#                aes(x = value, fill = variable), alpha = 0.95) +
#   geom_density(data = subset(long_data_delta, variable == "delta_opt"), 
#                aes(x = value, fill = variable), alpha = 0.8) +
#   labs(x = expression("Values of "* delta), y = "Posterior Density") +
#   scale_fill_manual(values = c("#D55E00", "grey20"), 
#                     name = "Prior: ",
#                     labels = c("Optimal", "Default")) +
#   geom_errorbarh(aes(xmin = hpdi_delta_opt[1], xmax = hpdi_delta_opt[2], y = 15 * 1.15, height = 0.15), 
#                  color = "#D55E00", alpha = 1, size = 0.5) +
#   geom_errorbarh(aes(xmin = hpdi_delta_unif[1], xmax = hpdi_delta_unif[2], y = 15 * 1.05, height = 0.15), 
#                  color = "grey20", alpha = 1, size = 0.5) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.position = "top")
# 
# # Print the plot
# print(plot_delta)



#   ____________________________________________________________________________
#   Plots with facet_wrap                                                   ####


data_combined <- bind_rows(
  mutate(long_data_age, category = "Age"),
  mutate(long_data_sex, category = "Sex"),
  mutate(long_data_perform, category = "Performance"),
  mutate(long_data_treat, category = "Treatment")
)


# Ensure 'variable' is a factor and its levels match the ones used in scale_fill_manual
data_combined$variable <- factor(data_combined$variable, levels = c("age_beta_opt", "age_beta_unif", "sex_beta_opt", "sex_beta_unif", "perform_beta_opt", "perform_beta_unif", "treat_beta_opt", "treat_beta_unif"))
data_combined$prior_type <- ifelse(grepl("opt", data_combined$variable), "CBF", "Uniform")


error_bars_df <- data.frame(
  category = rep(c("Age", "Sex", "Performance", "Treatment"), each = 2),
  xmin = c(hpdi_age_opt[1], hpdi_age_unif[1], hpdi_sex_opt[1], hpdi_sex_unif[1], 
           hpdi_perform_opt[1], hpdi_perform_unif[1], hpdi_treat_opt[1], hpdi_treat_unif[1]),
  xmax = c(hpdi_age_opt[2], hpdi_age_unif[2], hpdi_sex_opt[2], hpdi_sex_unif[2], 
           hpdi_perform_opt[2], hpdi_perform_unif[2], hpdi_treat_opt[2], hpdi_treat_unif[2]),
  y = c(41 * 1.15, 41 * 1.05, 2.0 * 1.15, 2.0 * 1.05, 
        1.5 * 1.15, 1.5 * 1.05, 1.8 * 1.15, 1.8 * 1.05),
  height = c(2, 2, 0.15, 0.15, 0.10, 0.10, 0.10, 0.10),
  color = rep(c("#8A0910", "#0072B2"), 4),
  alpha = rep(1, 8),
  size = rep(1.5, 8)
)

# Custom names for each category
facet_names <- c(Age = "Age", 
                 Sex = "Sex", 
                 Performance = "Performance Status", 
                 Treatment = "Treatment")


# Plotting "Uniform" densities first ensures they are in the background
plot_combined <- ggplot() +
  geom_density(data = filter(data_combined, prior_type == "Uniform"), 
               aes(x = value, fill = prior_type), alpha = 0.75) +
  geom_density(data = filter(data_combined, prior_type == "CBF"), 
               aes(x = value, fill = prior_type), alpha = 0.75) +
  scale_fill_manual(values = c("CBF" = "#8A0910", "Uniform" = "#0072B2")) +
  facet_wrap(~category, scales = "free", labeller = as_labeller(facet_names)) +
  geom_errorbarh(data = error_bars_df, aes(xmin = xmin, xmax = xmax, y = y, height = height), 
                 inherit.aes = FALSE, color = error_bars_df$color, alpha = error_bars_df$alpha, size = error_bars_df$size) +
  theme_bw(base_size = 18) +
  theme(legend.position = "top",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  labs(y = "Density", x = "Values") +
  guides(fill = guide_legend(title = expression(bold(paste("Prior for ", delta, ":")))))
# Print the combined plot with specified layering
print(plot_combined)


ggsave(filename = "post_param_melanoma.pdf",path = "Plots", plot = plot_combined,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)


