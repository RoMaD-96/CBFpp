#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "dplyr",
  "progressr",
  "ggplot2",
  "ggridges",
  "tidyverse",
  "tidyr",
  "ggpubr",
  "grid",
  "tibble",
  "modi")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))



source("R/CBF_Functions.R")

#   ____________________________________________________________________________
#   Load all the quantities of interest                                     ####


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

# Define the function to process the datasets in a specific directory
process_datasets <- function() {
  
  
  # Define the path to the directory
  path <- "R/Simulations/Gaussian/RData"
  
  # List all RData files in the specified directory
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  files <- gtools::mixedsort(files)
  hpdi <- c(0.75)
  # Initialize a list to store results
  results <- list()
  
  
  # Loop through each file
  for (file in files) {
    
    # Load the dataset
    load(file)
    
    # Assuming the data is loaded into 'bf_data'
    # Reshape the data into long format
    bf_data_long <- bf_data %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")
    bf_data_long <- as.data.frame(bf_data_long)
    
    for (j in 1:length(hpdi)) {
      optimal_model_index <- calculate_maximized_value(bf_data_long$Value, grid_beta$obs_bf, hpdi[j])
      post_theta_unif <- theta_m_post(theta_seq = theta_seq, theta = theta,
                                      theta_0 = theta_0, st_dev = se, st_dev_0 = se_0,
                                      eta = 1, nu = 1)
      post_delta_unif <- delta_m_post(delta_seq = delta_seq[2:499], theta = theta,
                                      theta_0 = theta_0, st_dev = se, st_dev_0 = se_0,
                                      eta = 1, nu = 1)
      sd_theta_unif <- sqrt(weighted.var(theta_seq, post_theta_unif))
      mean_theta_unif <- weighted.mean(theta_seq, post_theta_unif)
      mean_delta_unif <- weighted.mean(delta_seq[2:499], post_delta_unif)
      sd_delta_unif <- sqrt(weighted.var(delta_seq[2:499], post_delta_unif))
      
    
      post_theta_alt <- theta_m_post(theta_seq = theta_seq, theta = theta,
                                     theta_0 = theta_0, st_dev = se, st_dev_0 = se_0,
                                     eta = grid_beta$eta_val[optimal_model_index], nu = grid_beta$nu_val[optimal_model_index])
      post_delta_alt <- delta_m_post(delta_seq = delta_seq[2:499], theta = theta,
                                     theta_0 = theta_0, st_dev = se, st_dev_0 = se_0,
                                     eta = grid_beta$eta_val[optimal_model_index], nu = grid_beta$nu_val[optimal_model_index])
      sd_theta_alt <- sqrt(weighted.var(theta_seq, post_theta_alt))
      mean_theta_alt <- weighted.mean(theta_seq, post_theta_alt)
      mean_delta_alt <- weighted.mean(delta_seq[2:499], post_delta_alt)
      sd_delta_alt <- sqrt(weighted.var(delta_seq[2:499], post_delta_alt))
      
      
      
      if(optimal_model_index == "The model under the Null Hypothesis is optimal"){
        p <- c(0.25, 0.5, 0.75)
        beta_quant <- qbeta(p, 1, 1)
        mean_beta <- 0.5
        
        # Store the result
        results[[file]][[j]] <- list(opt_index = optimal_model_index, eta = 1,
                                     nu = 1, theta_dif = (true_theta-true_theta_0),
                                     first_q = round(beta_quant[1],2), median = round(beta_quant[2],2),
                                     third_q = round(beta_quant[3],2), mean = round(mean_beta,2), hpdi_val = hpdi[j],
                                     sd_theta_ref = sd_theta_unif, sd_theta_alter = sd_theta_unif,
                                     mean_theta_ref = mean_theta_unif, mean_theta_alter = mean_theta_unif,
                                     sd_delta_ref = sd_delta_unif, sd_delta_alter = sd_delta_unif,
                                     mean_delta_ref = mean_delta_unif, mean_delta_alter = mean_delta_unif
                                     )
      } else{
        p <- c(0.25,0.5,0.75)
        beta_quant <- qbeta(p, grid_beta$eta_val[optimal_model_index], grid_beta$nu_val[optimal_model_index])
        mean_beta <- grid_beta$eta_val[optimal_model_index]/(grid_beta$eta_val[optimal_model_index]+grid_beta$nu_val[optimal_model_index])

        # Store the result
        results[[file]][[j]] <- list(opt_index = optimal_model_index, eta = grid_beta$eta_val[optimal_model_index],
                                     nu = grid_beta$nu_val[optimal_model_index], theta_dif = (true_theta-true_theta_0),
                                     first_q = round(beta_quant[1],2), median = round(beta_quant[2],2),
                                     third_q = round(beta_quant[3],2), mean = round(mean_beta,2), hpdi_val = hpdi[j],
                                     sd_theta_ref = sd_theta_unif, sd_theta_alter = sd_theta_alt,
                                     mean_theta_ref = mean_theta_unif, mean_theta_alter = mean_theta_alt,
                                     sd_delta_ref = sd_delta_unif, sd_delta_alter = sd_delta_alt,
                                     mean_delta_ref = mean_delta_unif, mean_delta_alter = mean_delta_alt
        )
      }
        
      
    }
    # Calculate the optimal model index
    
  }
  return(results)
}

# To run the function

# Need to prespecify this
eta_null <- 1
nu_null <- 1
se <- 1
se_0 <- 1

results_opt_beta <- process_datasets()




#   ____________________________________________________________________________
#   Convert to dataframe                                                    ####


# Function to flatten the nested list into a format suitable for a DataFrame
flatten_results <- function(results_list) {
  flattened <- do.call(rbind, lapply(names(results_list), function(file_name) {
    do.call(rbind, lapply(results_list[[file_name]], function(hpdi_result) {
      cbind(data.frame(file = file_name), as.data.frame(t(unlist(hpdi_result))))
    }))
  }))
  return(flattened)
}

# Convert the nested list into a flattened DataFrame
results_df <- flatten_results(results_opt_beta)

# Ensure the column names are correct
correct_col_names <- c("file", "opt_index", "eta", "nu", "theta_dif",
                       "first_q", "median", "third_q", "mean", "hpdi_val",
                       "sd_theta_ref", "sd_theta_alter",
                       "mean_theta_ref", "mean_theta_alter",
                       "sd_delta_ref", "sd_delta_alter",
                       "mean_delta_ref", "mean_delta_alter")
names(results_df) <- correct_col_names
results_df$hpdi_val <- factor(results_df$hpdi_val, levels = c("0.75"))# Print the DataFrame
results_df$first_q <- as.numeric(results_df$first_q)
results_df$third_q <- as.numeric(results_df$third_q)
results_df$median <- as.numeric(results_df$median)
results_df$eta <- as.numeric(results_df$eta)
results_df$nu <- as.numeric(results_df$nu)
results_df$theta_dif <- as.numeric(results_df$theta_dif)
results_df$sd_theta_ref <- as.numeric(results_df$sd_theta_ref)
results_df$sd_theta_alter <- as.numeric(results_df$sd_theta_alter)
results_df$mean_theta_ref <- as.numeric(results_df$mean_theta_ref)
results_df$mean_theta_alter <- as.numeric(results_df$mean_theta_alter)
results_df$mean_delta_ref <- as.numeric(results_df$mean_delta_ref)
results_df$mean_delta_alter <- as.numeric(results_df$mean_delta_alter)
results_df$sd_delta_ref <- as.numeric(results_df$sd_delta_ref)
results_df$sd_delta_alter <- as.numeric(results_df$sd_delta_alter)

print(results_df)
results_df <- as_tibble(results_df)
#   ____________________________________________________________________________
#   Plots                                                                   ####

# Custom titles for facets
facet_titles <- setNames(c("75% HPDI"), levels(results_df$hpdi_val))
# Your ggplot code with corrected facet titles
gaussian_comp <- ggplot(results_df, aes(x = theta_dif, y = median, color = hpdi_val, group = hpdi_val)) +
  geom_point(size=2.5) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin = first_q, ymax = third_q), width = 0.1, size = 1.4) +
  geom_text(aes(y = third_q, label = paste("(",eta,",",nu,")")), vjust = -0.0, hjust = -0.2, size = 4.3, angle = 60, check_overlap = FALSE) +
  #facet_wrap(~ hpdi_val, ncol = 1, scales = "free_y", labeller = labeller(hpdi_val = facet_titles)) +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, angle = 30, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 22),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)) + # Adjust the angle and vertical alignment here
  labs(
    #title = expression("Gaussian with unknown mean " * mu),
    x = expression("Values of " * (mu[r]-mu[o])),
    y = expression("Values of  " *  pi[0](delta)[CBF]),
    color = "HPDI"
  ) +
  scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 0.2), guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = scales::number_format(accuracy = 0.1), limits = c(0, max(results_df$third_q+0.15)), guide = guide_axis(check.overlap = TRUE))+
  scale_color_manual(values = c("#8A0910","#0072B2",  "#009E20")) # Custom colors


print(gaussian_comp)


ggsave(filename = "gaussian_comp_mu_unknown.pdf",path = "Plots", plot = gaussian_comp,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)



df_long <- pivot_longer(results_df, cols = c(sd_theta_ref, sd_theta_alter), names_to = "sd_type", values_to = "value")

# Now create the plot
standard_dev_plot <- ggplot(df_long, aes(x = theta_dif, y = value, color = sd_type, shape = sd_type, group = sd_type)) +
  geom_point(position = position_dodge(width = 0), size = 3.8) +
  geom_line(position = position_dodge(width = 0), size = 1.5) +
  #facet_wrap(~ hpdi_val, scales = "free", labeller = labeller(hpdi_val = facet_titles)) +
  scale_color_manual(values = c("#8A0910", "#0072B2"), 
                     labels = c("CBF", "Uniform")) +
  scale_shape_manual(values = c(20, 18), 
                     labels = c("CBF", "Uniform")) +
  theme_bw(base_size = 16) +
  labs(x = expression("Values of " * (mu[r]-mu[o])), 
       y = expression("SD of " * pi(mu ~ "|" ~ y[0] ~ "," ~ y)), 
       color = expression(bold(paste("Prior for ", delta, ":"))), 
       shape = expression(bold(paste("Prior for ", delta, ":")))) +
  theme(
    axis.text.x = element_text(size = 18, angle = 30, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  ) +
  scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 0.2), 
                     guide = guide_axis(check.overlap = TRUE))

# Display the plot
print(standard_dev_plot)

ggsave(filename = "sd_gaussian_comp.pdf",path = "Plots", plot = standard_dev_plot,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)





df_long_mean <- pivot_longer(results_df, cols = c(mean_theta_ref, mean_theta_alter), names_to = "mean_type", values_to = "value")

# Now create the plot
mean_plot <- ggplot(df_long_mean, aes(x = theta_dif, y = value, color = mean_type, shape = mean_type, group = mean_type)) +
  geom_point(position = position_dodge(width = 0), size = 3.8) +
  geom_line(position = position_dodge(width = 0), size = 1.5) +
  #facet_wrap(~ hpdi_val, scales = "free", labeller = labeller(hpdi_val = facet_titles)) +
  scale_color_manual(values = c("#8A0910", "#0072B2"), 
                     labels = c("CBF", "Uniform")) +
  scale_shape_manual(values = c(20, 18), 
                     labels = c("CBF", "Uniform")) +
  theme_bw(base_size = 16) +
  labs(x = expression("Values of " * (mu[r]-mu[o])), 
       y = expression("Mean of " * pi(mu ~ "|" ~ y[0] ~ "," ~ y)), 
       color = expression(bold(paste("Prior for ", delta, ":"))), 
       shape = expression(bold(paste("Prior for ", delta, ":")))) +
  theme(
    axis.text.x = element_text(size = 18, angle = 30, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
  ) +
  scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 0.2), 
                     guide = guide_axis(check.overlap = TRUE))

# Display the plot
print(mean_plot)



#   ____________________________________________________________________________
#   ggarrange Combined Plot                                                 ####




# Combine the data for both plots
df_long_sd <- pivot_longer(results_df, cols = c(sd_theta_ref, sd_theta_alter), names_to = "type", values_to = "value")
df_long_sd <- df_long_sd %>%
  mutate(plot_type = "Standard Deviation")

df_long_mean <- pivot_longer(results_df, cols = c(mean_theta_ref, mean_theta_alter), names_to = "type", values_to = "value")
df_long_mean <- df_long_mean %>%
  mutate(plot_type = "Mean")

# Combine both data frames
combined_df <- bind_rows(df_long_sd, df_long_mean)

# Add the new column 'prior'
combined_df <- combined_df %>%
  mutate(prior = case_when(
    type %in% c("sd_theta_alter", "mean_theta_alter") ~ "CBF",
    type %in% c("sd_theta_ref", "mean_theta_ref") ~ "Uniform"
  ))

# # Define the levels for plot_type and convert it to a factor
# combined_df$plot_type <- factor(combined_df$plot_type, levels = c("Standard Deviation", "Mean"))
# 
# #Define facet labels
# facet_labels <- c(
#   "Standard Deviation" = expression("SD of " * pi(theta ~ "|" ~ y[0] ~ "," ~ y)),
#   "Mean" = expression("Mean of " * pi(theta ~ "|" ~ y[0] ~ "," ~ y))
# )
# 
# 
# levels(combined_df$plot_type) <- facet_labels

# Plot
combined_plot_theta <- ggplot(combined_df, aes(x = theta_dif, y = value, color = prior, shape = prior, group = interaction(type, plot_type))) +
  geom_point(position = position_dodge(width = 0), size = 3.8) +
  geom_line(position = position_dodge(width = 0), size = 1.5) +
  facet_wrap(~ plot_type, ncol = 1, scales = "free_y", labeller = labeller(plot_type = c("Standard Deviation" = "Standard Deviation", "Mean" = "Mean"))) +
  scale_color_manual(values = c("#8A0910", "#0072B2"), labels = c("CBF", "Uniform")) +
  scale_shape_manual(values = c(20, 18), labels = c("CBF", "Uniform")) +
  theme_bw(base_size = 16) +
  labs(x = element_blank(),
       y = element_blank(), 
       title= expression("Marginal Posterior for " * mu ),
       color = expression(bold("Prior for " * delta * ":")), 
       shape = expression(bold("Prior for " * delta * ":"))) +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, angle = 30, vjust = 1),
    strip.text.x = element_text(size = 18)
  ) +
  scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 0.4), guide = guide_axis(check.overlap = TRUE)) +
  guides(color = guide_legend(title.position = "top"), shape = guide_legend(title.position = "top"))

# Display the plot
print(combined_plot_theta)







# Combine the data for both plots
df_long_sd_delta <- pivot_longer(results_df, cols = c(sd_delta_ref, sd_delta_alter), names_to = "type", values_to = "value")
df_long_sd_delta <- df_long_sd_delta %>%
  mutate(plot_type = "Standard Deviation")

df_long_mean_delta <- pivot_longer(results_df, cols = c(mean_delta_ref, mean_delta_alter), names_to = "type", values_to = "value")
df_long_mean_delta <- df_long_mean_delta %>%
  mutate(plot_type = "Mean")

# Combine both data frames
combined_df <- bind_rows(df_long_sd_delta, df_long_mean_delta)

# Add the new column 'prior'
combined_df <- combined_df %>%
  mutate(prior = case_when(
    type %in% c("sd_delta_alter", "mean_delta_alter") ~ "CBF",
    type %in% c("sd_delta_ref", "mean_delta_ref") ~ "Uniform"
  ))

# # Define the levels for plot_type and convert it to a factor
# combined_df$plot_type <- factor(combined_df$plot_type, levels = c("Standard Deviation", "Mean"))
# 
# # Define facet labels
# facet_labels <- c(
#   "Standard Deviation" = expression("Standard Deviation of " * pi(delta ~ "|" ~ y[0] ~ "," ~ y)),
#   "Mean" = expression("Mean of " * pi(delta ~ "|" ~ y[0] ~ "," ~ y))
# )
# 
# levels(combined_df$plot_type) <- facet_labels

# Plot
combined_plot_delta <- ggplot(combined_df, aes(x = theta_dif, y = value, color = prior, shape = prior, group = interaction(type, plot_type))) +
  geom_point(position = position_dodge(width = 0), size = 3.8) +
  geom_line(position = position_dodge(width = 0), size = 1.5) +
  facet_wrap(~ plot_type, ncol = 1, scales = "free_y", labeller = labeller(plot_type = c("Standard Deviation" = "Standard Deviation", "Mean" = "Mean"))) +
  scale_color_manual(values = c("#8A0910", "#0072B2"), labels = c("CBF", "Uniform")) +
  scale_shape_manual(values = c(20, 18), labels = c("CBF", "Uniform")) +
  theme_bw(base_size = 16) +
  labs(x = element_blank(),
       y = element_blank(), 
       title= expression("Marginal Posterior for " * delta ),
       color = expression(bold("Prior for " * delta * ":")), 
       shape = expression(bold("Prior for " * delta * ":"))) +
  theme(
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, angle = 30, vjust = 1),
    strip.text.x = element_text(size = 18)
  ) +
  scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 0.4), guide = guide_axis(check.overlap = TRUE)) +
  guides(color = guide_legend(title.position = "top"), shape = guide_legend(title.position = "top"))

# Display the plot
print(combined_plot_delta)


plot_comb_delta_theta <- ggarrange(combined_plot_theta, combined_plot_delta, ncol = 2, nrow = 1, 
                                   common.legend = TRUE)
plot_comb_delta_theta <- annotate_figure(plot_comb_delta_theta, left = textGrob("Posterior Values", rot = 90, vjust = 0.8, gp = gpar(cex = 1.8)),
                                         bottom = textGrob(expression("Values of " * (mu[r]-mu[o])), vjust = -0.1,  gp = gpar(cex = 1.8)))


ggsave(filename = "combined_plot_gaussian.pdf",path = "Plots", plot = plot_comb_delta_theta,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)










plot_comb <- ggarrange(gaussian_comp, standard_dev_plot, ncol = 1, nrow = 2, 
                       common.legend = FALSE, heights = c(1,1.1))

ggsave(filename = "plot_comb_gaussian_comp.pdf",path = "Plots", plot = plot_comb,
       width = 15, height = 14, device='pdf', dpi=500, useDingbats = FALSE)

