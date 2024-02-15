#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "dplyr",
  "progressr",
  "ggplot2",
  "ggridges",
  "tidyverse",
  "tidyr",
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

# Define the function to process the datasets in a specific directory
process_datasets <- function() {
  # Define the path to the directory
  path <- "R/Simulations/Bernoulli/RData"
  
  # List all RData files in the specified directory
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  files <- gtools::mixedsort(files)
  hpdi <- c(0.65, 0.75, 0.95)
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
      
      
      norm_const_unif <- marg_lik(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                  p = q, q = q, eta = 1, nu = 1)
      post_theta_unif <- marg_post_theta(theta = theta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                         p = q, q = q, 
                                         eta = 1, nu = 1, norm_const)
      norm_const_alt <- marg_lik(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                 p = q, q = q, eta = grid_beta$eta_val[optimal_model_index], nu = grid_beta$nu_val[optimal_model_index])
      post_theta_alt <- marg_post_theta(theta = theta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                        p = q, q = q, 
                                        eta = grid_beta$eta_val[optimal_model_index], nu = grid_beta$nu_val[optimal_model_index], norm_const)
      
      sd_theta_unif <- sqrt(weighted.var(theta_seq, post_theta_unif))
      sd_theta_alt <- sqrt(weighted.var(theta_seq, post_theta_alt))
      
      
      if(optimal_model_index == "The model under the Null Hypothesis is optimal"){
        p <- c(0.25, 0.5, 0.75)
        beta_quant <- qbeta(p, 1, 1)
        mean_beta <- 0.5
        
        # Store the result
        results[[file]][[j]] <- list(opt_index = optimal_model_index, eta = 1,
                                     nu = 1, theta_dif = (y-y_0),
                                     first_q = round(beta_quant[1],2), median = round(beta_quant[2],2),
                                     third_q = round(beta_quant[3],2), mean = round(mean_beta,2), hpdi_val = hpdi[j],
                                     sd_theta_ref = sd_theta_unif, sd_theta_alter = sd_theta_unif)
      } else{
        p <- c(0.25,0.5,0.75)
        beta_quant <- qbeta(p, grid_beta$eta_val[optimal_model_index], grid_beta$nu_val[optimal_model_index])
        mean_beta <- grid_beta$eta_val[optimal_model_index]/(grid_beta$eta_val[optimal_model_index]+grid_beta$nu_val[optimal_model_index])
        
        # Store the result
        results[[file]][[j]] <- list(opt_index = optimal_model_index, eta = grid_beta$eta_val[optimal_model_index],
                                     nu = grid_beta$nu_val[optimal_model_index], theta_dif = (y-y_0),
                                     first_q = round(beta_quant[1],2), median = round(beta_quant[2],2),
                                     third_q = round(beta_quant[3],2), mean = round(mean_beta,2), hpdi_val = hpdi[j],
                                     sd_theta_ref = sd_theta_unif, sd_theta_alter = sd_theta_alt)
      }
        
      
    }
  }
  return(results)
}

# To run the function
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
correct_col_names <- c("file", "opt_index", "eta", "nu", "theta_dif", "first_q", "median", "third_q", "mean", "hpdi_val", "sd_theta_ref", "sd_theta_alter")
names(results_df) <- correct_col_names
results_df$hpdi_val <- factor(results_df$hpdi_val, levels = c("0.65", "0.75", "0.95"))# Print the DataFrame
results_df$first_q <- as.numeric(results_df$first_q)
results_df$third_q <- as.numeric(results_df$third_q)
results_df$median <- as.numeric(results_df$median)
results_df$eta <- as.numeric(results_df$eta)
results_df$nu <- as.numeric(results_df$nu)
results_df$theta_dif <- as.numeric(results_df$theta_dif)
results_df$sd_theta_ref <- as.numeric(results_df$sd_theta_ref)
results_df$sd_theta_alter <- as.numeric(results_df$sd_theta_alter)

print(results_df)
results_df <- as_tibble(results_df)
#   ____________________________________________________________________________
#   Plots                                                                   ####

# Custom titles for facets
facet_titles <- setNames(c("65% HPDI", "75% HPDI", "95% HPDI"), levels(results_df$hpdi_val))
# Your ggplot code with corrected facet titles
binomial_comp <- ggplot(results_df, aes(x = theta_dif, y = median, color = hpdi_val, group = hpdi_val)) +
  geom_point(size=2) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin = first_q, ymax = third_q), width = 0.27, size = 0.9) +
  geom_text(aes(y = third_q, label = paste("(", eta, ",", nu, ")")), vjust = -0.8, hjust = -0.1, angle = 45, size = 3.5, check_overlap = FALSE) +
  facet_wrap(~ hpdi_val, ncol = 1, scales = "free_y", labeller = labeller(hpdi_val = facet_titles)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 0, vjust = 1)) + # Adjust the angle and vertical alignment here
  labs(
    #title = expression("Gaussian with unknown mean " * mu),
    x = expression("Values of " * (y-y[0])),
    y = expression("Median optimal " * delta),
    color = "HPDI"
  ) +
  scale_x_continuous(breaks = seq(0, as.integer(max(results_df$theta_dif)), 1), guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = scales::number_format(accuracy = 0.1), limits = c(0, max(results_df$third_q+0.38)), guide = guide_axis(check.overlap = TRUE))+
  scale_color_manual(values = c("#8A0404", "#0072B2", "#009E20")) # Custom colors


print(binomial_comp)


ggsave(filename = "binomial_comp.pdf",path = "Calibrated_BF/Plots", plot = binomial_comp,
       width = 16, height = 8, device='pdf', dpi=500, useDingbats = FALSE)


df_long <- pivot_longer(results_df, cols = c(sd_theta_ref, sd_theta_alter), names_to = "sd_type", values_to = "value")
# Now create the plot
standard_dev_plot <- ggplot(df_long, aes(x = theta_dif, y = value, color = sd_type, shape = sd_type, group = sd_type)) +
  geom_point(position = position_dodge(width = 0), size = 3.8) +
  geom_line(position = position_dodge(width = 0), size = 1.5) +
  facet_wrap(~ hpdi_val, scales = "free", labeller = labeller(hpdi_val = facet_titles)) +
  scale_color_manual(values = c("#D55E00","grey20"), 
                     labels = c("Optimal", "Default")) +
  scale_shape_manual(values = c(20, 18), 
                     labels = c("Optimal", "Default")) +
  theme_minimal(base_size = 16) +
  labs(x = expression("Values of " * (y-y[0])), 
       y = expression("SD of " * pi(theta ~ "|" ~ y[0] ~ "," ~ y)), 
       color = "Prior:", 
       shape = "Prior:") +
  theme(
    axis.text.x = element_text(size = 12, angle = 40, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 16),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 2), 
                     guide = guide_axis(check.overlap = TRUE))

# Display the plot
print(standard_dev_plot)

ggsave(filename = "sd_binomial_comp.pdf",path = "Plots", plot = standard_dev_plot,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)



plot_comb <- ggarrange(binomial_comp, standard_dev_plot, ncol = 1, nrow = 2, 
                       common.legend = FALSE, heights = c(1.5,1))

ggsave(filename = "plot_comb_binomial_comp.pdf",path = "Plots", plot = plot_comb,
       width = 15, height = 14, device='pdf', dpi=500, useDingbats = FALSE)

