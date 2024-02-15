packages <- c(
  "rstan",
  "bridgesampling",
  "dplyr",
  "progressr",
  "ggplot2",
  "ggridges",
  "tidyverse",
  "tidyr",
  "tibble",
  "modi"
)

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

stan_model <- list("R/Simulations/Poisson/STAN/Poisson_Norm.stan")


# Define the function to process the datasets in a specific directory
process_datasets <- function() {
  # Define the path to the directory
  path <- "R/Simulations/Poisson/RData"
  source("R/CBF_Functions.R")
  
  # List all RData files in the specified directory
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  files <- gtools::mixedsort(files)
  hpdi <- c(0.65,0.75,0.95)
  # Initialize a list to store results
  results <- list()
  
  # Loop through each file
  for (file in files) {
    # Load the dataset
    load(file)
    stan_model <- list("R/Simulations/Poisson/STAN/Poisson_Norm.stan")
    source("R/CBF_Functions.R")
    # Assuming the data is loaded into 'bf_data'
    # Reshape the data into long format
    bf_data_long <- bf_data %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")
    bf_data_long <- as.data.frame(bf_data_long)
    
    data_stan_unif <- list(
      N0 = N_0,
      y0 = y_0,
      alpha0 = alpha_0,
      beta0 = beta_0,
      eta = 1,
      nu = 1,
      N = length(y),
      y = as.array(y))
    
    poi_stan_unif <- unnormalized_posterior_sampl(list(data_stan_unif),
                                                  stan_model,
                                                  iterations = 2000,
                                                  ref_rate = 0,
                                                  info_model = TRUE)
    
    poi_par_unif <- rstan::extract(poi_stan_unif[[1]])
    sd_lambda_unif <- sd(poi_par_unif[["lambda"]])
    
    for (j in 1:length(hpdi)) {
      optimal_model_index <- calculate_maximized_value(bf_data_long$Value, df_obs_bf$obs_bf, hpdi[j])
      stan_model <- list("R/Simulations/Poisson/STAN/Poisson_Norm.stan")
      source("R/CBF_Functions.R")
      
      data_stan <- list(
        N0 = N_0,
        y0 = y_0,
        alpha0 = alpha_0,
        beta0 = beta_0,
        eta = grid$x[optimal_model_index],
        nu = grid$y[optimal_model_index],
        N = length(y),
        y = as.array(y))
      
      poi_stan <- unnormalized_posterior_sampl(list(data_stan),
                                               stan_model,
                                               iterations = 2000,
                                               ref_rate = 0,
                                               info_model = TRUE)
      
      poi_par <- rstan::extract(poi_stan[[1]])
      sd_lambda_alt <- sd(poi_par[["lambda"]])
      
      p <- c(0.25,0.5,0.75)
      beta_quant <- qbeta(p, grid$x[optimal_model_index], grid$y[optimal_model_index])
      mean_beta <- grid$x[optimal_model_index]/(grid$x[optimal_model_index]+grid$y[optimal_model_index])
      
      # Store the result
      results[[file]][[j]] <- list(opt_index = optimal_model_index, eta = grid$x[optimal_model_index],
                                   nu = grid$y[optimal_model_index], lambda_dif = (true.lambda-true.lambda_0),
                                   first_q = beta_quant[1], median = beta_quant[2],
                                   third_q = beta_quant[3], mean = mean_beta, hpdi_val = hpdi[j],
                                   sd_lambda_ref = sd_lambda_unif, sd_lambda_alter = sd_lambda_alt)
    }
    # Calculate the optimal model index
    
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
correct_col_names <- c("file", "opt_index", "eta", "nu", "lambda_dif", "first_q", "median", "third_q", "mean", "hpdi_val", "sd_lambda_ref", "sd_lambda_alter" )
names(results_df) <- correct_col_names
results_df$hpdi_val <- factor(results_df$hpdi_val, levels = c("0.65", "0.75", "0.95"))# Print the DataFrame
print(results_df)


results_df <- results_df %>% 
  group_by(file, opt_index) %>%
  # Compute the minimum sd for each combination of model and opt_index
  mutate(min_sd = min(sd_lambda_alter)) %>%
  ungroup() %>%
  # Replace the sd with the calculated minimum sd 
  mutate(sd_lambda_alter = min_sd) %>%
  select(-min_sd)  

#   ____________________________________________________________________________
#   Plots                                                                   ####

# Custom titles for facets
facet_titles <- setNames(c("65% HPDI", "75% HPDI", "95% HPDI"), levels(results_df$hpdi_val))
# Your ggplot code with corrected facet titles
poisson_comp <- ggplot(results_df, aes(x = lambda_dif, y = median, color = hpdi_val, group = hpdi_val)) +
  geom_point(size=3) + # Adds points with specified size
  geom_line(size=1) + # Connects the points with lines
  geom_errorbar(aes(ymin = first_q, ymax = third_q), width = 0.01, size = 1.2) + # Adds error bars
  geom_text(aes(y = third_q, label = paste("(", eta, ",", nu, ")")), vjust = -0.8, hjust = 0.5, size = 3.5, check_overlap = TRUE) + # Adds text labels above error bars
  facet_wrap(~ hpdi_val, ncol = 1, scales = "free_y", labeller = labeller(hpdi_val = facet_titles)) + # Corrected facet labels
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 0, vjust = 1)) +
  labs(
    x = expression("Values of " * (lambda[c]-lambda[0])),
    y = expression("Median optimal " * delta),
    color = "HPDI"
  ) +
  scale_x_continuous(breaks = seq(0, max(results_df$lambda_dif), 0.05), labels = scales::number_format(accuracy = 0.05)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = scales::number_format(accuracy = 0.1), limits = c(0, max(results_df$third_q+0.1)))+
  scale_color_manual(values = c("#8A0404", "#0072B2", "#009E20")) # Custom colors

print(poisson_comp)


ggsave(filename = "poisson_comp.pdf",path = "Plots", plot = poisson_comp,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)

df_long <- pivot_longer(results_df, cols = c(sd_lambda_ref, sd_lambda_alter), names_to = "sd_type", values_to = "value")

# Now create the plot
standard_dev_plot <- ggplot(df_long, aes(x = lambda_dif, y = value, color = sd_type, shape = sd_type, group = sd_type)) +
  geom_point(position = position_dodge(width = 0), size = 3.8) +
  geom_line(position = position_dodge(width = 0), size = 1.5) +
  facet_wrap(~ hpdi_val, scales = "free", labeller = labeller(hpdi_val = facet_titles)) +
  scale_color_manual(values = c("#D55E00","grey20"), 
                     labels = c("Optimal", "Default")) +
  scale_shape_manual(values = c(20, 18), 
                     labels = c("Optimal", "Default")) +
  theme_minimal(base_size = 16) +
  labs(x = expression("Values of " * (lambda[c]-lambda[0])), 
       y = expression("SD of " * pi(lambda ~ "|" ~ y[0] ~ "," ~ y)), 
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
  scale_x_continuous(breaks = seq(0, max(results_df$lambda_dif), 0.05), 
                     guide = guide_axis(check.overlap = TRUE))

# Display the plot
print(standard_dev_plot)

ggsave(filename = "sd_poisson_comp.pdf",path = "Plots", plot = standard_dev_plot,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)


plot_comb <- ggarrange(poisson_comp, standard_dev_plot, ncol = 1, nrow = 2, 
                       common.legend = FALSE, heights = c(1.5,1))

ggsave(filename = "plot_comb_poisson_comp.pdf",path = "Plots", plot = plot_comb,
       width = 15, height = 13, device='pdf', dpi=500, useDingbats = FALSE)
