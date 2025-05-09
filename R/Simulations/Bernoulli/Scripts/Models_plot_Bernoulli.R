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
  "ggpubr",
  "grid",
  "modi")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

source("R/CBF_Functions.R")

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
  
  return(bayes_factor)
}

#   ____________________________________________________________________________
#   Loading all the quantities of interest                                  ####

# Function to process the datasets in a specific directory
process_datasets <- function() {
  
  # Define the path to the directory
  path <- "R/Simulations/Bernoulli/RData"
  
  # List all RData files in the specified directory
  files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  files <- gtools::mixedsort(files)
  hpdi <- c(0.75) # HPDI intervals

  results <- list()
  delta_post <- list()
  
  for (file in files) {

    load(file)
    
    # Reshape the data into long format
    bf_data_long <- bf_data %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")
    bf_data_long <- as.data.frame(bf_data_long)
    
    
    for (j in 1:length(hpdi)) {
      
      # Optimal CBF prior on delta
      optimal_model_index <- calculate_maximized_value(bf_data_long$Value, grid_beta$obs_bf, hpdi[j])
      
      norm_const_unif <- marg_lik(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                  p = q, q = q, eta = 1, nu = 1)
      post_theta_unif <- marg_post_theta(theta = theta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                         p = q, q = q, 
                                         eta = 1, nu = 1, norm_const)
      post_delta_unif <- marg_post_delta(delta = delta_seq[2:999], y = y, N = N, y_0 = y_0, N_0 = N_0,
                                         p = q, q = q, 
                                         eta = 1, nu = 1, norm_const)
      norm_const_alt <- marg_lik(delta = delta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                 p = q, q = q, eta = grid_beta$eta_val[optimal_model_index], nu = grid_beta$nu_val[optimal_model_index])
      post_theta_alt <- marg_post_theta(theta = theta_seq, y = y, N = N, y_0 = y_0, N_0 = N_0,
                                        p = q, q = q, 
                                        eta = grid_beta$eta_val[optimal_model_index], nu = grid_beta$nu_val[optimal_model_index], norm_const)
      post_delta_alt <- marg_post_delta(delta = delta_seq[2:999], y = y, N = N, y_0 = y_0, N_0 = N_0,
                                         p = q, q = q, 
                                        eta = grid_beta$eta_val[optimal_model_index], nu = grid_beta$nu_val[optimal_model_index], norm_const)
      
      sd_theta_unif <- sqrt(weighted.var(theta_seq, post_theta_unif))
      sd_theta_alt <- sqrt(weighted.var(theta_seq, post_theta_alt))
      sd_delta_unif <- sqrt(weighted.var(delta_seq[2:999], post_delta_unif))
      sd_delta_alt <- sqrt(weighted.var(delta_seq[2:999], post_delta_alt))
      
      mean_theta_unif <- weighted.mean(theta_seq, post_theta_unif)
      mean_theta_alt <- weighted.mean(theta_seq, post_theta_alt)
      mean_delta_unif <- weighted.mean(delta_seq[2:999], post_delta_unif)
      mean_delta_alt <- weighted.mean(delta_seq[2:999], post_delta_alt)
      
      
      if(optimal_model_index == "The model under the Null Hypothesis is optimal"){
        p <- c(0.25, 0.5, 0.75)
        beta_quant <- qbeta(p, 1, 1)
        mean_beta <- 0.5
        
        # Store the result
        results[[file]][[j]] <- list(opt_index = optimal_model_index, eta = 1,
                                     nu = 1, theta_dif = (y-y_0),
                                     first_q = round(beta_quant[1],2), median = round(beta_quant[2],2),
                                     third_q = round(beta_quant[3],2), mean = round(mean_beta,2), hpdi_val = hpdi[j],
                                     sd_theta_ref = sd_theta_unif, sd_theta_alter = sd_theta_unif,
                                     mean_theta_ref = mean_theta_unif, mean_theta_alter = mean_theta_unif,
                                     sd_delta_ref = sd_delta_unif, sd_delta_alter = sd_delta_unif,
                                     mean_delta_ref = mean_delta_unif, mean_delta_alter = mean_delta_unif)
      } else{
        p <- c(0.25,0.5,0.75)
        beta_quant <- qbeta(p, grid_beta$eta_val[optimal_model_index], grid_beta$nu_val[optimal_model_index])
        mean_beta <- grid_beta$eta_val[optimal_model_index]/(grid_beta$eta_val[optimal_model_index]+grid_beta$nu_val[optimal_model_index])
        
        # Store the result
        results[[file]][[j]] <- list(opt_index = optimal_model_index, eta = grid_beta$eta_val[optimal_model_index],
                                     nu = grid_beta$nu_val[optimal_model_index], theta_dif = (y-y_0),
                                     first_q = round(beta_quant[1],2), median = round(beta_quant[2],2),
                                     third_q = round(beta_quant[3],2), mean = round(mean_beta,2), hpdi_val = hpdi[j],
                                     sd_theta_ref = sd_theta_unif, sd_theta_alter = sd_theta_alt,
                                     mean_theta_ref = mean_theta_unif, mean_theta_alter = mean_theta_alt,
                                     sd_delta_ref = sd_delta_unif, sd_delta_alter = sd_delta_alt,
                                     mean_delta_ref = mean_delta_unif, mean_delta_alter = mean_delta_alt)
      }
    }
  }
  return(results)
}

results_opt_beta <- process_datasets()

#   ____________________________________________________________________________
#   Convert to dataframe                                                    ####

# Function to flatten the nested list into a format suitable for a dataframe
flatten_results <- function(results_list) {
  flattened <- do.call(rbind, lapply(names(results_list), function(file_name) {
    do.call(rbind, lapply(results_list[[file_name]], function(hpdi_result) {
      cbind(data.frame(file = file_name), as.data.frame(t(unlist(hpdi_result))))
    }))
  }))
  return(flattened)
}

# Convert the nested list into a flattened dataframe
results_df <- flatten_results(results_opt_beta)

correct_col_names <- c("file", "opt_index", "eta", "nu", "theta_dif",
                       "first_q", "median", "third_q", "mean", "hpdi_val",
                       "sd_theta_ref", "sd_theta_alter",
                       "mean_theta_ref", "mean_theta_alter",
                       "sd_delta_ref", "sd_delta_alter",
                       "mean_delta_ref", "mean_delta_alter")

names(results_df) <- correct_col_names
results_df$hpdi_val <- factor(results_df$hpdi_val, levels = c("0.75"))
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

# Adding RMAP and SAM prior simulations
load("R/Simulations/Bernoulli/RData/RMAP_SAM_prior/RMAP_prior.RData")
load("R/Simulations/Bernoulli/RData/RMAP_SAM_prior/SAM_prior.RData")

print(results_df)
results_df <- cbind(results_df, data_mix, data_SAM)
results_df <- as_tibble(results_df)

#   ____________________________________________________________________________
#   Plots                                                                   ####


##  ............................................................................
##  CBF prior on delta                                                      ####

facet_titles <- setNames(c("75% HPDI"), levels(results_df$hpdi_val))

binomial_comp <- ggplot(results_df, aes(x = theta_dif, y = median, color = hpdi_val, group = hpdi_val)) +
  geom_point(size=2) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin = first_q, ymax = third_q), width = 0.46, size = 1.4) +
  geom_text(aes(y = third_q, label = paste("(", eta, ",", nu, ")")), vjust = -0.8, hjust = -0.1, angle = 60, size = 4.3, check_overlap = FALSE) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 22, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, angle = 0, vjust = 1)) + 
  labs(
    x = expression("Values of " * (y-y[0])),
    y = expression("Values of  " *  pi[0](delta)[CBF]),
    color = "HPDI"
  ) +
  scale_x_continuous(breaks = seq(0, as.integer(max(results_df$theta_dif)), 1), guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = scales::number_format(accuracy = 0.1), limits = c(0, max(results_df$third_q+0.15)), guide = guide_axis(check.overlap = TRUE))+
  scale_color_manual(values = c("#8A0910", "#0072B2", "#009E20")) # Custom colors


print(binomial_comp)


ggsave(filename = "binomial_comp.pdf",path = "Plots", plot = binomial_comp,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)


##  ............................................................................
##  Posterior SD of theta                                                   ####

df_long <- pivot_longer(results_df, cols = c(sd_theta_ref, sd_theta_alter, 
                                             sd_theta_mix, sd_theta_SAM), 
                        names_to = "sd_type", values_to = "value")

df_long <- df_long %>%
  mutate(
    method = case_when(
      sd_type == "sd_theta_ref"       ~ "Unif.",
      sd_type == "sd_theta_alter"     ~ "CBF",
      sd_type == "sd_theta_mix"       ~ "RMAP",
      sd_type == "sd_theta_SAM"       ~ "SAM"
    )
  )

comparisons <- list(
  c("CBF", "RMAP"),
  c("CBF", "Unif."),
  c("CBF", "SAM")
)

# Function to create facet labels using plotmath
create_facet_label <- function(method1, method2) {
  if (method2 == "Unif.") {
    return(paste0("NPP[", method1, "] ~ 'vs' ~ NPP[", method2, "]"))
  } else {
    return(paste0("NPP[", method1, "] ~ 'vs' ~ ", method2))
  }
}

df_compare_list <- lapply(comparisons, function(comp) {
  method1 <- comp[1]  
  method2 <- comp[2]  
  
  # Facet label using plotmath
  facet_label <- create_facet_label(method1, method2)
  
  df_long %>%
    filter(method %in% comp) %>%
    mutate(facet = facet_label)
})

df_compare <- bind_rows(df_compare_list)

cb_palette <- c(
  "CBF"      = "#8A0910",  
  "RMAP"     = "#E69F00", 
  "Unif."    = "#0072B2",  
  "SAM"      = "#009E73"   
)

cb_shapes <- c(
  "CBF"      = 19,
  "RMAP"     = 17,
  "Unif."    = 18,
  "SAM"      = 15
)

label_expressions <- c(
  "CBF"     = expression(NPP["CBF"]),
  "Unif."   = expression(NPP["Unif."]),
  "RMAP"    = "RMAP",
  "SAM"     = "SAM"
)

# Reorder the 'method' factor to control legend order
df_compare <- df_compare %>%
  mutate(method = factor(method, levels = c("CBF", "Unif.", "RMAP", "SAM")))

standard_dev_plot <- ggplot(df_compare, aes(x = theta_dif, y = value, color = method, shape = method, group = method)) +
    geom_point( size = 3.8, alpha = 0.8) +
    geom_line( size = 1.2) +
    scale_color_manual(
    values = cb_palette,
    labels = label_expressions,
    name = expression(bold("Prior:"))
  ) +
    scale_shape_manual(
    values = cb_shapes,
    labels = label_expressions,
    name = expression(bold("Prior:"))
  ) +
    theme_bw(base_size = 16) +
    labs(x = element_blank(),
       y = element_blank(),
       title = expression("SD of " * pi(theta ~ "|" ~ y ~ "," ~ y[0] ~ "," ~ delta)),
         color = expression(bold(paste("Prior:"))),
         shape = expression(bold(paste("Prior:")))) +
    theme(
    axis.text.x = element_text(size = 18, angle = 0, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
    scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 2),
                       guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~ facet, labeller = label_parsed, ncol = 1)  

print(standard_dev_plot)


ggsave(filename = "sd_binomial_comp.pdf",path = "Plots", plot = standard_dev_plot,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)


##  ............................................................................
##  Posterior mean of theta                                                 ####

df_long_mean <- pivot_longer(results_df, cols = c(mean_theta_ref, mean_theta_alter,
                                                  mean_theta_mix, mean_theta_SAM),
                             names_to = "mean_type", values_to = "value")

df_long_mean <- df_long_mean %>%
  mutate(
    method = case_when(
      mean_type == "mean_theta_ref"       ~ "Unif.",
      mean_type == "mean_theta_alter"     ~ "CBF",
      mean_type == "mean_theta_mix"       ~ "RMAP",
      mean_type == "mean_theta_SAM"       ~ "SAM"
    )
  )

comparisons <- list(
  c("CBF", "RMAP"),
  c("CBF", "Unif."),
  c("CBF", "SAM")
)

df_compare_list <- lapply(comparisons, function(comp) {
  method1 <- comp[1]  
  method2 <- comp[2]  
  
  facet_label <- create_facet_label(method1, method2)
  
  df_long_mean %>%
    filter(method %in% comp) %>%
    mutate(facet = facet_label)
})

df_compare <- bind_rows(df_compare_list)

cb_palette <- c(
  "CBF"      = "#8A0910",  
  "RMAP"     = "#E69F00", 
  "Unif."    = "#0072B2",  
  "SAM"      = "#009E73"   
)

cb_shapes <- c(
  "CBF"      = 19,  
  "RMAP"     = 17,  
  "Unif."    = 18,  
  "SAM"      = 15 
)

label_expressions <- c(
  "CBF"     = expression(NPP["CBF"]),
  "Unif." = expression(NPP["Unif."]),
  "RMAP"    = "RMAP",
  "SAM"    = "SAM"
)

df_compare <- df_compare %>%
  mutate(method = factor(method, levels = c("CBF", "Unif.", "RMAP", "SAM")))

mean_plot <- ggplot(df_compare, aes(x = theta_dif, y = value, color = method, shape = method, group = method)) +
    geom_point( size = 3.8, alpha = 0.8) +
    geom_line( size = 1.2) +
    scale_color_manual(
    values = cb_palette,
    labels = label_expressions,
    name = expression(bold("Prior:"))
  ) +
    scale_shape_manual(
    values = cb_shapes,
    labels = label_expressions,
    name = expression(bold("Prior:"))
  ) +
    theme_bw(base_size = 16) +
    labs(x = element_blank(),
       y = element_blank(),
       title = expression("Mean of " * pi(theta ~ "|" ~ y ~ "," ~ y[0] ~ "," ~ delta)),
       color = expression(bold(paste("Prior:"))),
       shape = expression(bold(paste("Prior:")))) +
    theme(
    axis.text.x = element_text(size = 18, angle = 0, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  ) +
    scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 2),
                     guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~ facet, labeller = label_parsed, ncol = 1)  # Arrange facets in a single column for better vertical alignment

print(mean_plot)

ggsave(filename = "mean_binomial_comp.pdf",path = "Plots", plot = mean_plot,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)



##  ............................................................................
##  Combined plot posterior for theta                                       ####


plot_comb_theta <- ggarrange(standard_dev_plot,
                             mean_plot,
                             ncol = 2,
                             nrow = 1,
                             common.legend = TRUE)

plot_comb_theta <- annotate_figure(plot_comb_theta,
                                   left = textGrob(
                                    "Posterior Values",
                                    rot = 90,
                                    vjust = 0.8,
                                    gp = gpar(cex = 1.8)
                                   ),
                                   bottom = textGrob(
                                    expression("Values of " * (y - y[0])),
                                    vjust = -0.1,
                                    gp = gpar(cex = 1.8)
                                   )
                                 )


ggsave(filename = "combined_plot_binomial.pdf",path = "Plots", plot = plot_comb_theta,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)


##  ............................................................................
##  Plot posterior SD and mean of delta                                     ####

# Combine the data for both plots
df_long_sd_delta <- pivot_longer(results_df, cols = c(sd_delta_ref, sd_delta_alter), names_to = "type", values_to = "value")
df_long_sd_delta <- df_long_sd_delta %>%
  mutate(plot_type = "Standard Deviation")

df_long_mean_delta <- pivot_longer(results_df, cols = c(mean_delta_ref, mean_delta_alter), names_to = "type", values_to = "value")
df_long_mean_delta <- df_long_mean_delta %>%
  mutate(plot_type = "Mean")

combined_df <- bind_rows(df_long_sd_delta, df_long_mean_delta)

# Add the new column 'prior'
combined_df <- combined_df %>%
  mutate(prior = case_when(
    type %in% c("sd_delta_alter", "mean_delta_alter") ~ "CBF",
    type %in% c("sd_delta_ref", "mean_delta_ref") ~ "Uniform" ))

combined_df <- combined_df %>%
  mutate(stat_label = case_when(
    plot_type == "Mean" ~ "Mean~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~theta)",
    plot_type == "Standard Deviation"   ~ "SD~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~theta)"
  ))

# Reorder the factor levels: SD first, then Mean
combined_df$stat_label <- factor(combined_df$stat_label, 
                                 levels = c(
                                   "SD~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~theta)",
                                   "Mean~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~theta)"
                                 ))

combined_df <- combined_df %>%
  mutate(
    method = case_when(
      prior == "Uniform" ~ "Unif.",
      prior == "CBF"  ~ "CBF")
  )

cb_palette <- c(
  "CBF"   = "#8A0910",  
  "RMAP"  = "#E69F00",  
  "Unif." = "#0072B2",  
  "SAM"  = "#009E73"   
)

cb_shapes <- c(
  "CBF"   = 19,  
  "RMAP"  = 17,  
  "Unif." = 18,  
  "SAM"  = 15   
)

label_expressions <- c(
  "CBF"     = expression(NPP["CBF"]),
  "Unif."   = expression(NPP["Unif."]),
  "RMAP"    = "RMAP",
  "SAM"    = "SAM"
)

binomial_delta_plot <- ggplot(combined_df, aes(x = theta_dif, y = value, color = method, shape = method)) +
  geom_point(size = 3.8, alpha = 0.8) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = cb_palette,
    labels = label_expressions,
    name = expression(bold("Prior:"))
  ) +
    scale_shape_manual(
      values = cb_shapes,
      labels = label_expressions,
      name = expression(bold("Prior:"))
    ) +
    theme_bw(base_size = 16) +
    labs(
      x = expression("Values of " * (y - y[0])),
      y = "Posterior Values"
    ) +
    theme(
      axis.text.x = element_text(size = 18, angle = 0, vjust = 1),
      plot.title = element_text(hjust = 0.5, size = 22),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(size = 18),
      legend.position = "top",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 18),
      strip.text = element_text(size = 16, face = "bold"),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(breaks = seq(0, max(results_df$theta_dif), 2), guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~stat_label, scales = "free_y", labeller = label_parsed)

  print(binomial_delta_plot)

  ggsave(
    filename = "binomial_delta_plot.pdf", path = "Plots", plot = binomial_delta_plot,
    width = 15, height = 8, device = "pdf", dpi = 500, useDingbats = FALSE
  )
