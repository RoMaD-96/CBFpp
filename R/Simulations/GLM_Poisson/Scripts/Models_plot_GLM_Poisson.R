#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "dplyr",
  "ggplot2",
  "tidyverse",
  "tidyr",
  "ggpubr",
  "grid",
  "tibble",
  "repmix",
  "modi")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

source("R/CBF_Functions.R")
load("R/Simulations/GLM_Poisson/RData/GLM_Poisson_Sim_Models.RData")
load("R/Simulations/GLM_Poisson/RData/GLM_Poisson_Sim.RData")


#   ____________________________________________________________________________
#   Optimal CBF prior                                                       ####

cbf_prior <- data.frame()

p <- c(0.25,0.5,0.75)

for (model in seq_along(results_list)) {
  beta_quant <- qbeta(p,
                      results_list[[model]]$opt_bf$eta_H1,
                      results_list[[model]]$opt_bf$nu_H1)
  mean_beta <- with(results_list[[model]]$opt_bf, eta_H1 / (eta_H1 + nu_H1))

  row <- data.frame(
    q25        = beta_quant[1],
    q50        = beta_quant[2],
    q75        = beta_quant[3],
    eta        = results_list[[model]]$opt_bf$eta_H1,
    nu         = results_list[[model]]$opt_bf$nu_H1,
    hpdi       = factor("0.75", levels = c("0.75")),
    mean_beta  = mean_beta
  )

  cbf_prior <- rbind(cbf_prior, row)
}


# Different scenarios
beta1_curr <- c(0.50, 0.53, 0.54, 0.55, 0.57, 0.58, 0.59,
                0.60, 0.61, 0.62, 0.63, 0.64, 0.65)

diff_beta <- beta1_curr - 0.5

cbf_prior <- cbind(cbf_prior, beta1_curr, diff_beta)

#   ____________________________________________________________________________
#   Plots                                                                   ####


##  ............................................................................
##  CBF prior on delta                                                      ####

facet_titles <- setNames(c("75% HPDI"), levels(cbf_prior$hpdi))

GLM_Poisson_comp <- ggplot(cbf_prior, aes(x = diff_beta, y = q50, color = hpdi, group = hpdi)) +
  geom_point(size=2.5) +
  geom_line(size=1) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.006, size = 1.5) +
  geom_text(aes(y = q75, label = paste("(",eta,",",nu,")")), vjust = -0.0, hjust = -0.2, size = 5.3, angle = 60, check_overlap = FALSE) +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, angle = 0, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 22),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)) +
  labs(
    x = expression("Values of " * (beta["c"][","][1] - beta[0][","][1])),
    y = expression("Values of  " *  pi[0](delta)[CBF]),
    color = "HPDI"
  ) +
  scale_x_continuous(breaks = seq(0, max(cbf_prior$diff_beta), 0.01), guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = scales::number_format(accuracy = 0.1), limits = c(0, max(cbf_prior$q75+0.15)), guide = guide_axis(check.overlap = TRUE))+
  scale_color_manual(values = c("#8A0910","#0072B2",  "#009E20")) # Custom colors


print(GLM_Poisson_comp)

ggsave(filename = "GLM_Poisson_comp.pdf",path = "Plots", plot = GLM_Poisson_comp,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)


##  ............................................................................
##  Posterior SD of beta 1                                                  ####

sd_data <- data.frame()

for (model in 1:length(cbf_models)) {
  
  cbf_sd <- sd(cbf_models[[model]]$x1)
  unif_sd <- sd(unif_models[[model]]$x1)
  rmap_sd <- sd(rmap_models[[model]]$post.samples$x1)
  samp_sd <- sd(samp_models[[model]]$post.samples$x1)
  
  sd_vec <- c(cbf_sd, unif_sd, rmap_sd, samp_sd)
  
  sd_data <- rbind.data.frame(sd_data, sd_vec)
  
}


beta1_curr <- c(0.50, 0.53, 0.54, 0.55, 0.57, 0.58, 0.59,
                       0.60, 0.61, 0.62, 0.63, 0.64, 0.65)

diff_beta <- beta1_curr - 0.5

sd_data <- cbind(sd_data, beta1_curr, diff_beta)

colnames(sd_data) <- c("NPP_cbf", "NPP_unif", "RMAP", "SAM", "beta1_curr", "diff_beta")


df_long <- pivot_longer(sd_data, cols = c(NPP_cbf, NPP_unif,
                                             RMAP, SAM),
                        names_to = "sd_type", values_to = "value")

df_long <- df_long %>%
  mutate(
    method = case_when(
      sd_type == "NPP_unif" ~ "Unif.",
      sd_type == "NPP_cbf"  ~ "CBF",
      sd_type == "RMAP"     ~ "RMAP",
      sd_type == "SAM"     ~ "SAM"
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
  method1 <- comp[1]  # "CBF"
  method2 <- comp[2]  # e.g., "RMAP", "Unif.", "SAM"
  
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
  "Unif." = expression(NPP["Unif."]),
  "RMAP"    = "RMAP",
  "SAM"    = "SAM"
)

# Reorder to control legend order
df_compare <- df_compare %>%
  mutate(method = factor(method, levels = c("CBF", "Unif.", "RMAP", "SAM")))

standard_dev_plot <- ggplot(df_compare, aes(x = diff_beta, y = value, color = method, shape = method, group = method)) +
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
    labs(
    x = element_blank(),
    y = element_blank(),
    title = expression("SD of " * pi(beta[1] ~ "|" ~ y ~ "," ~ y[0] ~ "," ~ X ~ ","~ X[0] ~ "," ~ delta ~ "," ~ beta[-1]))
  ) +
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
    scale_x_continuous(breaks = seq(0, max(cbf_prior$diff_beta), 0.02), guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~ facet, labeller = label_parsed, ncol = 1)

print(standard_dev_plot)


ggsave(filename = "sd_GLM_Poisson.pdf",path = "Plots", plot = standard_dev_plot,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)


##  ............................................................................
##  Posterior mean of beta 1                                                ####

mean_data <- data.frame()

for (model in 1:length(cbf_models)) {
  
  cbf_mean <- mean(cbf_models[[model]]$x1)
  unif_mean <- mean(unif_models[[model]]$x1)
  rmap_mean <- mean(rmap_models[[model]]$post.samples$x1)
  samp_mean <- mean(samp_models[[model]]$post.samples$x1)
  
  mean_vec <- c(cbf_mean, unif_mean, rmap_mean, samp_mean)
  
  mean_data <- rbind.data.frame(mean_data, mean_vec)
  
}

beta1_curr <- c(0.50, 0.53, 0.54, 0.55, 0.57, 0.58, 0.59,
                0.60, 0.61, 0.62, 0.63, 0.64, 0.65)

diff_beta <- beta1_curr - 0.5

mean_data <- cbind(mean_data, beta1_curr, diff_beta)

colnames(mean_data) <- c("NPP_cbf", "NPP_unif", "RMAP", "SAM", "beta1_curr", "diff_beta")


df_long_mean <- pivot_longer(mean_data, cols = c(NPP_unif, NPP_cbf,
                                                 RMAP, SAM),
                             names_to = "mean_type", values_to = "value")

df_long_mean <- df_long_mean %>%
  mutate(
    method = case_when(
      mean_type == "NPP_unif" ~ "Unif.",
      mean_type == "NPP_cbf"  ~ "CBF",
      mean_type == "RMAP"     ~ "RMAP",
      mean_type == "SAM"     ~ "SAM"
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
  "Unif."   = expression(NPP["Unif."]),
  "RMAP"    = "RMAP",
  "SAM"     = "SAM"
)

# Reorder to control legend order
df_compare <- df_compare %>%
  mutate(method = factor(method, levels = c("CBF", "Unif.", "RMAP", "SAM")))

mean_plot <- ggplot(df_compare, aes(x = diff_beta, y = value, color = method, shape = method, group = method)) +
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
    labs(
    x = element_blank(),
    y = element_blank(),
    title  = expression("Mean of " * pi(beta[1] ~ "|" ~ y ~ "," ~ y[0] ~ "," ~ X ~ ","~ X[0] ~ "," ~ delta ~ "," ~ beta[-1]))
  ) +
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
    scale_x_continuous(breaks = seq(0, max(cbf_prior$diff_beta), 0.02), guide = guide_axis(check.overlap = TRUE)) +
    facet_wrap(~ facet, labeller = label_parsed, ncol = 1)  

print(mean_plot)

ggsave(filename = "mean_GLM_Poisson.pdf",path = "Plots", plot = mean_plot,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)

##  ............................................................................
##  Combined plot posterior for beta 1                                      ####

plot_comb_beta <- ggarrange(standard_dev_plot, mean_plot, ncol = 2, nrow = 1, 
                                   common.legend = TRUE)
plot_comb_beta <- annotate_figure(plot_comb_beta, left = textGrob("Posterior Values", rot = 90, vjust = 0.8, gp = gpar(cex = 1.8)),
                                         bottom = textGrob(expression("Values of " * (beta["c"][","][1] - beta[0][","][1])), vjust = -0.1,  gp = gpar(cex = 1.8)))


ggsave(filename = "combined_plot_GLM_Poisson.pdf",path = "Plots", plot = plot_comb_beta,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)

##  ............................................................................
##  Plot posterior SD and mean of delta                                     ####

delta_data <- data.frame()

for (model in seq_along(cbf_models)) {
  cbf_delta <- cbf_models[[model]]$a0_hist_1
  unif_delta <- unif_models[[model]]$a0_hist_1
  
  diff_beta_index <- rep(beta1_curr[model] - 0.5, length(cbf_delta))
  
  delta_vec <- data.frame(
    cbf_delta = cbf_delta,
    unif_delta = unif_delta,
    diff_beta_index = diff_beta_index
  )
  
  delta_data <- rbind(delta_data, delta_vec)
}

delta_long <- delta_data %>%
  pivot_longer(
    cols = c("cbf_delta", "unif_delta"),
    names_to = "model_type",
    values_to = "value"
  )

# Compute mean and sd by diff_beta_index and model_type
delta_summary <- delta_long %>%
  group_by(diff_beta_index, model_type) %>%
  summarize(
    mean_value = mean(value),
    sd_value = sd(value),
    .groups = "drop"
  )

# Reshape summary data into long format for mean and sd
delta_summary_long <- delta_summary %>%
  pivot_longer(
    cols = c(mean_value, sd_value),
    names_to = "stat",
    values_to = "stat_value"
  )

# SD first and Mean second
delta_summary_long <- delta_summary_long %>%
  mutate(stat_label = case_when(
    stat == "mean_value" ~ "Mean~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~X~\",\"~X[0]~\",\"~bold(beta))",
    stat == "sd_value"   ~ "SD~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~X~\",\"~X[0]~\",\"~bold(beta))"
  ))

delta_summary_long$stat_label <- factor(delta_summary_long$stat_label, 
                                        levels = c(
                                          "SD~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~X~\",\"~X[0]~\",\"~bold(beta))",
                                          "Mean~of~pi(delta~\"|\"~y~\",\"~y[0]~\",\"~X~\",\"~X[0]~\",\"~bold(beta))"
                                        ))

delta_summary_long <- delta_summary_long %>%
  mutate(
    method = case_when(
      model_type == "unif_delta" ~ "Unif.",
      model_type == "cbf_delta"  ~ "CBF",
      TRUE ~ "Other" 
    )
  )

cb_palette <- c(
  "CBF"   = "#8A0910",  
  "RMAP"  = "#E69F00",  
  "Unif." = "#0072B2",  
  "SAMP"  = "#009E73"   
)

cb_shapes <- c(
  "CBF"   = 19,  
  "RMAP"  = 17,  
  "Unif." = 18,  
  "SAMP"  = 15   
)

label_expressions <- c(
  "CBF"     = expression(NPP["CBF"]),
  "Unif."   = expression(NPP["Unif."]),
  "RMAP"    = "RMAP",
  "SAMP"    = "SAMP"
)


GLM_poisson_delta_plot <- ggplot(delta_summary_long, aes(x = diff_beta_index, y = stat_value, color = method, shape = method)) +
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
    x = expression("Values of " * (beta["c"][","][1] - beta[0][","][1])),
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
  scale_x_continuous(breaks = seq(0, max(cbf_prior$diff_beta), 0.02), guide = guide_axis(check.overlap = TRUE)) +
  facet_wrap(~ stat_label, scales = "free_y", labeller = label_parsed)

print(GLM_poisson_delta_plot)


ggsave(filename = "GLM_poisson_delta_plot.pdf",path = "Plots", plot = GLM_poisson_delta_plot,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)
