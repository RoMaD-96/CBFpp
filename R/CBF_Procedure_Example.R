library(ggplot2)
library(dplyr)
library(tidyverse)

# Create a data frame for the x values
x_values <- seq(-7, 7, length.out = 5000)

# Create a data frame with all distributions
df <- data.frame(x = x_values) %>%
  mutate(Dist1 = dnorm(x, mean = -1, sd = sqrt(1)),
         Dist2 = dnorm(x, mean = 0.5, sd = sqrt(1)),
         Dist3 = dnorm(x, mean = 2, sd = sqrt(1))) %>%
  pivot_longer(cols = starts_with("Dist"), names_to = "Distribution", values_to = "Density")

# Define the means and standard deviations
distributions <- data.frame(mean = c(-0.5, 0.5, 2), sd = sqrt(1), Distribution = c(paste0("LogBF[0][1]"),
                                                                                 paste0("LogBF[0][2]"), 
                                                                                 paste0("LogBF[0][3]")))

# Calculate the 75% HPDI for each distribution
distributions$lower <- qnorm(0.125, distributions$mean, distributions$sd)
distributions$upper <- qnorm(0.875, distributions$mean, distributions$sd)

# Points to be marked with a cross
cross_points <- data.frame(x = c(0.4, -0.2, 2.5), Distribution = c(paste0("LogBF[0][1]"),
                                                                  paste0("LogBF[0][2]"), 
                                                                  paste0("LogBF[0][3]")))

# Create a data frame for the x values
x_values <- seq(-5, 5, length.out = 1000)

# Create a data frame with all distributions
df <- expand.grid(x = x_values, Distribution = distributions$Distribution) %>%
  mutate(Density = case_when(
    Distribution == "LogBF[0][1]" ~ dnorm(x, mean = -0.5, sd = sqrt(1)),
    Distribution == "LogBF[0][2]" ~ dnorm(x, mean = 0.5, sd = sqrt(1)),
    Distribution == "LogBF[0][3]" ~ dnorm(x, mean = 2, sd = sqrt(1))
  ))


# Define a color palette
colors <-c("#7570B3", "#D95F02" , "#1B9E77") 
# Plot
cbf_example <- ggplot(df, aes(x = x, y = Density)) +
  geom_line(aes(color = Distribution)) +
  geom_area(data = df %>% filter(x > 0), aes(x = x, y = Density, fill = Distribution), alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 1.05) +
  geom_vline(data = distributions, aes(xintercept = lower, color = Distribution), linetype = "dashed", size = 1.05) +
  geom_vline(data = distributions, aes(xintercept = upper, color = Distribution), linetype = "dashed", size = 1.05) +
  geom_point(data = cross_points, aes(x = x, y = 0), shape = 4, size = 4.5, stroke = 1, color = "black") +
  facet_wrap(~ Distribution, scales = "free_y", nrow = 3, labeller = label_parsed) +
  scale_fill_manual(values = setNames(colors, distributions$Distribution)) +
  scale_color_manual(values = setNames(colors, distributions$Distribution)) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 18, angle = 0, vjust = 0),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18)
  ) +
  labs(title = "",
       x = "LogBF Values",
       y = "Density") +
  guides(fill = FALSE, color = FALSE) # Hide the legend

plot(cbf_example)

ggsave(filename = "CBF_example.pdf",path = "Plots", plot = cbf_example,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)
