#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "ggpubr",
  "dplyr",
  "ggplot2",
  "tidyverse",
  "tidyr",
  "grid",
  "tibble",
  "survival",
  "survminer",
  "hdbayes",
  "gridExtra"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#   ____________________________________________________________________________
#   Datasets                                                                ####

data("E2696")
data("E1694")
historical_data <- E2696
current_data <- E1694
historical_data$trial_name <- "E2696 Trial"
current_data$trial_name <- "E1694 Trial"

#   ____________________________________________________________________________
#   Survival curves                                                         ####

hd_surv <- survfit(Surv(
  historical_data$failtime,
  historical_data$failind
) ~ treatment, data = historical_data)
cd_surv <- survfit(Surv(
  current_data$failtime,
  current_data$failind
) ~ treatment, data = current_data)

# Historical dataset
time_cutoff <- 23
survs_hd <- summary(hd_surv, times = time_cutoff)$surv
labels_hd <- paste(round(survs_hd * 100), "%", sep = "")

plot_hd <- ggsurvplot(hd_surv,
  data = historical_data, combine = TRUE,
  facet.by = "trial_name",
  short.panel.labs = TRUE,
  pval = FALSE,
  conf.int = TRUE,
  onf.int.style = "step",
  palette = c("#4F7942", "#F46D43"),
  xlab = NULL,
  ylab = NULL,
  ggtheme = theme_bw(),
  risk.table = "percentage", # Add risk table
  risk.table.col = "strata", # Color the risk table by groups
  risk.table.height = 0.2, # Adjust risk table height
  linetype = 1, # Different line types for groups
  size = 1, # Adjust line size
  legend.labs = c("GMK", "GMK + IFN"), # Custom legend labels
  legend.title = "Treatment:",
  title = "",
  conf.int.style = "step",
  font.main = c(22),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(16),
  tables.theme = theme_cleantable(),
  risk.table.fontsize = 4.5,
  break.time.by = 5
) +
  theme_bw(base_size = 18) +
  scale_y_continuous(sec.axis = sec_axis(~., name = NULL, breaks = survs_hd, labels = labels_hd)) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    strip.text.x = element_text(size = 18)
  )

# Current dataset
time_cutoff <- 58
survs_cd <- summary(cd_surv, times = time_cutoff)$surv
labels_cd <- paste(round(survs_cd * 100), "%", sep = "")

plot_cd <- ggsurvplot_facet(cd_surv,
  data = current_data, combine = TRUE,
  facet.by = "trial_name",
  short.panel.labs = TRUE,
  pval = FALSE,
  conf.int = TRUE,
  onf.int.style = "step",
  palette = c("#4F7942", "#F46D43"),
  xlab = NULL,
  ylab = NULL,
  ggtheme = theme_bw(),
  risk.table = "percentage",
  risk.table.col = "strata",
  risk.table.height = 0.2,
  linetype = 1,
  size = 1,
  legend.labs = c("GMK", "IFN"),
  legend.title = "Treatment:",
  title = "",
  conf.int.style = "step",
  font.main = c(22),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(16),
  tables.theme = theme_cleantable(),
  risk.table.fontsize = 4.5,
  break.time.by = 5
) +
  theme_bw(base_size = 18) +
  scale_y_continuous(sec.axis = sec_axis(
    ~.,
    name = NULL,
    breaks = survs_cd,
    labels = labels_cd
  )) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    strip.text.x = element_text(size = 18)
  )

plot_cd

# Grid plot
melanoma_surv <- grid.arrange(plot_hd, plot_cd, nrow = 2)
melanoma_surv <- annotate_figure(melanoma_surv,
  left = textGrob("Survival Probability", rot = 90, vjust = 0.7, gp = gpar(cex = 1.8)),
  bottom = textGrob("Months", vjust = -0.1, gp = gpar(cex = 1.8))
)

ggsave(
  filename = "melanoma_surv.pdf", path = "Plots", plot = melanoma_surv,
  width = 15, height = 10, device = "pdf", dpi = 500, useDingbats = FALSE
)
