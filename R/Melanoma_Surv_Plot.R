
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
  "grid",
  "xtable",
  "tibble",
  "survival",
  "survminer",
  "hdbayes",
  "gridExtra"
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
#   Datasets                                                                ####


data("E2696")
data("E1694")
historical_data <- E2696
current_data <- E1694
historical_data$trial_name <- 'E2696 Trial'
current_data$trial_name <- 'E1694 Trial'


hd_surv <- survfit( Surv(historical_data$failtime, historical_data$failind) ~ treatment, data = historical_data)
cd_surv <- survfit( Surv(current_data$failtime, current_data$failind) ~ treatment, data = current_data)

time_cutoff <- 23
survs_hd = summary(hd_surv, times=time_cutoff)$surv
labels_hd = paste(round(survs_hd*100), '%', sep='')



plot_hd <- ggsurvplot(hd_surv, data = historical_data, combine = TRUE, 
           facet.by = "trial_name",
           short.panel.labs = TRUE,
           pval = FALSE, 
           conf.int = TRUE, 
           onf.int.style = "step",
           palette = c("#4F7942", "#F46D43"),
           xlab = NULL, 
           ylab = NULL,
           ggtheme = theme_bw(),         # Cleaner theme
           risk.table = "percentage",                 # Add risk table
           risk.table.col = "strata",         # Color the risk table by groups
           risk.table.height = 0.2,           # Adjust risk table height
           linetype = 1,               # Different line types for groups
           size = 1,                          # Adjust line size
           legend.labs = c("GMK", "GMK + IFN"), # Custom legend labels
           legend.title = "Treatment:",
           title = "",                        # Plot title
           conf.int.style = "step",
           font.main = c(22),
           font.x = c(20),
           font.y = c(20),
           font.tickslab = c(16),
           tables.theme = theme_cleantable(),
           risk.table.fontsize = 4.5,
           break.time.by = 5) +
   theme_bw(base_size = 18) +
  scale_y_continuous(sec.axis=sec_axis(~., name=NULL, breaks = survs_hd, labels= labels_hd))+
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

time_cutoff <- 58
survs_cd = summary(cd_surv, times=time_cutoff)$surv
labels_cd = paste(round(survs_cd*100), '%', sep='')



plot_cd <-
  ggsurvplot_facet(
    cd_surv,
    data = current_data,
    combine = TRUE,
    facet.by = "trial_name",
    short.panel.labs = TRUE,
    pval = FALSE,
    conf.int = TRUE,
    onf.int.style = "step",
    palette = c("#4F7942", "#F46D43"),
    # Improved color palette
    xlab = NULL,
    ylab = NULL,
    ggtheme = theme_bw(),
    # Cleaner theme
    risk.table = "percentage",
    # Add risk table
    risk.table.col = "strata",
    # Color the risk table by groups
    risk.table.height = 0.2,
    # Adjust risk table height
    linetype = 1,
    # Different line types for groups
    size = 1,
    # Adjust line size
    legend.labs = c("GMK", "IFN"),
    # Custom legend labels
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
  ) + theme_bw(base_size = 18) +
  scale_y_continuous(sec.axis=sec_axis(~., name=NULL, 
                                                  breaks = survs_cd,
                                                  labels= labels_cd))+
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

melanoma_surv <- grid.arrange(plot_hd, plot_cd, nrow = 2)

melanoma_surv <- annotate_figure(melanoma_surv, left = textGrob("Survival Probability", rot = 90, vjust = 0.7, gp = gpar(cex = 1.8)),
                                         bottom = textGrob("Months", vjust = -0.1,  gp = gpar(cex = 1.8)))



ggsave(filename = "melanoma_surv.pdf",path = "Plots", plot = melanoma_surv,
       width = 15, height = 10, device='pdf', dpi=500, useDingbats = FALSE)

# ggsurvplot_combine(
#   list(hd_surv, cd_surv),
#   list(historical_data, current_data),
#   pval = FALSE, 
#   conf.int = TRUE, 
#   onf.int.style = "step",
#   xlab = "Months", 
#   ylab = "Survival Probability",
#   ggtheme = theme_bw(),         # Cleaner theme
#   risk.table = "percentage",                 # Add risk table
#   risk.table.col = "strata",         # Color the risk table by groups
#   risk.table.height = 0.2,           # Adjust risk table height
#   linetype = 1,               # Different line types for groups
#   size = 1,                          # Adjust line size
#   # legend.labs = c("GMK", "GMK + IFN","GMK", "IFN"), # Custom legend labels
#   # legend.title = c("Treatment:"),
#   title = "",                        # Plot title
#   # surv.median.line = "hv",            # Add median survival lines
#   font.main = c(22),
#   font.x = c(20),
#   font.y = c(20),
#   font.tickslab = c(16),
#   tables.theme = theme_cleantable(),
#   risk.table.fontsize = 4.5,
#   break.time.by = 5)
# 
# 
# 
# 
# 
# 
# # Assuming historical_data and current_data are already loaded
# historical_data$dataset_type <- 'Historical'  # Add a new column to distinguish the data
# current_data$dataset_type <- 'Current'
# 
# # Combine the data frames
# combined_data <- rbind(historical_data, current_data)
# 
# # Convert treatment_group to a factor for faceting
# combined_data$dataset_type <- factor(combined_data$dataset_type, levels = c("Historical", "Historical"))
# 
# combined_surv <- survfit(Surv(failtime, failind) ~ treatment , data = combined_data)
# ggsurv <- ggsurvplot(combined_surv, data = combined_data, combine = TRUE, # Combine curves
#                      pval = FALSE, 
#                      conf.int = TRUE, 
#                      onf.int.style = "step",
#                      palette = c("deepskyblue3", "firebrick3"), # Improved color palette
#                      xlab = "Months", 
#                      ylab = "Survival Probability",
#                      ggtheme = theme_bw(),         # Cleaner theme
#                      risk.table = "percentage",                 # Add risk table
#                      risk.table.col = "strata",         # Color the risk table by groups
#                      risk.table.height = 0.2,           # Adjust risk table height
#                      linetype = 1,               # Different line types for groups
#                      size = 1,                          # Adjust line size
#                      legend.labs = c("GMK", "IFN"), # Custom legend labels
#                      legend.title = "Treatment:",
#                      title = "",                        # Plot title
#                      # surv.median.line = "hv",            # Add median survival lines
#                      font.main = c(22),
#                      font.x = c(20),
#                      font.y = c(20),
#                      font.tickslab = c(16),
#                      tables.theme = theme_cleantable(),
#                      risk.table.fontsize = 4.5,
#                      break.time.by = 5)
# 
# 
# # Load survminer for plotting
# library(survminer)
# ggsurvplot_facet(combined_surv, data = combined_data, facet.by = "dataset_type",
#                  pval = FALSE, conf.int = TRUE,
#                  xlab = "Months", ylab = "Survival Probability",
#                  ggtheme = theme_bw(), legend.labs = c("GMK", "IFN"), # Custom legend labels
#                 nrow = 2,
#                  palette = "simpsons", scales = "free_x")
