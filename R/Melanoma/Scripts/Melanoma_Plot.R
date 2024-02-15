#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "doFuture",
  "parallel",
  "doParallel",
  "foreach",
  "rstan",
  "bridgesampling",
  "future",
  "dplyr",
  "progressr",
  "patchwork",
  "ggplot2",
  "ggridges",
  "tidyr",
  "survival",
  "survminer",
  "readxl"
)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


#   ____________________________________________________________________________
#   Datasets                                                                ####

data <- drop_na(read_excel("1.Main_Work_With_Ioannis/Uniform_Reference/Samples_From_Alternative/Melanoma/Data/trial_e1684_e1690_Merged.xlsx", 
                           col_types = c("numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric", "numeric", "numeric", 
                                         "numeric")))

historical_data <- filter(data, study == 1684)[,-2]
current_data <- filter(data, study == 1690)[,-2]

## Standardization ##

log_age_hist <- log(historical_data$age)
log_age_current <- log(current_data$age)




combined_data <- bind_rows(
  historical_data %>% mutate(Group = 'E1684'),
  current_data %>% mutate(Group = 'E1690')
)

# Create survival objects
combined_data$surv_obj <- with(combined_data, Surv(survtime, scens))

# Fit survival models
surv_fit <- survfit(surv_obj ~ Group, data = combined_data)

# Plot
plot_surv <- ggsurvplot(
  surv_fit, 
  data = combined_data,
  pval = FALSE, 
  conf.int = TRUE, 
  palette = c("deepskyblue3", "firebrick3"), # Improved color palette
  xlab = "Years", 
  ylab = "Survival Probability",
  ggtheme = theme_minimal(),         # Cleaner theme
  risk.table = "percentage",                 # Add risk table
  risk.table.col = "strata",         # Color the risk table by groups
  risk.table.height = 0.2,           # Adjust risk table height
  linetype = 1,               # Different line types for groups
  size = 1,                          # Adjust line size
  legend.labs = c("E1684", "E1690"), # Custom legend labels
  legend.title = "Trial:",
  title = "",                        # Plot title
  surv.median.line = "hv",            # Add median survival lines
  font.main = c(22),
  font.x = c(20),
  font.y = c(20),
  font.tickslab = c(16),
  tables.theme = theme_cleantable(),
  risk.table.fontsize = 4.5,
  break.time.by = 1)

plot_surv$plot <- plot_surv$plot + 
  theme(legend.text = element_text(size = 13),
        legend.title = element_text(size = 13))
plot_surv$table <- plot_surv$table + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
plot_surv


grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave(filename = "melanoma_surv_plot.pdf",path = "1.Main_Work_With_Ioannis/Uniform_Reference/Plots", plot = plot_surv,
       width = 12, height = 8, device='pdf', dpi=500, useDingbats = FALSE)
