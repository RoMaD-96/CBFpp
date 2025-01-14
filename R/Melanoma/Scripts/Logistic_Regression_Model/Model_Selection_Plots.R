#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "RBesT",
  "SAMprior",
  "ggplot2",
  "ggpubr",
  "ggridges",
  "hdbayes",
  "dplyr",
  "ggplot2",
  "rethinking",
  "tidyr",
  "tidyverse"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

source("R/Melanoma/Scripts/Logistic_Regression_Model/Fun_Post_Pred.R", echo = TRUE)
source("R/CBF_Functions.R")

load("R/Melanoma/RData/bf_list.RData")
load("R/Melanoma/RData/obs_rep.RData")
load("R/Melanoma/RData/obs_bf_melanoma.RData")

#   ____________________________________________________________________________
#   Grid Search                                                             ####

eta <- round(seq(0.5, 6, by = 0.5), 2)
nu <- round(seq(0.5, 6, by = 0.5), 2)
grid <- expand.grid(x = eta, y = nu)
grid <- as.data.frame(grid)
model_id <- seq(1, nrow(grid))
grid <- cbind(grid, model_id)
model_id <- model_id[-14]
grid <- grid[-14, ]

#   ____________________________________________________________________________
#   Optimal Model                                                           ####

bf_data_rep <- as.data.frame(do.call(rbind, lapply(bf_list, unlist)))
colnames(bf_data_rep) <- paste0("BF_rep", 1:100)
print(bf_data_rep)

bf_data_long <- bf_data_rep %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = -row_id, names_to = "Column", values_to = "Value")

df_obs_bf <- filter(df_obs_bf, obs_bf >= 0)
optimal_model_index <- calculate_maximized_value(bf_data_long$Value, df_obs_bf$obs_bf, 0.75)

grid[optimal_model_index,]

#   ____________________________________________________________________________
#   Datasets                                                                ####

set.seed(442)

data("E2696")
data("E1694")
historical_data <- E2696
current_data <- E1694
all_data <- list(current_data, historical_data)

formula <- failind ~ age + treatment + sex + perform
p <- length(attr(terms(formula), "term.labels"))
family <- binomial('logit')

#   ____________________________________________________________________________
#   Approximate normalizing constant                                        ####

a0 <- seq(0, 1, length.out = 21)

if (requireNamespace("parallel")) {
  if (instantiate::stan_cmdstan_exists()) {
    logncfun <- function(a0, ...) {
      hdbayes::glm.npp.lognc(
        formula = formula, family = family, histdata = historical_data, a0 = a0, ...
      )
    }
    cl <- parallel::makeCluster(10)
    parallel::clusterSetRNGStream(cl, 123)
    parallel::clusterExport(cl, varlist = c('formula', 'family', 'historical_data'))
    time.npp.1 <- system.time(
      a0.lognc <- parLapply(
        cl = cl, X = a0, fun = logncfun, iter_warmup = 1000,
        iter_sampling = 1000, chains = 4, refresh = 0
      )
    )
    parallel::stopCluster(cl)
  }
  a0.lognc <- data.frame(do.call(rbind, a0.lognc))
}

#   ____________________________________________________________________________
#   Model Specification                                                     ####


##  ............................................................................
##  Normalized Power Prior (NPP) mixture of beta                            ####

model_npp_mix <- glm.npp_mix(
  formula = formula, family = family, data.list = all_data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = 0.5, 
  a0.shape2 = 6,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Normalized Power Prior (NPP) CBF                                        ####

model_opt_bf <- glm.npp(
  formula = formula, family = family, data.list = all_data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = 5, 
  a0.shape2 = 0.5,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Normalized Power Prior (NPP) MSE Shen et. al. (2023)                    ####

model_mse <- glm.npp(
  formula = formula, family = family, data.list = all_data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = 6, 
  a0.shape2 = 0.5,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Normalized Power Prior (NPP) Uniform                                    ####

model_unif <- glm.npp(
  formula = formula, family = family, data.list = all_data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = 1, a0.shape2 = 1, 
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Normalized Power Prior (NPP) Jeffreys                                   ####

model_jeffreys <- glm.npp(
  formula = formula, family = family, data.list = all_data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = 0.5, a0.shape2 = 0.5, 
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Normalized Power Prior (NPP) Beta(2,2)                                  ####

model_beta_2 <- glm.npp(
  formula = formula, family = family, data.list = all_data,
  a0.lognc = a0.lognc$a0,
  lognc = matrix(a0.lognc$lognc, ncol = 1),
  a0.shape1 = 2, a0.shape2 = 2,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Robust Meta Analytic Prior (RMAP)                                       ####

model_RMAP <- glm.rmap(
  formula = formula, family = family, data.list = all_data,
  w = 0.5,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Commensurate Prior (CP)                                                 ####

model_CP <- glm.commensurate(
  formula = formula,
  family = family, data.list = all_data,
  p.spike = 0.3,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

##  ............................................................................
##  Self-Adapting Mixture Prior (SAMP)                                      ####

N_0 <- nrow(historical_data)
y_0 <- nrow(historical_data[historical_data$failind == 1,])
p_0 <- y_0 / N_0

N <- nrow(current_data)
y <- nrow(current_data[current_data$failind == 1,])

alpha_h <- y_0 + 1
beta_h <- N_0 + 1 - y_0
alpha_v <- 1
beta_v <- 1

prior_historical <- mixbeta(inf = c(1, alpha_h, beta_h))
wSAM <- SAM_weight(if.prior = prior_historical, delta = 0.1, n = N, r = y)
print(wSAM)

model_SAMP <- glm.rmap(
  formula = formula, family = family, data.list = all_data,
  w = wSAM,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  refresh = 100, seed = 433
)

#   ____________________________________________________________________________
#   Table                                                                   ####

data_parameter_post <- data.frame(
  # Age Parameter Estimates
  age_beta_opt = model_opt_bf$age,
  age_beta_unif = model_unif$age,
  age_beta_jeffreys = model_jeffreys$age,
  age_beta_2 = model_beta_2$age,
  age_beta_RMAP = model_RMAP$post.samples$age,
  age_beta_SAMP = model_SAMP$post.samples$age,
  age_beta_CP = model_CP$age,
  age_beta_npp_mix = model_npp_mix$age,  
  
  # Treatment Parameter Estimates
  treat_beta_opt = model_opt_bf$treatment,
  treat_beta_unif = model_unif$treatment,
  treat_beta_jeffreys = model_jeffreys$treatment,
  treat_beta_2 = model_beta_2$treatment,
  treat_beta_RMAP = model_RMAP$post.samples$treatment,
  treat_beta_SAMP = model_SAMP$post.samples$treatment,
  treat_beta_CP = model_CP$treatment,
  treat_beta_npp_mix = model_npp_mix$treatment,  
  
  # Sex Parameter Estimates
  sex_beta_opt = model_opt_bf$sex,
  sex_beta_unif = model_unif$sex,
  sex_beta_jeffreys = model_jeffreys$sex,
  sex_beta_2 = model_beta_2$sex,
  sex_beta_RMAP = model_RMAP$post.samples$sex,
  sex_beta_SAMP = model_SAMP$post.samples$sex,
  sex_beta_CP = model_CP$sex,
  sex_beta_npp_mix = model_npp_mix$sex,  
  
  # Performance Parameter Estimates
  perform_beta_opt = model_opt_bf$perform,
  perform_beta_unif = model_unif$perform,
  perform_beta_jeffreys = model_jeffreys$perform,
  perform_beta_2 = model_beta_2$perform,
  perform_beta_RMAP = model_RMAP$post.samples$perform,
  perform_beta_SAMP = model_SAMP$post.samples$perform,
  perform_beta_CP = model_CP$perform,
  perform_beta_npp_mix = model_npp_mix$perform,  
  
  # Delta Parameter Estimates
  delta_opt = model_opt_bf$a0_hist_1,
  delta_unif = model_unif$a0_hist_1,
  delta_opt_jeffreys = model_jeffreys$a0_hist_1, 
  delta_beta_2 = model_beta_2$a0_hist_1,
  delta_npp_mix = model_npp_mix$a0_hist_1  
)


# HPDI 
# Calculate HPDIs for Age Parameters
hpdi_age_opt <- HPDI(data_parameter_post$age_beta_opt, prob = 0.95)
hpdi_age_unif <- HPDI(data_parameter_post$age_beta_unif, prob = 0.95)
hpdi_age_jeffreys <- HPDI(data_parameter_post$age_beta_jeffreys, prob = 0.95)
hpdi_age_beta_2 <- HPDI(data_parameter_post$age_beta_2, prob = 0.95)
hpdi_age_RMAP <- HPDI(data_parameter_post$age_beta_RMAP, prob = 0.95)
hpdi_age_SAMP <- HPDI(data_parameter_post$age_beta_SAMP, prob = 0.95)
hpdi_age_CP <- HPDI(data_parameter_post$age_beta_CP, prob = 0.95)
hpdi_age_npp_mix <- HPDI(data_parameter_post$age_beta_npp_mix, prob = 0.95)  

# Calculate HPDIs for Treatment Parameters
hpdi_treat_opt <- HPDI(data_parameter_post$treat_beta_opt, prob = 0.95)
hpdi_treat_unif <- HPDI(data_parameter_post$treat_beta_unif, prob = 0.95)
hpdi_treat_jeffreys <- HPDI(data_parameter_post$treat_beta_jeffreys, prob = 0.95)
hpdi_treat_beta_2 <- HPDI(data_parameter_post$treat_beta_2, prob = 0.95)
hpdi_treat_RMAP <- HPDI(data_parameter_post$treat_beta_RMAP, prob = 0.95)
hpdi_treat_SAMP <- HPDI(data_parameter_post$treat_beta_SAMP, prob = 0.95)
hpdi_treat_CP <- HPDI(data_parameter_post$treat_beta_CP, prob = 0.95)
hpdi_treat_npp_mix <- HPDI(data_parameter_post$treat_beta_npp_mix, prob = 0.95)  

# Calculate HPDIs for Sex Parameters
hpdi_sex_opt <- HPDI(data_parameter_post$sex_beta_opt, prob = 0.95)
hpdi_sex_unif <- HPDI(data_parameter_post$sex_beta_unif, prob = 0.95)
hpdi_sex_jeffreys <- HPDI(data_parameter_post$sex_beta_jeffreys, prob = 0.95)
hpdi_sex_beta_2 <- HPDI(data_parameter_post$sex_beta_2, prob = 0.95)
hpdi_sex_RMAP <- HPDI(data_parameter_post$sex_beta_RMAP, prob = 0.95)
hpdi_sex_SAMP <- HPDI(data_parameter_post$sex_beta_SAMP, prob = 0.95)
hpdi_sex_CP <- HPDI(data_parameter_post$sex_beta_CP, prob = 0.95)
hpdi_sex_npp_mix <- HPDI(data_parameter_post$sex_beta_npp_mix, prob = 0.95)  

# Calculate HPDIs for Performance Parameters
hpdi_perform_opt <- HPDI(data_parameter_post$perform_beta_opt, prob = 0.95)
hpdi_perform_unif <- HPDI(data_parameter_post$perform_beta_unif, prob = 0.95)
hpdi_perform_jeffreys <- HPDI(data_parameter_post$perform_beta_jeffreys, prob = 0.95)
hpdi_perform_beta_2 <- HPDI(data_parameter_post$perform_beta_2, prob = 0.95)
hpdi_perform_RMAP <- HPDI(data_parameter_post$perform_beta_RMAP, prob = 0.95)
hpdi_perform_SAMP <- HPDI(data_parameter_post$perform_beta_SAMP, prob = 0.95)
hpdi_perform_CP <- HPDI(data_parameter_post$perform_beta_CP, prob = 0.95)
hpdi_perform_npp_mix <- HPDI(data_parameter_post$perform_beta_npp_mix, prob = 0.95)  


# Dataframe for Latex Table
data_parameter_post_tab <- data.frame(
  prior = c("Beta(5,0.5)", 
            "Beta(1,1)", 
            "Beta(0.5,0.5)", 
            "Beta(2,2)", 
            "RMAP", 
            "SAM", 
            "CP", 
            "NPP_Mix"),  # Added NPP_Mix
  
  # Mean Estimates for Age
  mean_age = format(round(c(
    mean(data_parameter_post$age_beta_opt),
    mean(data_parameter_post$age_beta_unif),
    mean(data_parameter_post$age_beta_jeffreys),
    mean(data_parameter_post$age_beta_2),
    mean(data_parameter_post$age_beta_RMAP),
    mean(data_parameter_post$age_beta_SAMP),
    mean(data_parameter_post$age_beta_CP),
    mean(data_parameter_post$age_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # Mean Estimates for Treatment
  mean_treat = format(round(c(
    mean(data_parameter_post$treat_beta_opt),
    mean(data_parameter_post$treat_beta_unif),
    mean(data_parameter_post$treat_beta_jeffreys),
    mean(data_parameter_post$treat_beta_2),
    mean(data_parameter_post$treat_beta_RMAP),
    mean(data_parameter_post$treat_beta_SAMP),
    mean(data_parameter_post$treat_beta_CP),
    mean(data_parameter_post$treat_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # Mean Estimates for Sex
  mean_sex = format(round(c(
    mean(data_parameter_post$sex_beta_opt),
    mean(data_parameter_post$sex_beta_unif),
    mean(data_parameter_post$sex_beta_jeffreys),
    mean(data_parameter_post$sex_beta_2),
    mean(data_parameter_post$sex_beta_RMAP),
    mean(data_parameter_post$sex_beta_SAMP),
    mean(data_parameter_post$sex_beta_CP),
    mean(data_parameter_post$sex_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # Mean Estimates for Performance
  mean_perform = format(round(c(
    mean(data_parameter_post$perform_beta_opt),
    mean(data_parameter_post$perform_beta_unif),
    mean(data_parameter_post$perform_beta_jeffreys),
    mean(data_parameter_post$perform_beta_2),
    mean(data_parameter_post$perform_beta_RMAP),
    mean(data_parameter_post$perform_beta_SAMP),
    mean(data_parameter_post$perform_beta_CP),
    mean(data_parameter_post$perform_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # Standard Deviation for Age
  sd_age = format(round(c(
    sd(data_parameter_post$age_beta_opt),
    sd(data_parameter_post$age_beta_unif),
    sd(data_parameter_post$age_beta_jeffreys),
    sd(data_parameter_post$age_beta_2),
    sd(data_parameter_post$age_beta_RMAP),
    sd(data_parameter_post$age_beta_SAMP),
    sd(data_parameter_post$age_beta_CP),
    sd(data_parameter_post$age_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # Standard Deviation for Treatment
  sd_treat = format(round(c(
    sd(data_parameter_post$treat_beta_opt),
    sd(data_parameter_post$treat_beta_unif),
    sd(data_parameter_post$treat_beta_jeffreys),
    sd(data_parameter_post$treat_beta_2),
    sd(data_parameter_post$treat_beta_RMAP),
    sd(data_parameter_post$treat_beta_SAMP),
    sd(data_parameter_post$treat_beta_CP),
    sd(data_parameter_post$treat_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # Standard Deviation for Sex
  sd_sex = format(round(c(
    sd(data_parameter_post$sex_beta_opt),
    sd(data_parameter_post$sex_beta_unif),
    sd(data_parameter_post$sex_beta_jeffreys),
    sd(data_parameter_post$sex_beta_2),
    sd(data_parameter_post$sex_beta_RMAP),
    sd(data_parameter_post$sex_beta_SAMP),
    sd(data_parameter_post$sex_beta_CP),
    sd(data_parameter_post$sex_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # Standard Deviation for Performance
  sd_perform = format(round(c(
    sd(data_parameter_post$perform_beta_opt),
    sd(data_parameter_post$perform_beta_unif),
    sd(data_parameter_post$perform_beta_jeffreys),
    sd(data_parameter_post$perform_beta_2),
    sd(data_parameter_post$perform_beta_RMAP),
    sd(data_parameter_post$perform_beta_SAMP),
    sd(data_parameter_post$perform_beta_CP),
    sd(data_parameter_post$perform_beta_npp_mix)  
  ), 3), nsmall = 3),
  
  # HPDI for Age Parameters
  hpdi_age = c(
    paste("(", round(hpdi_age_opt[1], 3), ",", round(hpdi_age_opt[2], 3), ")"),
    paste("(", round(hpdi_age_unif[1], 3), ",", round(hpdi_age_unif[2], 3), ")"),
    paste("(", round(hpdi_age_jeffreys[1], 3), ",", round(hpdi_age_jeffreys[2], 3), ")"),
    paste("(", round(hpdi_age_beta_2[1], 3), ",", round(hpdi_age_beta_2[2], 3), ")"),
    paste("(", round(hpdi_age_RMAP[1], 3), ",", round(hpdi_age_RMAP[2], 3), ")"),
    paste("(", round(hpdi_age_SAMP[1], 3), ",", round(hpdi_age_SAMP[2], 3), ")"),
    paste("(", round(hpdi_age_CP[1], 3), ",", round(hpdi_age_CP[2], 3), ")"),
    paste("(", round(hpdi_age_npp_mix[1], 3), ",", round(hpdi_age_npp_mix[2], 3), ")")  
  ),
  
  # HPDI for Treatment Parameters
  hpdi_treat = c(
    paste("(", round(hpdi_treat_opt[1], 3), ",", round(hpdi_treat_opt[2], 3), ")"),
    paste("(", round(hpdi_treat_unif[1], 3), ",", round(hpdi_treat_unif[2], 3), ")"),
    paste("(", round(hpdi_treat_jeffreys[1], 3), ",", round(hpdi_treat_jeffreys[2], 3), ")"),
    paste("(", round(hpdi_treat_beta_2[1], 3), ",", round(hpdi_treat_beta_2[2], 3), ")"),
    paste("(", round(hpdi_treat_RMAP[1], 3), ",", round(hpdi_treat_RMAP[2], 3), ")"),
    paste("(", round(hpdi_treat_SAMP[1], 3), ",", round(hpdi_treat_SAMP[2], 3), ")"),
    paste("(", round(hpdi_treat_CP[1], 3), ",", round(hpdi_treat_CP[2], 3), ")"),
    paste("(", round(hpdi_treat_npp_mix[1], 3), ",", round(hpdi_treat_npp_mix[2], 3), ")")  
  ),
  
  # HPDI for Sex Parameters
  hpdi_sex = c(
    paste("(", round(hpdi_sex_opt[1], 3), ",", round(hpdi_sex_opt[2], 3), ")"),
    paste("(", round(hpdi_sex_unif[1], 3), ",", round(hpdi_sex_unif[2], 3), ")"),
    paste("(", round(hpdi_sex_jeffreys[1], 3), ",", round(hpdi_sex_jeffreys[2], 3), ")"),
    paste("(", round(hpdi_sex_beta_2[1], 3), ",", round(hpdi_sex_beta_2[2], 3), ")"),
    paste("(", round(hpdi_sex_RMAP[1], 3), ",", round(hpdi_sex_RMAP[2], 3), ")"),
    paste("(", round(hpdi_sex_SAMP[1], 3), ",", round(hpdi_sex_SAMP[2], 3), ")"),
    paste("(", round(hpdi_sex_CP[1], 3), ",", round(hpdi_sex_CP[2], 3), ")"),
    paste("(", round(hpdi_sex_npp_mix[1], 3), ",", round(hpdi_sex_npp_mix[2], 3), ")")  
  ),
  
  # HPDI for Performance Parameters
  hpdi_perform = c(
    paste("(", round(hpdi_perform_opt[1], 3), ",", round(hpdi_perform_opt[2], 3), ")"),
    paste("(", round(hpdi_perform_unif[1], 3), ",", round(hpdi_perform_unif[2], 3), ")"),
    paste("(", round(hpdi_perform_jeffreys[1], 3), ",", round(hpdi_perform_jeffreys[2], 3), ")"),
    paste("(", round(hpdi_perform_beta_2[1], 3), ",", round(hpdi_perform_beta_2[2], 3), ")"),
    paste("(", round(hpdi_perform_RMAP[1], 3), ",", round(hpdi_perform_RMAP[2], 3), ")"),
    paste("(", round(hpdi_perform_SAMP[1], 3), ",", round(hpdi_perform_SAMP[2], 3), ")"),
    paste("(", round(hpdi_perform_CP[1], 3), ",", round(hpdi_perform_CP[2], 3), ")"),
    paste("(", round(hpdi_perform_npp_mix[1], 3), ",", round(hpdi_perform_npp_mix[2], 3), ")")  
  )
)

#   ____________________________________________________________________________
#   Plots                                                                   ####


##  ............................................................................
##  Posterior densities                                                     ####

long_data_age <- gather(
  data_parameter_post,
  key = "variable",
  value = "value",
  age_beta_opt, age_beta_unif, age_beta_RMAP, age_beta_SAMP, age_beta_CP, age_beta_npp_mix
)

long_data_sex <- gather(
  data_parameter_post,
  key = "variable",
  value = "value",
  sex_beta_opt, sex_beta_unif, sex_beta_RMAP, sex_beta_SAMP, sex_beta_CP, sex_beta_npp_mix
)

long_data_perform <- gather(
  data_parameter_post,
  key = "variable",
  value = "value",
  perform_beta_opt, perform_beta_unif, perform_beta_RMAP, perform_beta_SAMP, perform_beta_CP, perform_beta_npp_mix
)

long_data_treat <- gather(
  data_parameter_post,
  key = "variable",
  value = "value",
  treat_beta_opt, treat_beta_unif, treat_beta_RMAP, treat_beta_SAMP, treat_beta_CP, treat_beta_npp_mix
)

data_combined <- bind_rows(
  mutate(long_data_age, category = "Age"),
  mutate(long_data_sex, category = "Sex"),
  mutate(long_data_perform, category = "Performance"),
  mutate(long_data_treat, category = "Treatment")
)

data_combined$prior_type <- case_when(
  grepl("opt", data_combined$variable) ~ "CBF",
  grepl("unif", data_combined$variable) ~ "Uniform",
  grepl("RMAP", data_combined$variable) ~ "RMAP",
  grepl("SAMP", data_combined$variable) ~ "SAM",
  grepl("CP", data_combined$variable) ~ "CP",
  grepl("npp_mix", data_combined$variable) ~ "NPP_Mix"
)

error_bars_df <- data.frame(
  category = rep(c("Age", "Sex", "Performance", "Treatment"), each = 6),
  prior_type = rep(c("CBF", "Uniform", "RMAP", "SAM", "CP", "NPP_Mix"), times = 4),
  xmin = c(
    hpdi_age_opt[1], hpdi_age_unif[1], hpdi_age_RMAP[1], hpdi_age_SAMP[1], hpdi_age_CP[1], hpdi_age_npp_mix[1],
    hpdi_sex_opt[1], hpdi_sex_unif[1], hpdi_sex_RMAP[1], hpdi_sex_SAMP[1], hpdi_sex_CP[1], hpdi_sex_npp_mix[1],
    hpdi_perform_opt[1], hpdi_perform_unif[1], hpdi_perform_RMAP[1], hpdi_perform_SAMP[1], hpdi_perform_CP[1], hpdi_perform_npp_mix[1],
    hpdi_treat_opt[1], hpdi_treat_unif[1], hpdi_treat_RMAP[1], hpdi_treat_SAMP[1], hpdi_treat_CP[1], hpdi_treat_npp_mix[1]
  ),
  xmax = c(
    hpdi_age_opt[2], hpdi_age_unif[2], hpdi_age_RMAP[2], hpdi_age_SAMP[2], hpdi_age_CP[2], hpdi_age_npp_mix[2],
    hpdi_sex_opt[2], hpdi_sex_unif[2], hpdi_sex_RMAP[2], hpdi_sex_SAMP[2], hpdi_sex_CP[2], hpdi_sex_npp_mix[2],
    hpdi_perform_opt[2], hpdi_perform_unif[2], hpdi_perform_RMAP[2], hpdi_perform_SAMP[2], hpdi_perform_CP[2], hpdi_perform_npp_mix[2],
    hpdi_treat_opt[2], hpdi_treat_unif[2], hpdi_treat_RMAP[2], hpdi_treat_SAMP[2], hpdi_treat_CP[2], hpdi_treat_npp_mix[2]
  )
)

facet_names <- c(
  Age = "Age", 
  Sex = "Sex", 
  Performance = "Performance Status", 
  Treatment = "Treatment"
)

fill_colors <- c(
  "CBF"      = "#8A0404",
  "NPP_Mix"  = "#AA4499",
  "Uniform"  = "#0072B2",  
  "RMAP"     = "#E69F00",    
  "SAM"      = "#009E20",   
  "CP"       = "#484F59"   
)

error_bar_colors <- c(
  "CBF"      = "#4B0000",
  "NPP_Mix"  = "#CC79A7",
  "Uniform"  = "#004C7A",  
  "RMAP"     = "#A65C00",  
  "SAM"      = "#006014",  
  "CP"       = "#1C232C"   
)

desired_order <- c("CBF", "NPP_Mix", "Uniform", "RMAP", "SAM", "CP")

data_combined$prior_type <- factor(data_combined$prior_type, levels = desired_order)
error_bars_df$prior_type <- factor(error_bars_df$prior_type, levels = desired_order)

label_expressions <- c(
  "CBF"     = expression(NPP["CBF"]),
  "Uniform" = expression(NPP["Unif."]),
  "NPP_Mix" = expression(NPP["Mix."]),
  "RMAP"    = "RMAP",
  "SAM"     = "SAM",
  "CP"      = "CP"
)

plot_combined <- ggplot(data_combined, aes(x = value, y = prior_type, fill = prior_type)) +
  geom_density_ridges(alpha = 0.4, scale = 0.8, rel_min_height = 0.01) +
  scale_fill_manual(values = fill_colors, labels = label_expressions) +
  facet_wrap(~category, scales = "free", labeller = as_labeller(facet_names)) +
  geom_errorbarh(
    data = error_bars_df, 
    aes(xmin = xmin, xmax = xmax, y = prior_type, color = prior_type),
    inherit.aes = FALSE,
    alpha = 1.0,
    height = 0.3,
    size = 1.2,
    show.legend = FALSE
  ) +
  geom_point(
    data = error_bars_df,
    aes(x = (xmin + xmax) / 2, y = prior_type, color = prior_type),
    inherit.aes = FALSE,
    size = 2,
    shape = 21,
    fill = "white",
    stroke = 1.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = error_bar_colors) +
  scale_y_discrete(labels = label_expressions, expand = c(0.3, 0)) +
  labs(y = "Prior Distribution", x = "Values") +
  guides(fill = guide_legend(title = expression(bold(paste("Prior: ")))) ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 18),
    strip.text.x = element_text(size = 18),
    strip.placement = "outside",
    strip.switch.pad.grid = unit(0.9, "in")
  )

print(plot_combined)

ggsave(
  filename = "post_param_melanoma.pdf",
  path = "Plots",
  plot = plot_combined,
  width = 15,
  height = 10,
  device = 'pdf',
  dpi = 500,
  useDingbats = FALSE
)


##  ............................................................................
##  SD differences in %                                                     ####

priors <- c("opt", "unif", "jeffreys", "RMAP", "SAM", "CP", "npp_mix")
categories <- c("age", "sex", "perform", "treat")
sd_diff_all <- data.frame()

for (cat in categories) {
  sd_opt <- sd(data_parameter_post[[paste0(cat, "_beta_opt")]])
  sd_unif <- sd(data_parameter_post[[paste0(cat, "_beta_unif")]])    
  sd_jeffreys <- sd(data_parameter_post[[paste0(cat, "_beta_jeffreys")]])
  sd_RMAP <- sd(data_parameter_post[[paste0(cat, "_beta_RMAP")]])
  sd_SAMP <- sd(data_parameter_post[[paste0(cat, "_beta_SAMP")]])
  sd_CP <- sd(data_parameter_post[[paste0(cat, "_beta_CP")]])
  sd_mix <- sd(data_parameter_post[[paste0(cat, "_beta_npp_mix")]])
  
  sd_ref <- sd_unif
  
  pct_diff_opt <- (sd_opt / sd_ref - 1) * 100
  pct_diff_unif <- (sd_unif / sd_ref - 1) * 100  
  pct_diff_jeffreys <- (sd_jeffreys / sd_ref - 1) * 100
  pct_diff_RMAP <- (sd_RMAP / sd_ref - 1) * 100
  pct_diff_SAMP <- (sd_SAMP / sd_ref - 1) * 100
  pct_diff_CP <- (sd_CP / sd_ref - 1) * 100
  pct_diff_mix <- (sd_mix / sd_ref - 1) * 100
  
  temp_df <- data.frame(
    category = toupper(cat),
    SDPercentDiff = c(
      pct_diff_opt,
      pct_diff_unif,
      pct_diff_jeffreys,
      pct_diff_RMAP,
      pct_diff_SAMP,
      pct_diff_CP,
      pct_diff_mix
    ),
    prior_type = c("CBF", "Uniform", "Jeffreys", "RMAP", "SAM", "CP", "NPP_Mix")
  )
  
  sd_diff_all <- bind_rows(sd_diff_all, temp_df)
}

prior_label <- c(
  "CBF"      = expression(NPP["CBF"]),
  "Jeffreys" = expression(NPP["Jeffr."]),
  "RMAP"     = "RMAP",
  "SAM"      = "SAM",
  "CP"       = "CP",
  "NPP_Mix"  = expression(NPP["Mix"])
)

facet_names <- c(
  "AGE" = "Age", 
  "SEX" = "Sex", 
  "PERFORM" = "Performance Status", 
  "TREAT" = "Treatment"
)

fill_colors <- c(
  "CBF"      = "#8A0404",  
  "Uniform"  = "#0072B2",  
  "Jeffreys" = "#56B4E9",  
  "RMAP"     = "#E69F00",  
  "SAM"      = "#009E20",  
  "CP"       = "#484F59",  
  "NPP_Mix"  = "#CC79A7"   
)

sd_diff_all <- sd_diff_all %>%
  group_by(category) %>%
  arrange(desc(SDPercentDiff)) %>%
  ungroup()

desired_order <- c("CBF", "NPP_Mix", "Uniform", "Jeffreys", "RMAP", "SAM", "CP")
sd_diff_all$prior_type <- factor(sd_diff_all$prior_type, levels = desired_order)

sd_diff_filtered <- sd_diff_all %>%
  filter(prior_type != "Uniform") %>%
  mutate(
    color_var = case_when(
      prior_type == "CP" ~ "CP",
      SDPercentDiff > 0 ~ "GREY",
      TRUE ~ as.character(prior_type)
    )
  )

fill_colors_with_grey <- c(
  "CBF"      = "#8A0404", 
  "Jeffreys" = "#56B4E9",
  "RMAP"     = "#E69F00",
  "SAM"      = "#009E20",
  "CP"       = "#484F59", 
  "NPP_Mix"  = "#CC79A7",
  "GREY"     = "#484F59"
)

prior_label <- c(
  "CBF"      = expression(NPP["CBF"]),
  "Jeffreys" = expression(NPP["Jeffr."]),
  "RMAP"     = "RMAP",
  "SAM"      = "SAM",
  "CP"       = "CP",
  "NPP_Mix"  = expression(NPP["Mix"])
)

lollipop <- ggplot(data = sd_diff_filtered, aes(x = prior_type, y = SDPercentDiff, color = color_var)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 1) +
  geom_segment(aes(x = prior_type, xend = prior_type, y = 0, yend = SDPercentDiff), size = 1.5) +
  geom_point(size = 4) +
  facet_wrap(~ category, scales = "free_y", labeller = as_labeller(facet_names)) +
  scale_color_manual(
    values = fill_colors_with_grey,
    labels = prior_label,
    breaks = c("CBF", "Jeffreys", "RMAP", "SAM", "CP", "NPP_Mix")
  ) +
  scale_x_discrete(labels = prior_label) +
  coord_flip() +
  labs(
    x = "Prior Type",
    y = "SD Difference (%)",
    color = "Prior: "
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    strip.text = element_text(size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "lollipop_melanoma.pdf",
  path = "Plots",
  plot = lollipop,
  width = 15,
  height = 10,
  device = "pdf",
  dpi = 500,
  useDingbats = FALSE
)