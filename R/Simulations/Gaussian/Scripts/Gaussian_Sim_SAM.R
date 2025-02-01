#   ____________________________________________________________________________
#   Libraries                                                               ####

packages <- c(
  "repmix",
  "SAMprior",
  "modi"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))

#   ____________________________________________________________________________
#   Starting Simulations                                                    ####

par_seq <- seq(0, 6, by = 0.2)

theta_seq <- seq(-20, 20, by = 0.01)

# Historical Data Parameters
true_theta_0 <- 0
true_se_0 <- 1

# Create empty lists to store the mean and sd for each iteration
mean_list <- numeric(length(par_seq))
sd_list <- numeric(length(par_seq))
w_SAM_list <- numeric(length(par_seq))


for (M in 1:length(par_seq)) {
  cat("Doing theta =", par_seq[M], "\n")

  ##  ............................................................................
  ##  Data                                                                    ####

  set.seed(4231)

  # Current Data Parameters
  true_theta <- par_seq[M]
  true_se <- 1

  ##  ............................................................................
  ##  Model and Grids Specification                                           ####

  prior_historical <- mixnorm(inf = c(1, 0, 1))

  w_SAM <- SAM_weight(
    if.prior = prior_historical,
    theta.h = 0,
    delta = 3, ## Clinically significant difference
    n = 1,
    m = true_theta,
    sigma = true_se
  )

  post_theta_mix <- thetaposteriormix(
    theta = theta_seq, tr = true_theta,
    to = true_theta_0, sr = true_se, so = true_se_0,
    w = w_SAM, m = 0, v = 100
  )

  sd_theta_sam_mix <- sqrt(weighted.var(theta_seq, post_theta_mix))
  mean_theta_sam_mix <- weighted.mean(theta_seq, post_theta_mix)

  # Store results in the lists
  mean_list[M] <- mean_theta_sam_mix
  sd_list[M] <- sd_theta_sam_mix
  w_SAM_list[M] <- w_SAM
}

SAM_dataset <- data.frame(
  mean_theta_sam_mix = mean_list,
  sd_theta_sam_mix = sd_list,
  w_SAM = w_SAM_list
)

print(SAM_dataset)

save(SAM_dataset, file = "R/Simulations/Gaussian/RData/SAM_prior/SAM_prior.RData")
