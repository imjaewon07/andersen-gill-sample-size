# ============================================================
# 00_functions.R
# Shared functions for Vm estimation and sample size simulations
# Paper: Sample size calculation for Andersen-Gill models with
#        estimated between-subject heterogeneity
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ------------------------------------------------------------
# Data generation: counts from a mixed Poisson process
#   tau_i ~ Unif(tau0(1-delta), tau0(1+delta))
#   m_i   ~ frailty distribution with E[m_i] = 1, Var(m_i) = Vm
#   N_i | tau_i, m_i ~ Poisson(lambda0 m_i tau_i)
#
# The estimator targets Vm = Var(m_i) without distributional
# assumptions on the frailty; `frailty` selects the mixing law
# used in the DGP:
#   "gamma"     Gamma(1/Vm, 1/Vm)      -> NB marginal (sampled directly)
#   "lognormal" LN(-s2/2, s2), s2 = log(1+Vm)
#   "invgauss"  IG(mean = 1, shape = 1/Vm)  (Var = mean^3/shape)
# ------------------------------------------------------------

# Inverse-Gaussian sampler (Michael, Schucany & Haas, 1976); no extra packages
rinvgauss_ <- function(n, mean = 1, shape = 1){
  nu <- rnorm(n)
  y  <- nu^2
  x  <- mean + mean^2*y/(2*shape) -
    mean/(2*shape)*sqrt(4*mean*shape*y + mean^2*y^2)
  u  <- runif(n)
  ifelse(u <= mean/(mean + x), x, mean^2/x)
}

r_counts_fast <- function(n, lambda0, Vm, tau0, delta,
                          frailty = c("gamma", "lognormal", "invgauss")){
  stopifnot(delta >= 0 && delta <= 1)
  frailty <- match.arg(frailty)
  tau <- runif(n, min = tau0*(1-delta), max = tau0*(1+delta))
  
  if (Vm <= 0) {
    N <- rpois(n, lambda = lambda0 * tau)
  } else if (frailty == "gamma") {
    # NB marginal: no need to draw the frailty explicitly
    eta <- 1/Vm
    p <- eta / (eta + lambda0 * tau)
    N <- rnbinom(n, size = eta, prob = p)
  } else {
    m <- switch(frailty,
                lognormal = { s2 <- log(1 + Vm); rlnorm(n, -s2/2, sqrt(s2)) },
                invgauss  = rinvgauss_(n, mean = 1, shape = 1/Vm)
    )
    N <- rpois(n, lambda = lambda0 * m * tau)
  }
  tibble(i = seq_len(n), N = N, tau = tau)
}

# ------------------------------------------------------------
# Orthogonalized method-of-moments estimator of Vm
#   psi = (N - mu)^2 - N - Vm (2 mu N - mu^2), mu = lambda_hat * tau
# ------------------------------------------------------------
estimate_vm_orth_mom <- function(N, tau, truncate_nonneg = FALSE){
  stopifnot(length(N) == length(tau))
  
  lambda_hat <- sum(N) / sum(tau)
  mu_hat <- lambda_hat * tau
  
  num <- sum((N - mu_hat)^2 - N)
  den <- sum(2 * mu_hat * N - mu_hat^2)
  
  Vm_hat <- num / den
  if (truncate_nonneg) Vm_hat <- pmax(Vm_hat, 0)
  
  list(
    lambda_hat = lambda_hat,
    Vm_hat = Vm_hat,
    mu_hat = mu_hat,
    numerator = num,
    denominator = den
  )
}

# Sandwich SE for Vm_hat.
# df_correct = TRUE applies a sqrt(n/(n-1)) small-sample correction;
# default FALSE to match the manuscript's reported SE.
var_vm_sandwich_orth <- function(N, tau, lambda_hat, Vm_hat, df_correct = FALSE){
  mu_hat <- lambda_hat * tau
  
  psi_hat <- (N - mu_hat)^2 - N -
    Vm_hat * (2 * mu_hat * N - mu_hat^2)
  
  den <- sum(2 * mu_hat * N - mu_hat^2)
  
  se_Vm <- sqrt(sum(psi_hat^2)) / abs(den)
  if (df_correct) {
    n <- length(N)
    se_Vm <- se_Vm * sqrt(n / (n - 1))
  }
  
  list(
    se_Vm = se_Vm,
    var_Vm = se_Vm^2,
    psi_hat = psi_hat,
    denominator = den
  )
}

# ------------------------------------------------------------
# Sample size planning (superiority / NI) and approximate power
# ------------------------------------------------------------
vbeta_plugin <- function(beta, p1 = 0.5, p0 = 0.5,
                         lambda0, tau0, delta, Vm_hat){
  term1 <- (1/(p1*exp(beta)) + 1/p0) * 1/(lambda0*tau0)
  term2 <- (1/p1 + 1/p0) * Vm_hat * (1 + delta^2/3)
  term1 + term2
}

n_plan_supNI <- function(beta_alt, M0, alpha = 0.05, power = 0.90,
                         p1 = 0.5, p0 = 0.5,
                         lambda0, tau0, delta, Vm_hat){
  z1 <- qnorm(1 - alpha/2)
  z2 <- qnorm(power)
  Vb <- vbeta_plugin(beta = beta_alt, p1 = p1, p0 = p0,
                     lambda0 = lambda0, tau0 = tau0, delta = delta,
                     Vm_hat = Vm_hat)
  n <- ((z1 + z2)^2 * Vb) / ((log(M0) - beta_alt)^2)
  ceiling(n)
}

power_approx <- function(n, beta_alt, M0, alpha,
                         p1 = 0.5, p0 = 0.5,
                         lambda0, tau0, delta, Vm_true){
  z1 <- qnorm(1 - alpha/2)
  Vb_true <- (1/(p1*exp(beta_alt)) + 1/p0) * 1/(lambda0*tau0) +
    (1/p1 + 1/p0) * Vm_true * (1 + delta^2/3)
  zpow <- sqrt(n) * abs(log(M0) - beta_alt) / sqrt(Vb_true) - z1
  pnorm(zpow)
}

# ------------------------------------------------------------
# Trial data generation
# ------------------------------------------------------------
# Two-arm counts, fixed allocation round(n_total * p1)
gen_trial_counts <- function(n_total, p1 = 0.5,
                             lambda0, Vm_true, tau0, delta, beta_true){
  n1 <- round(n_total * p1)
  n0 <- n_total - n1
  
  dat1 <- r_counts_fast(n = n1, lambda0 = lambda0*exp(beta_true),
                        Vm = Vm_true, tau0 = tau0, delta = delta)
  dat0 <- r_counts_fast(n = n0, lambda0 = lambda0,
                        Vm = Vm_true, tau0 = tau0, delta = delta)
  tibble(
    N   = c(dat1$N, dat0$N),
    tau = c(dat1$tau, dat0$tau),
    x   = c(rep(1, n1), rep(0, n0))
  )
}

# Linked interim/final counts for the same subjects
# (independent-increments construction); r_frac = information fraction
gen_trial_counts_linked <- function(
    n_total, p1 = 0.5,
    lambda0, Vm_true, tau0, delta, beta_true, r_frac
){
  n_total <- as.integer(ceiling(n_total))
  n1 <- round(n_total * p1)
  x_arm <- c(rep(1L, n1), rep(0L, n_total - n1))
  
  if (Vm_true > 0) {
    frailty_m <- rgamma(n_total, shape = 1 / Vm_true, scale = Vm_true)
  } else {
    frailty_m <- rep(1, n_total)
  }
  
  tau_full_obs    <- runif(n_total,
                           min = tau0 * (1 - delta),
                           max = tau0 * (1 + delta))
  tau_interim_obs <- r_frac * tau_full_obs
  
  rate_i <- lambda0 * exp(beta_true * x_arm) * frailty_m
  
  count_interim <- rpois(n_total, lambda = rate_i * tau_interim_obs)
  count_extra   <- rpois(n_total, lambda = rate_i * (tau_full_obs - tau_interim_obs))
  count_full    <- count_interim + count_extra
  
  tibble(
    id = seq_len(n_total),
    x_arm = x_arm,
    frailty_m = frailty_m,
    tau_interim_obs = tau_interim_obs,
    tau_full_obs = tau_full_obs,
    count_interim = count_interim,
    count_full = count_full
  )
}

# ------------------------------------------------------------
# Final Wald test: closed-form two-group Poisson MLE + HC0 sandwich.
# Algebraically identical to
#   glm(N ~ x + offset(log(tau)), family = poisson) + vcovHC(type = "HC0"):
# the two-group model is saturated, so lambda_g = sum(N_g)/sum(tau_g) and
#   Var(beta_hat) = sum_g sum_i (N_i - lambda_g tau_i)^2 / (sum_i N_gi)^2.
# ------------------------------------------------------------
final_wald_test <- function(dat, M0 = 1, alpha = 0.05){
  i1 <- dat$x == 1
  N1 <- dat$N[i1];  T1 <- dat$tau[i1]
  N0 <- dat$N[!i1]; T0 <- dat$tau[!i1]
  
  lam1 <- sum(N1) / sum(T1)
  lam0 <- sum(N0) / sum(T0)
  beta_hat <- log(lam1 / lam0)
  
  v1 <- sum((N1 - lam1 * T1)^2) / (sum(N1))^2
  v0 <- sum((N0 - lam0 * T0)^2) / (sum(N0))^2
  se_beta <- sqrt(v1 + v0)
  
  zcrit <- qnorm(1 - alpha/2)
  z_stat <- (beta_hat - log(M0)) / se_beta
  
  tibble(
    beta_hat = beta_hat,
    se_beta  = se_beta,
    z_stat   = z_stat,
    reject   = as.numeric(abs(z_stat) > zcrit)
  )
}

# ------------------------------------------------------------
# Pipeline: historical estimation -> planning -> approximate power
# ------------------------------------------------------------
run_pipeline_once_sep2 <- function(
    n_hist,
    lambda0_H, tau0_H,
    lambda0_T, tau0_T,
    delta,
    Vm_hist,
    Vm_true,
    beta_alt = log(0.8), M0 = 1, alpha = 0.05, power_target = 0.80,
    p1 = 0.5, truncate_nonneg = FALSE
){
  hist <- r_counts_fast(n = n_hist, lambda0 = lambda0_H,
                        Vm = Vm_hist, tau0 = tau0_H, delta = delta)
  
  est    <- estimate_vm_orth_mom(hist$N, hist$tau,
                                 truncate_nonneg = truncate_nonneg)
  Vm_hat <- est$Vm_hat
  seVm   <- var_vm_sandwich_orth(hist$N, hist$tau,
                                 est$lambda_hat, Vm_hat)$se_Vm
  
  plan_n <- function(V) n_plan_supNI(
    beta_alt = beta_alt, M0 = M0, alpha = alpha,
    power = power_target, p1 = p1, p0 = 1 - p1,
    lambda0 = lambda0_T, tau0 = tau0_T, delta = delta,
    Vm_hat = V
  )
  
  n_plan       <- plan_n(Vm_hat)
  n_oracle     <- plan_n(Vm_true)
  n_plusHalfSE <- plan_n(Vm_hat + 0.5 * seVm)
  n_plusSE     <- plan_n(Vm_hat + seVm)
  
  pwr <- function(n) power_approx(
    n = n, beta_alt = beta_alt, M0 = M0, alpha = alpha,
    p1 = p1, p0 = 1 - p1,
    lambda0 = lambda0_T, tau0 = tau0_T, delta = delta,
    Vm_true = Vm_true
  )
  
  tibble(
    lambda0_H = lambda0_H, tau0_H = tau0_H,
    lambda0_T = lambda0_T, tau0_T = tau0_T,
    delta = delta,
    Vm_hist = Vm_hist, Vm_true = Vm_true,
    Vm_hat = Vm_hat, se_Vm = seVm,
    n_planned = n_plan,
    n_oracle = n_oracle,
    n_plusHalfSE = n_plusHalfSE,
    n_plusSE = n_plusSE,
    power_approx_plan = pwr(n_plan),
    power_approx_plusHalfSE = pwr(n_plusHalfSE),
    power_approx_plusSE = pwr(n_plusSE)
  )
}

# ------------------------------------------------------------
# Pipeline: interim control-arm re-estimation -> approximate power
# (unlinked interim data; power evaluated analytically)
# ------------------------------------------------------------
run_pipeline_once_interim_control <- function(
    n_init,
    lambda0_I, tau0_I,
    lambda0_T, tau0_T,
    delta,
    Vm_true,
    beta_alt = log(0.8), M0 = 1, alpha = 0.05, power_target = 0.80,
    p1 = 0.5, truncate_nonneg = FALSE
){
  interim <- gen_trial_counts(
    n_total = n_init, p1 = p1,
    lambda0 = lambda0_I, Vm_true = Vm_true,
    tau0 = tau0_I, delta = delta, beta_true = beta_alt
  )
  
  ctrl <- interim %>% dplyr::filter(x == 0)
  
  est    <- estimate_vm_orth_mom(ctrl$N, ctrl$tau,
                                 truncate_nonneg = truncate_nonneg)
  Vm_hat <- est$Vm_hat
  seVm   <- var_vm_sandwich_orth(ctrl$N, ctrl$tau,
                                 est$lambda_hat, Vm_hat)$se_Vm
  
  plan_n <- function(V) n_plan_supNI(
    beta_alt = beta_alt, M0 = M0, alpha = alpha,
    power = power_target, p1 = p1, p0 = 1 - p1,
    lambda0 = lambda0_T, tau0 = tau0_T, delta = delta,
    Vm_hat = V
  )
  
  n_update            <- plan_n(Vm_hat)
  n_oracle            <- plan_n(Vm_true)
  n_update_plusHalfSE <- plan_n(Vm_hat + 0.5 * seVm)
  n_update_plusSE     <- plan_n(Vm_hat + seVm)
  
  n_final            <- max(n_init, n_update)
  n_final_plusHalfSE <- max(n_init, n_update_plusHalfSE)
  n_final_plusSE     <- max(n_init, n_update_plusSE)
  
  pwr <- function(n) power_approx(
    n = n, beta_alt = beta_alt, M0 = M0, alpha = alpha,
    p1 = p1, p0 = 1 - p1,
    lambda0 = lambda0_T, tau0 = tau0_T, delta = delta,
    Vm_true = Vm_true
  )
  
  tibble(
    n_init = n_init,
    lambda0_I = lambda0_I, tau0_I = tau0_I,
    lambda0_T = lambda0_T, tau0_T = tau0_T,
    delta = delta, Vm_true = Vm_true,
    Vm_hat_interim = Vm_hat,
    se_Vm_interim  = seVm,
    n_oracle = n_oracle,
    n_update = n_update,
    n_final  = n_final,
    n_update_plusHalfSE = n_update_plusHalfSE,
    n_update_plusSE     = n_update_plusSE,
    n_final_plusHalfSE  = n_final_plusHalfSE,
    n_final_plusSE      = n_final_plusSE,
    power_final            = pwr(n_final),
    power_final_plusHalfSE = pwr(n_final_plusHalfSE),
    power_final_plusSE     = pwr(n_final_plusSE)
  )
}

# ------------------------------------------------------------
# Pipeline: full two-stage design with linked interim/final data
# and empirical final Wald test (used for type I error)
# ------------------------------------------------------------
run_pipeline_with_final_test <- function(
    n_hist,
    lambda0_H, tau0_H,
    Vm_hist,
    lambda0_T, tau0_T,
    Vm_true,
    beta_true,
    r_frac,
    delta,
    beta_alt = log(0.8),
    M0 = 1,
    alpha = 0.05,
    power_target = 0.80,
    p1 = 0.5,
    truncate_nonneg = FALSE
){
  # Stage 1: initial planning from historical control
  hist <- r_counts_fast(n = n_hist, lambda0 = lambda0_H,
                        Vm = Vm_hist, tau0 = tau0_H, delta = delta)
  
  est_hist    <- estimate_vm_orth_mom(hist$N, hist$tau,
                                      truncate_nonneg = truncate_nonneg)
  Vm_hat_hist <- est_hist$Vm_hat
  
  plan_n <- function(V) as.integer(ceiling(n_plan_supNI(
    beta_alt = beta_alt, M0 = M0, alpha = alpha,
    power = power_target, p1 = p1, p0 = 1 - p1,
    lambda0 = lambda0_T, tau0 = tau0_T, delta = delta,
    Vm_hat = V
  )))
  
  n_init <- plan_n(Vm_hat_hist)
  
  # Stage 2: initial trial subjects (linked interim/final counts)
  trial_init <- gen_trial_counts_linked(
    n_total = n_init, p1 = p1,
    lambda0 = lambda0_T, Vm_true = Vm_true,
    tau0 = tau0_T, delta = delta,
    beta_true = beta_true, r_frac = r_frac
  )
  
  ctrl <- trial_init %>% dplyr::filter(x_arm == 0)
  
  est_int    <- estimate_vm_orth_mom(ctrl$count_interim, ctrl$tau_interim_obs,
                                     truncate_nonneg = truncate_nonneg)
  Vm_hat_int <- est_int$Vm_hat
  seVm_int   <- var_vm_sandwich_orth(ctrl$count_interim, ctrl$tau_interim_obs,
                                     est_int$lambda_hat, Vm_hat_int)$se_Vm
  
  # Stage 3: updated total sample sizes
  n_update            <- plan_n(Vm_hat_int)
  n_update_plusHalfSE <- plan_n(Vm_hat_int + 0.5 * seVm_int)
  n_update_plusSE     <- plan_n(Vm_hat_int + seVm_int)
  
  n_final            <- max(n_init, n_update)
  n_final_plusHalfSE <- max(n_init, n_update_plusHalfSE)
  n_final_plusSE     <- max(n_init, n_update_plusSE)
  
  # Extra subjects generated once for the largest final n, so the three
  # analysis sets are nested and comparable
  n_final_max <- max(n_final, n_final_plusHalfSE, n_final_plusSE)
  n_extra <- n_final_max - n_init
  
  if (n_extra > 0) {
    trial_extra <- gen_trial_counts_linked(
      n_total = n_extra, p1 = p1,
      lambda0 = lambda0_T, Vm_true = Vm_true,
      tau0 = tau0_T, delta = delta,
      beta_true = beta_true, r_frac = r_frac
    )
    trial_extra$id <- n_init + seq_len(nrow(trial_extra))
    trial_all <- dplyr::bind_rows(trial_init, trial_extra)
  } else {
    trial_all <- trial_init
  }
  
  make_final_dat <- function(trial_data, n_keep) {
    idx <- trial_data$id <= as.integer(ceiling(n_keep))
    data.frame(
      N   = trial_data$count_full[idx],
      tau = trial_data$tau_full_obs[idx],
      x   = trial_data$x_arm[idx]
    )
  }
  
  # Stage 4/5: final full-follow-up Wald tests
  test_main <- final_wald_test(make_final_dat(trial_all, n_final),
                               M0 = M0, alpha = alpha)
  test_half <- final_wald_test(make_final_dat(trial_all, n_final_plusHalfSE),
                               M0 = M0, alpha = alpha)
  test_se   <- final_wald_test(make_final_dat(trial_all, n_final_plusSE),
                               M0 = M0, alpha = alpha)
  
  tibble(
    n_hist = n_hist,
    r_frac = r_frac,
    lambda0_H = lambda0_H, tau0_H = tau0_H,
    lambda0_T = lambda0_T, tau0_T = tau0_T,
    delta = delta,
    Vm_hist = Vm_hist, Vm_true = Vm_true,
    beta_true = beta_true,
    Vm_hat_hist = Vm_hat_hist,
    Vm_hat_interim = Vm_hat_int,
    se_Vm_interim = seVm_int,
    n_init = n_init,
    n_update = n_update,
    n_final = n_final,
    n_update_plusHalfSE = n_update_plusHalfSE,
    n_update_plusSE = n_update_plusSE,
    n_final_plusHalfSE = n_final_plusHalfSE,
    n_final_plusSE = n_final_plusSE,
    beta_hat_final = test_main$beta_hat,
    se_beta_final = test_main$se_beta,
    reject_final = test_main$reject,
    reject_final_plusHalfSE = test_half$reject,
    reject_final_plusSE = test_se$reject
  )
}

# ------------------------------------------------------------
# Shared summary helper for planning studies
# ------------------------------------------------------------
summarise_plan <- function(df, ...){
  df %>%
    group_by(...) %>%
    summarise(
      Vm_hat_mean = mean(Vm_hat, na.rm = TRUE),
      n_oracle = median(n_oracle, na.rm = TRUE),
      
      n_mean = round(mean(n_planned, na.rm = TRUE), 0),
      n_sd   = round(sd(n_planned, na.rm = TRUE), 1),
      power_mean_plan   = mean(power_approx_plan, na.rm = TRUE),
      power_sd_plan     = sd(power_approx_plan, na.rm = TRUE),
      power_prob80_plan = mean(power_approx_plan >= 0.80, na.rm = TRUE),
      power_prob75_plan = mean(power_approx_plan >= 0.75, na.rm = TRUE),
      
      n_plusSE_mean = round(mean(n_plusSE, na.rm = TRUE), 0),
      n_plusSE_sd   = round(sd(n_plusSE, na.rm = TRUE), 1),
      power_mean_plusSE   = mean(power_approx_plusSE, na.rm = TRUE),
      power_sd_plusSE     = sd(power_approx_plusSE, na.rm = TRUE),
      power_prob80_plusSE = mean(power_approx_plusSE >= 0.80, na.rm = TRUE),
      power_prob75_plusSE = mean(power_approx_plusSE >= 0.75, na.rm = TRUE),
      
      n_plusHalfSE_mean = round(mean(n_plusHalfSE, na.rm = TRUE), 0),
      n_plusHalfSE_sd   = round(sd(n_plusHalfSE, na.rm = TRUE), 1),
      power_mean_plusHalfSE   = mean(power_approx_plusHalfSE, na.rm = TRUE),
      power_sd_plusHalfSE     = sd(power_approx_plusHalfSE, na.rm = TRUE),
      power_prob80_plusHalfSE = mean(power_approx_plusHalfSE >= 0.80, na.rm = TRUE),
      power_prob75_plusHalfSE = mean(power_approx_plusHalfSE >= 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}