# ============================================================
# 05_type1_error.R  (parallel version)
# Empirical type I error of the final robust Wald test when the
# same trial's control arm is used for interim Vm re-estimation.
# Linked interim/final counts per subject. Full grid simulated;
# the manuscript reports a subset (Table 6 / Figure 4).
#
# Parallelized with furrr. First: install.packages("furrr")
# ============================================================

source("R/00_functions.R")
if (!dir.exists("results")) stopifnot(dir.create("results"))

library(furrr)
plan(multisession, workers = max(1, parallelly::availableCores() - 1))

lambda0_H_ <- 0.004
tau0_H_    <- 200
n_hist_    <- 400

lambda0_T_ <- 0.01
tau0_T_    <- 400
delta_     <- 1.0
Vm_true_   <- 2.0

scale_vm <- c(0.5, 0.7, 0.8, 0.9)
r_vals   <- c(0.1, 0.3, 0.5, 0.7, 0.9)

K <- 20000

beta_true_null <- 0
# beta_true_null <- log(M0)  # use instead when evaluating the NI boundary

grid_type1 <- tidyr::crossing(
  s_vm = scale_vm,
  r_frac = r_vals,
  iter = seq_len(K)
)

# furrr manages per-element RNG streams via furrr_options(seed = TRUE),
# so no set.seed() here; results are reproducible for a fixed worker
# setup but will not match a sequential purrr run draw-for-draw.
res_type1 <- furrr::future_pmap_dfr(
  grid_type1,
  function(s_vm, r_frac, iter) {
    run_pipeline_with_final_test(
      n_hist    = n_hist_,
      lambda0_H = lambda0_H_,
      tau0_H    = tau0_H_,
      Vm_hist   = Vm_true_ * s_vm,
      
      lambda0_T = lambda0_T_,
      tau0_T    = tau0_T_,
      Vm_true   = Vm_true_,
      beta_true = beta_true_null,
      r_frac    = r_frac,
      
      delta     = delta_,
      beta_alt  = log(0.8),
      M0        = 1,
      alpha     = 0.05,
      power_target = 0.80,
      p1 = 0.5
    ) %>%
      dplyr::mutate(s_vm = s_vm, iter = iter)
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

plan(sequential)   # release the workers

saveRDS(res_type1, "results/05_type1_raw.rds")

summary_type1 <- res_type1 %>%
  dplyr::group_by(s_vm, r_frac) %>%
  dplyr::summarise(
    n_init_mean  = mean(n_init, na.rm = TRUE),
    n_final_mean = mean(n_final, na.rm = TRUE),
    
    Vm_hat_hist_mean    = mean(Vm_hat_hist, na.rm = TRUE),
    Vm_hat_interim_mean = mean(Vm_hat_interim, na.rm = TRUE),
    
    type1_main       = mean(reject_final, na.rm = TRUE),
    type1_plusHalfSE = mean(reject_final_plusHalfSE, na.rm = TRUE),
    type1_plusSE     = mean(reject_final_plusSE, na.rm = TRUE),
    
    mc_se_main       = sqrt(type1_main * (1 - type1_main) / dplyr::n()),
    mc_se_plusHalfSE = sqrt(type1_plusHalfSE * (1 - type1_plusHalfSE) / dplyr::n()),
    mc_se_plusSE     = sqrt(type1_plusSE * (1 - type1_plusSE) / dplyr::n()),
    
    .groups = "drop"
  )

write_csv(summary_type1, "results/05_type1_summary.csv")
print(summary_type1, n = Inf, width = Inf)