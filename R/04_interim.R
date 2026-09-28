# ============================================================
# 04_interim.R  (extended r grid for Table 5 AND Figure 3)
# Interim control-arm re-estimation of Vm with limited follow-up.
# Initial n sourced from 03_misspec.R runs with the misspecification
# factors reported in Table 5 (s = 0.5, 0.7).
# r_frac = information fraction r at the interim look.
# Requires: results/03_misspec_raw.rds  (extended-grid version)
# ============================================================

source("R/00_functions.R")
if (!dir.exists("results")) stopifnot(dir.create("results"))
set.seed(2025)

res_misspec <- readRDS("results/03_misspec_raw.rds")

lambda0_T <- 0.01
tau0_T_   <- 400
delta_    <- 1.0
Vm_true_  <- 2.0

r_vals       <- seq(0.1, 0.9, by = 0.1)   # fine grid for Figure 3
s_vm_interim <- c(0.5, 0.7)               # Table 5 blocks
power_target <- 0.8
# Grid: 2 x 9 x 5000 = 90k interim runs, each generating a full trial
# cohort -- heavier per run than Study 3. Parallel: library(furrr);
# plan(multisession); use furrr::future_pmap_dfr(...,
# .options = furrr_options(seed = TRUE)) in place of pmap_dfr below.

stage1_source <- res_misspec %>%
  dplyr::filter(tau0_T == tau0_T_, Vm_true == Vm_true_,
                s_vm %in% s_vm_interim) %>%
  dplyr::transmute(s_vm, iter, n_init = n_planned)

grid2 <- tidyr::crossing(stage1_source, r_frac = r_vals) %>%
  dplyr::mutate(tau0_I = tau0_T_ * r_frac)

res_interim <- purrr::pmap_dfr(
  grid2,
  function(s_vm, iter, n_init, r_frac, tau0_I) {
    run_pipeline_once_interim_control(
      n_init    = n_init,
      lambda0_I = lambda0_T,
      tau0_I    = tau0_I,
      lambda0_T = lambda0_T,
      tau0_T    = tau0_T_,
      delta     = delta_,
      Vm_true   = Vm_true_,
      beta_alt  = log(0.8),
      M0        = 1,
      alpha     = 0.05,
      power_target = power_target,
      p1 = 0.5
    ) %>%
      dplyr::mutate(s_vm = s_vm, iter = iter,
                    r_frac = r_frac, tau0_I = tau0_I)
  }
)

saveRDS(res_interim, "results/04_interim_raw.rds")

summary_interim <- res_interim %>%
  dplyr::group_by(s_vm, r_frac, tau0_I, Vm_true) %>%
  dplyr::summarise(
    n_init_mean = mean(n_init, na.rm = TRUE),
    n_init_sd   = sd(n_init, na.rm = TRUE),
    
    Vm_hat_interim_mean = mean(Vm_hat_interim, na.rm = TRUE),
    se_Vm_interim_mean  = mean(se_Vm_interim, na.rm = TRUE),
    
    n_oracle = median(n_oracle, na.rm = TRUE),
    
    n_update_mean = round(mean(n_update, na.rm = TRUE), 0),
    n_update_sd   = round(sd(n_update, na.rm = TRUE), 1),
    n_final_mean  = round(mean(n_final, na.rm = TRUE), 0),
    n_final_sd    = round(sd(n_final, na.rm = TRUE), 1),
    
    power_mean_final   = mean(power_final, na.rm = TRUE),
    power_sd_final     = sd(power_final, na.rm = TRUE),
    power_prob80_final = mean(power_final >= 0.80, na.rm = TRUE),
    power_prob75_final = mean(power_final >= 0.75, na.rm = TRUE),
    
    n_final_plusSE_mean = round(mean(n_final_plusSE, na.rm = TRUE), 0),
    n_final_plusSE_sd   = round(sd(n_final_plusSE, na.rm = TRUE), 1),
    power_mean_plusSE   = mean(power_final_plusSE, na.rm = TRUE),
    power_sd_plusSE     = sd(power_final_plusSE, na.rm = TRUE),
    power_prob80_plusSE = mean(power_final_plusSE >= 0.80, na.rm = TRUE),
    power_prob75_plusSE = mean(power_final_plusSE >= 0.75, na.rm = TRUE),
    
    n_final_plusHalfSE_mean = round(mean(n_final_plusHalfSE, na.rm = TRUE), 0),
    n_final_plusHalfSE_sd   = round(sd(n_final_plusHalfSE, na.rm = TRUE), 1),
    power_mean_plusHalfSE   = mean(power_final_plusHalfSE, na.rm = TRUE),
    power_sd_plusHalfSE     = sd(power_final_plusHalfSE, na.rm = TRUE),
    power_prob80_plusHalfSE = mean(power_final_plusHalfSE >= 0.80, na.rm = TRUE),
    power_prob75_plusHalfSE = mean(power_final_plusHalfSE >= 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

write_csv(summary_interim, "results/04_interim_summary.csv")

# Table 5 rows (subset of the fine r grid)
table5 <- summary_interim %>%
  dplyr::filter(round(r_frac, 1) %in% c(0.1, 0.3, 0.5, 0.7, 0.9))
write_csv(table5, "results/04_interim_table5.csv")
print(table5, n = Inf, width = Inf)