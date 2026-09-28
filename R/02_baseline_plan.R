# ============================================================
# 02_baseline_plan.R
# Baseline case: historical Vm estimation -> plug-in sample size
# -> power reliability (correctly specified historical Vm)
# ============================================================

source("R/00_functions.R")
dir.create("results", showWarnings = FALSE)
set.seed(2025)

lambda0_T   <- 0.01
tau0_T_vals <- c(100, 200, 400)
delta_      <- 0.10
Vm_vals <- c(0.2, 2.0)

n_hist_vals  <- 400
K            <- 5000
power_target <- 0.8

grid <- tidyr::expand_grid(tau0_T = tau0_T_vals, n_hist = n_hist_vals, Vm_ = Vm_vals)

res_baseline <- purrr::pmap_dfr(grid, function(tau0_T, n_hist, Vm_) {
  purrr::map_dfr(seq_len(K), ~ run_pipeline_once_sep2(
    n_hist    = n_hist,
    lambda0_H = 0.004,
    tau0_H    = 200,
    lambda0_T = lambda0_T,
    tau0_T    = tau0_T,
    delta     = delta_,
    Vm_hist   = Vm_,
    Vm_true   = Vm_,
    beta_alt  = log(0.8),
    M0        = 1,
    alpha     = 0.05,
    power_target = power_target,
    p1 = 0.5
  )) %>%
    dplyr::mutate(tau0_T = tau0_T, n_hist = n_hist)
})

saveRDS(res_baseline, "results/02_baseline_raw.rds")

summary_baseline <- summarise_plan(res_baseline, tau0_T, n_hist, Vm_true)
write_csv(summary_baseline, "results/02_baseline_summary.csv")
print(summary_baseline, n = Inf, width = Inf)
