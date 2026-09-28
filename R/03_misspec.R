# ============================================================
# 03_misspec.R  (extended grid for Table 4 AND Figure 2)
# Sensitivity of planned sample size and power to historical
# misspecification of Vm: Vm_hist = s * Vm_true.
#
# Grid:
#   tau0_T  in {100, 200, 400}
#   Vm_true in {0.2, 2.0}
#   s       on a fine grid for the figure; the table rows
#           (s = 0.5, 0.7, 0.8, 0.9, 1.1) are a subset.
#
# Raw results feed 04_interim.R (which filters tau0_T == 400,
# s in {0.5, 0.7}) and the figure scripts.
# ============================================================

source("R/00_functions.R")
if (!dir.exists("results")) stopifnot(dir.create("results"))
set.seed(2025)

# ---- knobs -----------------------------------------------------------
K          <- 5000                       # replicates per cell
s_grid     <- seq(0.5, 1.2, by = 0.05)   # fine grid for Figure 2
tau0_T_vals <- c(100, 200, 400)
Vm_true_vals <- c(0.2, 2.0)
# Full grid: 3 x 2 x 15 x K cells. With K = 5000 this is ~450k pipeline
# runs; expect a long run. For a quick pass, lower K or coarsen s_grid.
# Parallel: library(furrr); plan(multisession); replace pmap_dfr below
# with furrr::future_pmap_dfr(..., .options = furrr_options(seed = TRUE)).

lambda0_T <- 0.01
delta_    <- 1.0
n_hist_   <- 400

grid <- tidyr::expand_grid(
  tau0_T  = tau0_T_vals,
  Vm_true = Vm_true_vals,
  s_vm    = s_grid
)

res_misspec <- purrr::pmap_dfr(grid, function(tau0_T, Vm_true, s_vm) {
  purrr::map_dfr(seq_len(K), function(iter) {
    run_pipeline_once_sep2(
      n_hist    = n_hist_,
      lambda0_H = 0.004,
      tau0_H    = 200,
      lambda0_T = lambda0_T,
      tau0_T    = tau0_T,
      delta     = delta_,
      Vm_hist   = Vm_true * s_vm,
      Vm_true   = Vm_true,
      beta_alt  = log(0.8),
      M0        = 1,
      alpha     = 0.05,
      power_target = 0.80,
      p1 = 0.5
    ) %>%
      dplyr::mutate(tau0_T = tau0_T, s_vm = s_vm, iter = iter)
  })
})

saveRDS(res_misspec, "results/03_misspec_raw.rds")

summary_misspec <- summarise_plan(res_misspec, tau0_T, Vm_true, s_vm)
write_csv(summary_misspec, "results/03_misspec_summary.csv")

# Table 4 rows (subset of the fine grid)
table4 <- summary_misspec %>%
  dplyr::filter(s_vm %in% c(0.5, 0.7, 0.8, 0.9, 1.1))
write_csv(table4, "results/03_misspec_table4.csv")
print(table4, n = Inf, width = Inf)