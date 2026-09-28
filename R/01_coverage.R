# ============================================================
# 01_coverage.R
# Finite-sample bias, SE accuracy, and CI coverage of Vm_hat.
# Part A: main grid under gamma frailty (NB marginal).
# Part B: frailty-robustness check -- the estimator targets
#   Vm = Var(m_i) without distributional assumptions, so the
#   Table 1 subset is repeated with lognormal and inverse-
#   Gaussian frailties matched to the same Vm.
# ============================================================

source("R/00_functions.R")
dir.create("results", showWarnings = FALSE)
set.seed(2025)

n_rep <- 5000

run_once <- function(n, lambda0, Vm, tau0, delta, frailty = "gamma"){
  dat <- r_counts_fast(n, lambda0, Vm, tau0, delta, frailty = frailty)
  est <- estimate_vm_orth_mom(dat$N, dat$tau)
  se  <- var_vm_sandwich_orth(dat$N, dat$tau, est$lambda_hat, est$Vm_hat)
  tibble(frailty = frailty,
         n = n, lambda0 = lambda0, Vm_true = Vm, tau0 = tau0, delta = delta,
         lambda_hat = est$lambda_hat, Vm_hat = est$Vm_hat,
         se_Vm = se$se_Vm)
}

simulate_grid <- function(grid, n_rep){
  map_dfr(seq_len(nrow(grid)), function(j){
    pars <- grid[j, ]
    map_dfr(seq_len(n_rep), ~ run_once(pars$n, pars$lambda0, pars$Vm,
                                       pars$tau0, pars$delta, pars$frailty))
  })
}

summarise_coverage <- function(sims){
  z <- qnorm(0.975)
  sims %>%
    mutate(
      ci_lo = pmax(Vm_hat - z * se_Vm, 0),   # only the lower limit truncated at 0
      ci_hi = Vm_hat + z * se_Vm,
      covered = (Vm_true >= ci_lo) & (Vm_true <= ci_hi)
    ) %>%
    group_by(frailty, tau0, lambda0, Vm_true, delta, n) %>%
    summarise(
      emp_mean    = mean(Vm_hat, na.rm = TRUE),
      emp_sd      = sd(Vm_hat, na.rm = TRUE),
      mean_est_se = mean(se_Vm, na.rm = TRUE),
      coverage    = round(mean(covered, na.rm = TRUE), 2),
      .groups = "drop"
    )
}

# ---- Part A: main grid, gamma frailty --------------------------------
grid_main <- expand.grid(
  tau0    = c(200, 400),
  lambda0 = c(0.002, 0.005),
  Vm      = c(0.01, 0.2, 0.8, 2.0),
  delta   = 0.1,            # Var(q) = delta^2/3
  n       = c(200, 400, 800, 1200),
  KEEP.OUT.ATTRS = FALSE
) %>% as_tibble() %>% mutate(frailty = "gamma")

sims_main <- simulate_grid(grid_main, n_rep)
saveRDS(sims_main, "results/01_coverage_raw.rds")

coverage_main <- summarise_coverage(sims_main)
write_csv(coverage_main, "results/01_coverage_summary.csv")
print(coverage_main, n = Inf)

# ---- Part B: frailty robustness on the Table 1 subset ----------------
# Same Vm, lognormal and inverse-Gaussian mixing.
grid_robust <- expand.grid(
  tau0    = 400,
  lambda0 = c(0.002, 0.005),
  Vm      = c(0.8, 2.0),
  delta   = 0.1,
  n       = c(400, 800),
  frailty = c("lognormal", "invgauss"),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>% as_tibble()

sims_robust <- simulate_grid(grid_robust, n_rep)
saveRDS(sims_robust, "results/01_coverage_frailty_robust_raw.rds")

coverage_robust <- summarise_coverage(sims_robust)
write_csv(coverage_robust, "results/01_coverage_frailty_robust_summary.csv")
print(coverage_robust, n = Inf)
