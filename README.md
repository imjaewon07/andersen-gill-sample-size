# Simulation code

Supplementary code for "Sample size calculation for Andersen–Gill models
with estimated between-subject heterogeneity."

## Layout

```
R/00_functions.R      shared functions (DGP, orthogonalized MOM estimator,
                      sandwich SE, planning formula, Wald test, pipelines)
R/01_coverage.R       Study 1/2: bias, SE accuracy, CI coverage of Vm_hat
                      (Part A: gamma frailty grid; Part B: lognormal and
                      inverse-Gaussian frailty robustness on the Table 1 subset)
R/02_baseline_plan.R  baseline planning: correctly specified historical Vm
R/03_misspec.R        Study 3: historical misspecification (Vm_hist = s Vm_true)
R/04_interim.R        Study 4: interim control-arm re-estimation
                      (reads results/03_misspec_raw.rds; Table 5 uses s = 0.5, 0.7)
R/05_type1_error.R    Study 5: empirical type I error with linked
                      interim/final data
results/              .rds (raw) and .csv (summary) outputs
```

## Run order

From the project root:

```r
source("R/01_coverage.R")
source("R/02_baseline_plan.R")
source("R/03_misspec.R")      # must run before 04
source("R/04_interim.R")
source("R/05_type1_error.R")  # slowest; see parallel note in the script
```

More scenarios are simulated than are reported; the manuscript tables and
figures use subsets of the summary CSVs.

## Notes

- All randomness is controlled by `set.seed(2025)` at the top of each script.
  Scripts are independent given the stated dependency (03 -> 04).
- `final_wald_test()` is a closed-form equivalent of
  `glm(N ~ x + offset(log(tau)), family = poisson)` with an HC0 sandwich;
  the equivalence is exact because the two-group model is saturated.
- `var_vm_sandwich_orth(df_correct = TRUE)` applies an optional
  sqrt(n/(n-1)) small-sample correction (off by default, matching the
  manuscript).
- Frailty distributions: the DGP supports gamma (NB marginal, sampled
  directly), lognormal (sigma^2 = log(1+Vm)), and inverse-Gaussian
  (shape = 1/Vm) mixing, all with E[m]=1 and Var(m)=Vm; the inverse-
  Gaussian sampler is self-contained (Michael-Schucany-Haas), so no
  extra packages are needed.
- Requires: tidyverse. Optional for parallel runs: furrr.
