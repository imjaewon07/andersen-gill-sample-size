# ============================================================
# make_figures_base.R  (v3)
# Figures 1-4 in base R graphics. Color plus lty/pch (B/W-safe).
# Changes from v2:
#   - Fig 2: top row = mean planned sample size (oracle reference
#     lines 1859 / 346); smooth lines (no point markers).
#   - Fig 3: interim Vm-hat row removed -> 2 x 2 layout.
# Run from the project root (the folder containing results/).
# Output: figures/figure{1,2,3,4}_260903.pdf  (vector)
# ============================================================

if (!dir.exists("figures")) dir.create("figures")

COL <- c("#4477AA", "#BB5566", "#228833")   # plug-in, +0.5 SE, +1.0 SE
LTY <- c(1, 2, 4)
PCH <- c(16, 17, 15)
LEG <- c("Plug-in", "+0.5 SE", "+SE")

panel3 <- function(x, Y, ylab, ylim = NULL, href = NULL, vref = NULL,
                   main = NULL, xlab = "", xaxt_at = NULL, type = "b") {
  ## Y: matrix with 3 columns (plug-in, +0.5 SE, +1.0 SE)
  if (is.null(ylim)) ylim <- range(Y, href, na.rm = TRUE)
  matplot(x, Y, type = type, lty = LTY, pch = PCH, col = COL,
          lwd = if (type == "l") 1.6 else 1.3,
          cex = 0.8, xlab = xlab, ylab = ylab, ylim = ylim,
          main = main, xaxt = if (is.null(xaxt_at)) "s" else "n",
          cex.main = 1.0, font.main = 1, las = 1)
  if (!is.null(xaxt_at)) axis(1, at = xaxt_at)
  if (!is.null(href)) abline(h = href, lty = 3, col = "grey45")
  if (!is.null(vref)) abline(v = vref, lty = 3, col = "grey45")
  box()
}

bottom_legend <- function() {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
      new = TRUE)
  plot.new()
  legend("bottom", LEG, lty = LTY, pch = PCH, col = COL,
         horiz = TRUE, bty = "n", lwd = 1.3, cex = 0.95)
}

## ------------------------------------------------------------
## Figure 1 -- Study 2: reliability and planned n vs tau_T
## ------------------------------------------------------------
b <- read.csv("results/02_baseline_summary.csv")

pdf("figures/figure1_260903.pdf", width = 8, height = 6.4)
par(mfrow = c(2, 2), mar = c(4, 4.4, 2.5, 1), oma = c(2.5, 0, 0, 0))
for (vt in c(2.0, 0.2)) {
  d <- b[b$Vm_true == vt, ]
  d <- d[order(d$tau0_T), ]
  panel3(d$tau0_T,
         cbind(d$power_prob75_plan, d$power_prob75_plusHalfSE,
               d$power_prob75_plusSE),
         ylab = expression(Pr(power >= 0.75)),
         ylim = c(0.6, 1), href = 0.90, xaxt_at = c(100, 200, 400),
         xlab = expression(tau[T]),
         main = bquote(V[m]^{(T)} == .(vt)))
}
for (vt in c(2.0, 0.2)) {
  d <- b[b$Vm_true == vt, ]
  d <- d[order(d$tau0_T), ]
  panel3(d$tau0_T,
         cbind(d$n_mean, d$n_plusHalfSE_mean, d$n_plusSE_mean),
         ylab = "Mean planned sample size",
         xaxt_at = c(100, 200, 400), xlab = expression(tau[T]))
}
bottom_legend()
dev.off()

## ------------------------------------------------------------
## Figure 2 -- Study 3: misspecification (tau_T = 400, fine s grid)
## Top row: mean planned sample size (oracle lines 1859 / 346).
## Bottom row: Pr(power >= 0.75) with the 0.90 criterion line.
## Smooth lines (type = "l"): the s grid is fine (step 0.05).
## ------------------------------------------------------------
m <- read.csv("results/03_misspec_summary.csv")
m <- m[m$tau0_T == 400, ]
oracle_n <- c("2" = 1859, "0.2" = 346)

pdf("figures/figure2_260903.pdf", width = 8, height = 6.4)
par(mfrow = c(2, 2), mar = c(4, 4.6, 2.5, 1), oma = c(2.5, 0, 0, 0))
for (vt in c(2.0, 0.2)) {
  d <- m[m$Vm_true == vt, ]
  d <- d[order(d$s_vm), ]
  panel3(d$s_vm,
         cbind(d$n_mean, d$n_plusHalfSE_mean, d$n_plusSE_mean),
         ylab = "Mean planned sample size",
         href = oracle_n[as.character(vt)], vref = 1,
         xlab = expression(s[V[m]]), type = "l",
         main = bquote(V[m]^{(T)} == .(vt)))
}
for (vt in c(2.0, 0.2)) {
  d <- m[m$Vm_true == vt, ]
  d <- d[order(d$s_vm), ]
  panel3(d$s_vm,
         cbind(d$power_prob75_plan, d$power_prob75_plusHalfSE,
               d$power_prob75_plusSE),
         ylab = expression(Pr(power >= 0.75)),
         ylim = c(0, 1), href = 0.90, vref = 1,
         xlab = expression(s[V[m]]), type = "l")
}
bottom_legend()
dev.off()

## ------------------------------------------------------------
## Figure 3 -- Study 4: interim re-estimation vs r  (2 x 2)
## rows: mean final n (ref 1859), Pr(>=0.75) (ref 0.90);
## cols: s = 0.5, 0.7
## ------------------------------------------------------------
it <- read.csv("results/04_interim_summary.csv")
it <- it[order(it$s_vm, it$r_frac), ]

pdf("figures/figure3_260903.pdf", width = 8, height = 6.4)
par(mfrow = c(2, 2), mar = c(4, 4.6, 2.5, 1), oma = c(2.5, 0, 0, 0))
for (s in c(0.5, 0.7)) {
  d <- it[it$s_vm == s, ]
  panel3(d$r_frac,
         cbind(d$n_final_mean, d$n_final_plusHalfSE_mean,
               d$n_final_plusSE_mean),
         ylab = "Mean final sample size",
         ylim = c(1750, 2300), href = 1859,
         xlab = "Information fraction r",
         main = bquote(s[V[m]] == .(s)))
}
for (s in c(0.5, 0.7)) {
  d <- it[it$s_vm == s, ]
  panel3(d$r_frac,
         cbind(d$power_prob75_final, d$power_prob75_plusHalfSE,
               d$power_prob75_plusSE),
         ylab = expression(Pr(power >= 0.75)),
         ylim = c(0.55, 1), href = 0.90,
         xlab = "Information fraction r")
}
bottom_legend()
dev.off()

## ------------------------------------------------------------
## Figure 4 -- Study 5: empirical type I error with MC error bars
## `off` sets the spacing between the three rules WITHIN each
## r cluster.
## ------------------------------------------------------------
t1 <- read.csv("results/05_type1_summary.csv")
t1 <- t1[t1$s_vm %in% c(0.5, 0.7), ]

pdf("figures/figure4_260903.pdf", width = 8, height = 4.2)
par(mfrow = c(1, 2), mar = c(4, 4.4, 2.5, 1), oma = c(2.5, 0, 0, 0))
off <- c(-0.035, 0, 0.035)
for (s in c(0.5, 0.7)) {
  d <- t1[t1$s_vm == s, ]
  d <- d[order(d$r_frac), ]
  est <- cbind(d$type1_main, d$type1_plusHalfSE, d$type1_plusSE)
  se  <- cbind(d$mc_se_main, d$mc_se_plusHalfSE, d$mc_se_plusSE)
  plot(NULL, xlim = c(0.03, 0.97), ylim = c(0.040, 0.062),
       xlab = "Information fraction r", ylab = "Empirical type I error",
       main = bquote(s[V[m]] == .(s)), cex.main = 1.0, font.main = 1,
       las = 1, xaxt = "n")
  axis(1, at = c(0.1, 0.3, 0.5, 0.7, 0.9))
  abline(h = 0.05, lty = 3, col = "grey45")
  for (j in 1:3) {
    x <- d$r_frac + off[j]
    arrows(x, est[, j] - 1.96 * se[, j], x, est[, j] + 1.96 * se[, j],
           angle = 90, code = 3, length = 0.02, col = COL[j], lwd = 0.9)
    points(x, est[, j], pch = PCH[j], cex = 0.9, col = COL[j])
  }
  box()
}
bottom_legend()
dev.off()

message("Figures written to figures/ (PDF, vector)")