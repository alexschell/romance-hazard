# R/04_flexible_hazard.R
# Flexible hazard rate estimation via two complementary approaches:
#
#   A. mgcv::gam() — "Poisson trick" on interval-expanded counting process data
#      Treats event counts as Poisson with log(exposure) offset; models log-hazard
#      as an arbitrary smooth function of time using penalised regression splines.
#      Pros: maximum flexibility (any mgcv smooth), easy to add covariates.
#      Cons: requires manual interval expansion (done in 02_preprocess.R).
#
#   B. rstpm2::stpm2() — Royston-Parmar flexible parametric survival model
#      (skipped: rstpm2 not available in offline environment; install from CRAN when online)
#
# Note: run A only for now.
#
# Produces:
#   - Summary output to console
#   - output/figures/hazard_gam.png
#   - output/figures/survival_gam.png
#   - output/figures/hazard_stpm2.png
#   - output/figures/survival_stpm2.png

library(mgcv)
# library(rstpm2)  # install from CRAN when online
# library(gratia)  # install from CRAN when online
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

source(here("R", "05_viz_templates.R"))

SURV_FILE     <- here("data", "processed", "surv_df.rds")
EXPANDED_FILE <- here("data", "processed", "expanded.rds")
FIG_DIR       <- here("output", "figures")

surv_df  <- readRDS(SURV_FILE)
expanded <- readRDS(EXPANDED_FILE)

# ══════════════════════════════════════════════════════════════════════════════
# A. GAM / Poisson trick (mgcv)
# ══════════════════════════════════════════════════════════════════════════════
# Model: log h(t) = f(t)   where f() is a penalised cubic regression spline
# The offset log(exposure) converts from event-count scale to rate scale.

cat("\n── Fitting GAM flexible hazard (mgcv) ──\n")

gam_fit <- mgcv::gam(
  event ~ s(t_mid, bs = "cr", k = 10) + offset(log(exposure)),
  family = poisson(link = "log"),
  data   = expanded,
  method = "REML"   # REML smoothness selection
)

summary(gam_fit)

# Check effective degrees of freedom of the smooth (higher = more wiggly)
cat(sprintf("\nSmooth edf: %.2f  (k = 10; increase k if edf ≈ k - 1)\n",
            summary(gam_fit)$edf))

# ── Predict hazard over a fine time grid ──
t_grid <- seq(0, max(expanded$t_end), length.out = 200)
pred_gam <- data.frame(t_mid = t_grid, exposure = 1)  # exposure = 1 → hazard rate

pred_gam <- pred_gam |>
  mutate(
    fit   = predict(gam_fit, newdata = pred_gam, type = "response"),
    se    = predict(gam_fit, newdata = pred_gam, type = "response", se.fit = TRUE)$se.fit,
    lower = pmax(fit - 1.96 * se, 0),
    upper = fit + 1.96 * se,
    time  = t_grid
  )

# ── Derive survival from GAM hazard (numerical integration) ──
# S(t) = exp(-∫₀ᵗ h(u) du)  approximated by cumulative trapezoid sum
dt <- diff(c(0, pred_gam$time))
pred_gam <- pred_gam |>
  mutate(
    cum_hazard = cumsum(fit * dt),
    survival   = exp(-cum_hazard)
  )

p_haz_gam <- plot_hazard(
  t        = pred_gam$time,
  hazard   = pred_gam$fit,
  lower    = pred_gam$lower,
  upper    = pred_gam$upper,
  title    = "Flexible Hazard Rate (GAM, Poisson trick)",
  xlab     = "Duration (months)"
)
save_fig(p_haz_gam, "hazard_gam.png", width = 8, height = 5)

p_surv_gam <- plot_survival(
  t        = pred_gam$time,
  survival = pred_gam$survival,
  title    = "Survival Function (derived from GAM hazard)",
  xlab     = "Duration (months)"
)
save_fig(p_surv_gam, "survival_gam.png", width = 8, height = 5)

# Optional: visualise GAM smooth directly (log-hazard scale)
# gratia::draw(gam_fit) — available once gratia installed from CRAN

# ══════════════════════════════════════════════════════════════════════════════
# B. Royston-Parmar (rstpm2) — skipped, install rstpm2 from CRAN when online
# ══════════════════════════════════════════════════════════════════════════════

message("\nFlexible hazard analysis complete. Figures saved to ", FIG_DIR)
