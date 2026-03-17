# R/03_km_analysis.R
# Kaplan-Meier survival estimation and non-parametric summary.
#
# Produces:
#   - Overall KM curve
#   - Nelson-Aalen cumulative hazard
#   - Saved figures to output/figures/

library(survival)
# library(ggsurvfit)  # install from CRAN when online; using survminer instead
library(ggplot2)
library(dplyr)
library(here)

source(here("R", "05_viz_templates.R"))

SURV_FILE  <- here("data", "processed", "surv_df.rds")
FIG_DIR    <- here("output", "figures")

df <- readRDS(SURV_FILE)

# ── Overall KM fit ─────────────────────────────────────────────────────────────
km_fit <- survfit(Surv(duration_months, event) ~ 1, data = df)

summary_km <- summary(km_fit, times = c(3, 6, 12, 24, 36, 60))
cat("\n── KM survival estimates at key timepoints (months) ──\n")
print(data.frame(
  months   = summary_km$time,
  survival = round(summary_km$surv, 3),
  lower    = round(summary_km$lower, 3),
  upper    = round(summary_km$upper, 3),
  n_risk   = summary_km$n.risk,
  n_event  = summary_km$n.event
))

# Median survival time
cat(sprintf(
  "\nMedian relationship duration: %.1f months (95%% CI: %.1f – %.1f)\n",
  summary(km_fit)$table["median"],
  summary(km_fit)$table["0.95LCL"],
  summary(km_fit)$table["0.95UCL"]
))

# ── KM curve plot ──────────────────────────────────────────────────────────────
p_km <- plot_km(
  fit   = km_fit,
  title = "Relationship Survival",
  xlab  = "Duration (months)"
)
save_fig(p_km, "km_curve.png", width = 8, height = 5)

# ── Nelson-Aalen cumulative hazard ─────────────────────────────────────────────
na_fit <- survfit(Surv(duration_months, event) ~ 1, data = df, type = "fleming")

p_cumhaz <- plot_cumhaz(
  fit   = na_fit,
  title = "Cumulative Hazard (Nelson-Aalen)",
  xlab  = "Duration (months)"
)
save_fig(p_cumhaz, "cumhaz_na.png", width = 8, height = 5)

message("\nKM analysis complete. Figures saved to ", FIG_DIR)
