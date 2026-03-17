# R/05_viz_templates.R
# Reusable ggplot2 chart functions for survival analysis outputs.
# Source this file at the top of any analysis script.
#
# Note: uses survminer for KM curves (ggsurvfit not available offline).
#       gratia not available; GAM smooths plotted manually.
#
# Functions:
#   theme_rh()          — project ggplot2 theme
#   save_fig()          — ggsave wrapper with consistent defaults
#   plot_km()           — Kaplan-Meier survival curve (survminer)
#   plot_cumhaz()       — Nelson-Aalen cumulative hazard curve
#   plot_hazard()       — Smooth hazard rate curve with CI band
#   plot_survival()     — Parametric survival curve with CI band

library(ggplot2)
library(survminer)
library(survival)
library(here)

FIG_DIR <- here("output", "figures")

# ── Project theme ──────────────────────────────────────────────────────────────
theme_rh <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = element_text(colour = "grey40"),
      axis.title       = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

# ── Save helper ────────────────────────────────────────────────────────────────
save_fig <- function(plot, filename, width = 8, height = 5, dpi = 150) {
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(FIG_DIR, filename)
  if (inherits(plot, "ggsurvplot")) {
    # ggsurvplot objects must be printed to a graphics device directly
    png(path, width = width, height = height, units = "in", res = dpi)
    print(plot)
    dev.off()
  } else {
    ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  }
  message("Saved: ", path)
  invisible(path)
}

# ── Kaplan-Meier curve (survminer) ─────────────────────────────────────────────
plot_km <- function(fit, title = "Kaplan-Meier Survival Curve",
                    xlab = "Time", conf_int = TRUE) {
  ggsurvplot(
    fit,
    conf.int       = conf_int,
    risk.table     = TRUE,
    risk.table.height = 0.25,
    xlab           = xlab,
    ylab           = "Survival probability",
    title          = title,
    palette        = "#2166AC",
    ggtheme        = theme_rh(),
    risk.table.fontsize = 3
  )
}

# ── Cumulative hazard curve ────────────────────────────────────────────────────
plot_cumhaz <- function(fit, title = "Cumulative Hazard", xlab = "Time") {
  km_df <- data.frame(
    time   = fit$time,
    cumhaz = fit$cumhaz,
    lower  = pmax(fit$cumhaz - 1.96 * fit$std.err, 0),
    upper  = fit$cumhaz + 1.96 * fit$std.err
  )

  ggplot(km_df, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#D6604D", alpha = 0.15) +
    geom_step(aes(y = cumhaz), colour = "#D6604D", linewidth = 0.9) +
    labs(title = title, x = xlab, y = "Cumulative hazard H(t)") +
    theme_rh()
}

# ── Smooth hazard rate curve ───────────────────────────────────────────────────
plot_hazard <- function(t, hazard, lower = NULL, upper = NULL,
                        title = "Hazard Rate", xlab = "Time") {
  df <- data.frame(t = t, hazard = hazard)

  p <- ggplot(df, aes(x = t, y = hazard))

  if (!is.null(lower) && !is.null(upper)) {
    df$lower <- lower
    df$upper <- upper
    p <- p + geom_ribbon(data = df,
                         aes(ymin = pmax(lower, 0), ymax = upper),
                         fill = "#4DAC26", alpha = 0.15)
  }

  p +
    geom_line(colour = "#4DAC26", linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    labs(title = title, x = xlab, y = "Hazard rate h(t)") +
    theme_rh()
}

# ── Parametric survival curve ──────────────────────────────────────────────────
plot_survival <- function(t, survival, lower = NULL, upper = NULL,
                          title = "Survival Function", xlab = "Time") {
  df <- data.frame(t = t, survival = survival)

  p <- ggplot(df, aes(x = t, y = survival))

  if (!is.null(lower) && !is.null(upper)) {
    df$lower <- lower
    df$upper <- upper
    p <- p + geom_ribbon(data = df,
                         aes(ymin = pmax(lower, 0), ymax = pmin(upper, 1)),
                         fill = "#762A83", alpha = 0.15)
  }

  p +
    geom_line(colour = "#762A83", linewidth = 0.9) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    labs(title = title, x = xlab, y = "Survival probability S(t)") +
    theme_rh()
}
