# R/05_viz_templates.R
# Reusable ggplot2 chart functions for survival analysis outputs.
# Source this file at the top of any analysis script.
#
# Functions:
#   theme_rh()          — project ggplot2 theme
#   save_fig()          — ggsave wrapper with consistent defaults
#   plot_km()           — Kaplan-Meier survival curve with CI ribbon
#   plot_cumhaz()       — Nelson-Aalen cumulative hazard curve
#   plot_hazard()       — Smooth hazard rate curve with CI band
#   plot_survival()     — Parametric survival curve with CI band

library(ggplot2)
library(ggsurvfit)
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
#' Save a ggplot to output/figures/ with consistent DPI.
save_fig <- function(plot, filename, width = 8, height = 5, dpi = 150) {
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(FIG_DIR, filename)
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved: ", path)
  invisible(path)
}

# ── Kaplan-Meier curve ─────────────────────────────────────────────────────────
#' Plot a KM survival curve using ggsurvfit.
#'
#' @param fit    survfit object (from survival::survfit())
#' @param title  plot title
#' @param xlab   x-axis label (default "Time")
#' @param conf_int logical — show confidence interval ribbon
plot_km <- function(fit, title = "Kaplan-Meier Survival Curve",
                    xlab = "Time", conf_int = TRUE) {
  p <- survfit2(fit) |>
    ggsurvfit(linewidth = 0.9, color = "#2166AC") +
    labs(title = title, x = xlab, y = "Survival probability") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
    theme_rh()

  if (conf_int) {
    p <- p + add_confidence_interval(fill = "#2166AC", alpha = 0.15)
  }

  p <- p + add_risktable(
    risktable_stats = c("n.risk", "n.event"),
    theme = theme_risktable_default(axis.text.y.size = 9)
  )

  p
}

# ── Cumulative hazard curve ────────────────────────────────────────────────────
#' Plot Nelson-Aalen cumulative hazard.
#'
#' @param fit    survfit object fitted with type = "fleming"
#' @param title  plot title
#' @param xlab   x-axis label
plot_cumhaz <- function(fit, title = "Cumulative Hazard",
                         xlab = "Time") {
  km_df <- data.frame(
    time   = fit$time,
    cumhaz = fit$cumhaz,
    lower  = fit$cumhaz - 1.96 * fit$std.err,
    upper  = fit$cumhaz + 1.96 * fit$std.err
  )

  ggplot(km_df, aes(x = time)) +
    geom_ribbon(aes(ymin = pmax(lower, 0), ymax = upper),
                fill = "#D6604D", alpha = 0.15) +
    geom_step(aes(y = cumhaz), colour = "#D6604D", linewidth = 0.9) +
    labs(title = title, x = xlab, y = "Cumulative hazard H(t)") +
    theme_rh()
}

# ── Smooth hazard rate curve ───────────────────────────────────────────────────
#' Plot an estimated smooth hazard rate h(t) with confidence band.
#'
#' @param t      numeric vector of time points
#' @param hazard numeric vector of hazard estimates
#' @param lower  lower CI bound (optional)
#' @param upper  upper CI bound (optional)
#' @param title  plot title
#' @param xlab   x-axis label
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
#' Plot a model-estimated survival curve S(t) with confidence band.
#'
#' @param t        numeric vector of time points
#' @param survival numeric vector of survival estimates
#' @param lower    lower CI bound (optional)
#' @param upper    upper CI bound (optional)
#' @param title    plot title
#' @param xlab     x-axis label
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
