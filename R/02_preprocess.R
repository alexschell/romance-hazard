# R/02_preprocess.R
# Transform relationship duration data into two forms:
#   1. Standard survival format (one row per relationship) — for KM and rstpm2
#   2. Counting process / interval-expanded format — for GAM-based flexible hazard
#
# The "Poisson trick": expand each relationship into monthly intervals, then
# model event counts with a GAM + Poisson family and log(exposure) offset.
# This lets mgcv::gam() estimate an arbitrarily flexible baseline hazard h(t).
#
# Input:  data/processed/relationships.rds
# Output: data/processed/surv_df.rds       (standard survival, one row per case)
#         data/processed/expanded.rds      (counting process, one row per interval)

library(dplyr)
library(tidyr)
library(survival)
library(here)

INPUT_FILE    <- here("data", "processed", "relationships.rds")
SURV_FILE     <- here("data", "processed", "surv_df.rds")
EXPANDED_FILE <- here("data", "processed", "expanded.rds")

# Interval width in months for the Poisson expansion
INTERVAL_WIDTH <- 1

df <- readRDS(INPUT_FILE)

# ── 1. Standard survival format ────────────────────────────────────────────────
# Attach a survival object for convenience; downstream scripts use raw columns.
surv_df <- df |>
  mutate(
    surv_obj = Surv(duration_months, event)
  )

saveRDS(surv_df, SURV_FILE)
message("Saved standard survival data: ", nrow(surv_df), " rows → ", SURV_FILE)

# ── 2. Counting process / interval expansion ────────────────────────────────────
# For each relationship, create one row per monthly interval up to event/censoring.
# Each row records:
#   t_start  — interval left endpoint (months)
#   t_end    — interval right endpoint (months)
#   t_mid    — midpoint (used as the time covariate in GAM)
#   exposure — actual time at risk within interval (< 1 for the final interval)
#   event    — 1 if the relationship ended within this interval, else 0

expand_to_intervals <- function(row, width = INTERVAL_WIDTH) {
  dur   <- row$duration_months
  ev    <- row$event

  # Interval breakpoints covering [0, dur]
  breaks <- seq(0, ceiling(dur / width) * width, by = width)
  n      <- length(breaks) - 1

  intervals <- tibble(
    id      = row$id,
    t_start = breaks[seq_len(n)],
    t_end   = breaks[seq_len(n) + 1]
  ) |>
    mutate(
      # Clip the final interval to the actual duration
      t_end    = pmin(t_end, dur),
      exposure = t_end - t_start,
      t_mid    = (t_start + t_end) / 2,
      # Event only occurs in the last interval AND only if not censored
      event    = 0L
    )

  if (ev == 1) {
    intervals$event[nrow(intervals)] <- 1L
  }

  # Drop zero-width intervals (can arise if duration is an exact multiple of width)
  intervals |> filter(exposure > 0)
}

# Add a row id if not present
if (!"id" %in% names(df)) {
  df <- df |> mutate(id = row_number())
}

expanded <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
  expand_to_intervals(df[i, ])
}))

# Carry covariate columns from the original data
covariate_cols <- setdiff(names(df), c("id", "start_date", "end_date",
                                        "end_date_effective", "event",
                                        "duration_days", "duration_months",
                                        "surv_obj"))
if (length(covariate_cols) > 0) {
  expanded <- expanded |>
    left_join(df |> select(id, all_of(covariate_cols)), by = "id")
}

saveRDS(expanded, EXPANDED_FILE)
message(sprintf(
  "Saved counting-process data: %d intervals from %d relationships → %s",
  nrow(expanded), n_distinct(expanded$id), EXPANDED_FILE
))

# Sanity check: total events should match original
stopifnot(sum(expanded$event) == sum(df$event))
message("Sanity check passed: event counts match.")
