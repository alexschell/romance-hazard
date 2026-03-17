# R/01_load_data.R
# Load and validate raw relationship data from CSV.
#
# Supports two CSV formats (auto-detected by column names):
#
#   Format A — direct duration (current data):
#     duration_days : numeric, duration of relationship in days
#     event         : 1 = ended, 0 = censored (ongoing)
#
#   Format B — date pairs:
#     start_date    : parseable date string
#     end_date      : parseable date string (NA/empty = censored / ongoing)
#
# In either format, additional covariate columns are carried through.
#
# Output: data/processed/relationships.rds

library(dplyr)
library(lubridate)
library(here)

RAW_FILE    <- here("data", "raw", "relationships.csv")
OUTPUT_FILE <- here("data", "processed", "relationships.rds")

# --- Load ---
raw <- read.csv(RAW_FILE, stringsAsFactors = FALSE, na.strings = c("", "NA", "N/A"))

# --- Detect format and derive duration_days + event ---
if ("duration_days" %in% names(raw) && "event" %in% names(raw)) {
  # Format A: direct duration
  message("Detected format: duration_days + event")

  df <- raw |>
    mutate(
      duration_days = as.numeric(duration_days),
      event         = as.integer(event)
    )

  if (any(df$duration_days <= 0, na.rm = TRUE)) {
    warning(sprintf(
      "%d rows have duration_days <= 0; these may cause issues in log-scale models.",
      sum(df$duration_days <= 0, na.rm = TRUE)
    ))
  }

} else if ("start_date" %in% names(raw)) {
  # Format B: date pairs
  message("Detected format: start_date / end_date")

  df <- raw |>
    mutate(
      start_date         = ymd(start_date),
      end_date           = ymd(end_date),
      end_date_effective = if_else(is.na(end_date), today(), end_date),
      event              = as.integer(!is.na(end_date)),
      duration_days      = as.numeric(end_date_effective - start_date)
    ) |>
    filter(duration_days >= 0)

} else {
  stop(
    "CSV must contain either:\n",
    "  (A) columns 'duration_days' and 'event', or\n",
    "  (B) column 'start_date' (and optionally 'end_date')"
  )
}

# --- Add duration in months ---
df <- df |>
  mutate(
    duration_months = duration_days / 30.4375,
    id = row_number()
  )

message(sprintf(
  "Loaded %d relationships: %d ended, %d censored.",
  nrow(df), sum(df$event), sum(1L - df$event)
))
message(sprintf(
  "Duration range: %.0f – %.0f days (%.1f – %.1f months)",
  min(df$duration_days), max(df$duration_days),
  min(df$duration_months), max(df$duration_months)
))

# --- Save ---
dir.create(here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(df, OUTPUT_FILE)
message("Saved to ", OUTPUT_FILE)
