# R/01_load_data.R
# Load and validate raw relationship data from CSV.
#
# Expected CSV columns (adjust COLUMN_MAP below to match your actual headers):
#   - start_date  : date relationship began (parseable date string)
#   - end_date    : date relationship ended (NA or empty = censored / ongoing)
#   - [optional]  : any covariate columns to carry through
#
# Output: data/processed/relationships.rds

library(dplyr)
library(lubridate)
library(here)

# --- Configuration ---
RAW_FILE    <- here("data", "raw", "relationships.csv")
OUTPUT_FILE <- here("data", "processed", "relationships.rds")

# Map your actual CSV column names to the canonical names used downstream.
# Edit the right-hand side to match your headers.
COLUMN_MAP <- c(
  start_date = "start_date",
  end_date   = "end_date"
  # add covariates here, e.g.:
  # met_online = "met_online"
)

# --- Load ---
raw <- read.csv(RAW_FILE, stringsAsFactors = FALSE, na.strings = c("", "NA", "N/A"))

# --- Rename to canonical names ---
raw <- raw |>
  rename(any_of(setNames(COLUMN_MAP, names(COLUMN_MAP))))

# --- Parse dates ---
raw <- raw |>
  mutate(
    start_date = ymd(start_date),
    end_date   = ymd(end_date)
  )

# --- Validate ---
problems <- list()

if (any(is.na(raw$start_date))) {
  problems <- c(problems, sprintf(
    "%d rows have unparseable start_date", sum(is.na(raw$start_date))
  ))
}

future_starts <- raw$start_date > today()
if (any(future_starts, na.rm = TRUE)) {
  problems <- c(problems, sprintf(
    "%d rows have start_date in the future", sum(future_starts, na.rm = TRUE)
  ))
}

end_before_start <- !is.na(raw$end_date) & raw$end_date < raw$start_date
if (any(end_before_start, na.rm = TRUE)) {
  problems <- c(problems, sprintf(
    "%d rows have end_date before start_date", sum(end_before_start, na.rm = TRUE)
  ))
}

if (length(problems) > 0) {
  warning("Data quality issues:\n", paste("-", problems, collapse = "\n"))
}

# --- Derive core fields ---
df <- raw |>
  mutate(
    # Use today as the censoring date for ongoing relationships
    end_date_effective = if_else(is.na(end_date), today(), end_date),
    # event = 1 if relationship ended (not censored)
    event = as.integer(!is.na(end_date)),
    # Duration in days and months
    duration_days   = as.numeric(end_date_effective - start_date),
    duration_months = duration_days / 30.4375
  ) |>
  filter(duration_days >= 0)

message(sprintf(
  "Loaded %d relationships: %d ended, %d censored (ongoing).",
  nrow(df), sum(df$event), sum(1 - df$event)
))

# --- Save ---
saveRDS(df, OUTPUT_FILE)
message("Saved to ", OUTPUT_FILE)
