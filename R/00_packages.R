# R/00_packages.R
# Run this script ONCE after cloning to bootstrap the renv environment.
#
# Steps:
#   1. Install renv if not present
#   2. renv::init() creates .Rprofile, renv/, and renv.lock
#   3. Install project packages
#   4. renv::snapshot() locks versions
#
# After this, collaborators just need: renv::restore()

# --- Bootstrap renv ---
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::init()

# --- Core survival analysis ---
install.packages("survival")    # KM, Cox PH, counting process Surv()
install.packages("rstpm2")      # Royston-Parmar flexible parametric models
install.packages("flexsurv")    # Additional flexible parametric families

# --- Flexible hazard via GAM (Poisson trick) ---
install.packages("mgcv")        # GAMs; base R but listed explicitly for renv
install.packages("gratia")      # GAM visualization: smooth terms, derivatives

# --- Visualization ---
install.packages("ggplot2")
install.packages("ggsurvfit")   # ggplot2-native survival curves (modern survminer alt)

# --- Data wrangling ---
install.packages("dplyr")
install.packages("tidyr")
install.packages("lubridate")   # date arithmetic for duration computation
install.packages("here")        # relative paths anchored to project root

# --- Lock environment ---
renv::snapshot()

message("Environment ready. Commit renv.lock to version control.")
