# Reproducible pipeline entry point -----------------------------------------
# First-time setup:
#   install.packages(c("renv", "targets", "tarchetypes"))
#   renv::restore()
# Run analysis:
#   targets::tar_make()

source(file.path("R", "00_setup.R"))
install_missing_packages()
load_packages()

library(targets)
library(tarchetypes)

tar_option_set(
  packages = required_packages,
  error = "stop",
  memory = "transient",
  garbage_collection = TRUE
)

list(
  tar_target(input_check, check_input_files(), cue = tar_cue(mode = "always")),
  tar_target(data_preparation, {
    source(file.path("R", "01_data_preparation.R"), local = TRUE)
    TRUE
  }),
  tar_target(flavonoid_intake_construction, {
    source(file.path("R", "02_flavonoid_intake_construction.R"), local = TRUE)
    TRUE
  }),
  tar_target(survival_and_covariates, {
    source(file.path("R", "03_survival_and_covariates.R"), local = TRUE)
    TRUE
  }),
  tar_target(food_source_and_spline_analysis, {
    source(file.path("R", "04_food_source_and_spline_analysis.R"), local = TRUE)
    TRUE
  }),
  tar_target(cox_models_tables_interactions, {
    source(file.path("R", "05_cox_models_tables_interactions.R"), local = TRUE)
    TRUE
  })
)
