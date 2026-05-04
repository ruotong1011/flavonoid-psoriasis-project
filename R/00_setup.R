# Project setup -------------------------------------------------------------
# This file is sourced by _targets.R before the pipeline is run.

options(stringsAsFactors = FALSE)

required_packages <- c(
  "dplyr", "tidyr", "stringr", "purrr", "tibble", "readr", "survival",
  "tableone", "knitr", "ggplot2", "lubridate", "ggsurvfit", "forestmodel",
  "mediation", "rms", "broom", "targets", "tarchetypes"
)

install_missing_packages <- function(pkgs = required_packages) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

load_packages <- function(pkgs = required_packages) {
  invisible(lapply(pkgs, library, character.only = TRUE))
}

paths <- list(
  data_raw       = file.path("data", "raw"),
  data_processed = file.path("data", "processed"),
  tables         = file.path("results", "tables"),
  figures        = file.path("results", "figures")
)

for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# Input file names expected in data/raw. Do not commit UK Biobank raw data.
input_files <- list(
  participants_recall_over1 = file.path(paths$data_raw, "participants.recall.over1.csv"),
  participants_recall_over2 = file.path(paths$data_raw, "participants.recall.over2.csv"),
  demographics              = file.path(paths$data_raw, "demographics.csv"),
  incidence                 = file.path(paths$data_raw, "4907incidence.csv"),
  control                   = file.path(paths$data_raw, "484895control.csv"),
  average_24h_foods          = file.path(paths$data_raw, "average_24hrecall_individual_foods.csv"),
  flavonoid_database        = file.path(paths$data_raw, "flavonoid_ukbiobank.csv")
)

check_input_files <- function(files = input_files) {
  missing <- files[!file.exists(unlist(files))]
  if (length(missing) > 0) {
    stop(
      "Missing input files in data/raw. See data/README_data.md. Missing: ",
      paste(names(missing), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
  invisible(file)
}
