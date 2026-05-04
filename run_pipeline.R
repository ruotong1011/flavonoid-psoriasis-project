# Run this file from the project root
install.packages(c("renv", "targets", "tarchetypes"), repos = "https://cloud.r-project.org")
renv::restore(prompt = FALSE)
targets::tar_make()
