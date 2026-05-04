# Minimal renv activation file placeholder.
# After cloning, run: install.packages("renv"); renv::init(bare = TRUE); renv::snapshot()
if (requireNamespace("renv", quietly = TRUE)) {
  try(renv::load(), silent = TRUE)
}
