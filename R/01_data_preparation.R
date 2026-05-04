# 01_data_preparation.R
# Auto-refactored from the original analysis notebook: github.Rmd.
# Run via targets::tar_make() from the project root.
# NOTE: UK Biobank raw data are not included in this repository.

source(file.path("R", "00_setup.R"))
load_packages()

# ---- Rmd chunk ----


# Load necessary libraries
library(psych)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)
library(tableone)
library(knitr)
library(ggplot2)
library(lubridate)
library(ggsurvfit)
library(forestmodel)
library(lubridate)
library(mediation)

# ---- Rmd chunk ----
# Load and filter psoriasis demographic data
participants.recall.over1 <- read.csv(input_files$participants_recall_over1)
demographics <- read.csv(input_files$demographics)
filtered.incidence <- read.csv(input_files$incidence)
p.demographics <- demographics %>%
  filter(Participant.ID %in% filtered.incidence$Participant.ID) %>%
  mutate(group = "psoriasis") # Filter participants based on incidence data

# Load and filter non-psoriasis demographic data
demographics <- read.csv(input_files$demographics)
filtered.control <- read.csv(input_files$control)
n.demographics <- demographics %>%
  filter(Participant.ID %in% filtered.control$Participant.ID) %>%
  mutate(group = "non-psoriasis") # Filter participants based on control data

# Combine psoriasis and non-psoriasis demographic data
demographics <- rbind(p.demographics, n.demographics) %>%
  filter(Participant.ID %in% participants.recall.over1$Participant.ID)
demographics <- demographics[order(demographics$Participant.ID), ]

# Load non-psoriasis 24-hour recall data
p.24h.foods.average <- read.csv(input_files$average_24h_foods) %>% 
  filter(Participant.ID %in% demographics$Participant.ID)
p.24h.foods.average[, 1:201][is.na(p.24h.foods.average[, 1:201])] <- 0

p.24h.foods.average <- p.24h.foods.average[, c(2:201, 1)]

p.24h.foods.average <- cbind(
  p.24h.foods.average[, 1:139], 
  Vegetarian.alternatives.intake = 0, 
  p.24h.foods.average[, 140:ncol(p.24h.foods.average)]
)

# Display the updated column names to confirm the change
names(p.24h.foods.average)

participants.recall.over2 <- read.csv(input_files$participants_recall_over2)

# load polyphenol and flavonoid checkout dataset
flavonoid.checkout <- read.csv(file.path(paths$data_raw, "flavonoid_ukbiobank.csv"))
new.row <- as.data.frame(matrix(0, ncol = ncol(flavonoid.checkout)))
colnames(new.row) <- colnames(flavonoid.checkout)
flavonoid.checkout <- rbind(new.row, flavonoid.checkout)
# Convert columns 9 to 223 in flavonoid.checkout to numeric
flavonoid.checkout[881, 95] <- 0
flavonoid.checkout[, 9:223] <- lapply(flavonoid.checkout[, 9:223], function(x) as.numeric(as.character(x)))
