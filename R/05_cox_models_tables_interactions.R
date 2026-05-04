# 05_cox_models_tables_interactions.R
# Auto-refactored from the original analysis notebook: github.Rmd.
# Run via targets::tar_make() from the project root.
# NOTE: UK Biobank raw data are not included in this repository.

source(file.path("R", "00_setup.R"))
load_packages()

# ---- Rmd chunk ----
# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
polyphenolwithhealth$Survival.Obj <- with(polyphenolwithhealth, Surv(survival.time, event))

# Loop through columns 3 to 89 (assumed dietary scores)
for (i in c(3:8)) {
  
  # Run univariate Cox model
  model.uni <- coxph(Survival.Obj ~ polyphenolwithhealth[,i] + Sex + Age.group + Ethnic.background + multimorbidity + IMD.quartile + Qualification, data = polyphenolwithhealth)
  
  # Extract variable name
  x <- names(polyphenolwithhealth)[i]
  y <- "Survival.Obj"
  
  # Extract model statistics
  coefficients <- coef(model.uni)[1]
  se <- summary(model.uni)$coefficients[1, "se(coef)"]
  z.value <- summary(model.uni)$coefficients[1, "z"]
  p.value <- summary(model.uni)$coefficients[1, "Pr(>|z|)"]
  ci.lower <- exp(confint(model.uni)[1, 1])
  ci.upper <- exp(confint(model.uni)[1, 2])
  hr <- exp(coefficients)
  
  # Format HR and CI
  formatted.hr <- sprintf("%.3f", hr)
  formatted.ci.lower <- sprintf("%.3f", ci.lower)
  formatted.ci.upper <- sprintf("%.3f", ci.upper)
  CI <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
  
  # Number of non-missing observations
  n <- sum(!is.na(polyphenolwithhealth[, i]))
  
  # Calculate person-years
  person.years <- sum(polyphenolwithhealth$survival.time[!is.na(polyphenolwithhealth[, i])], na.rm = TRUE)
  
  # Calculate number of cases and incidence
  cases <- sum(polyphenolwithhealth$group == "psoriasis" & !is.na(polyphenolwithhealth[, i]))
  incidence <- cases / n * 100
  cases.total <- paste0(cases, "/", n, " (", sprintf("%.3f", incidence), ")")

  # Store the results
  results <- data.frame(
    Nutrient = x,
    case = cases,
    total = n,
    incidence = incidence,
    Case.Total = cases.total,
    Person.Years = round(person.years, 1),
    HR.95.CI = CI,
    P.value = signif(p.value, 3)
  )
  
  # Append to the main results data frame
  dietary.score.results <- rbind(dietary.score.results, results)
}

# Print the results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid dose model 2.2.csv")

# ---- Rmd chunk ----
dietaryscorewithquartile <- polyphenolwithhealth %>%
  mutate(across(3:8, ~ cut(., breaks = quantile(., probs = 0:4 / 4, na.rm = TRUE), 
                             include.lowest = TRUE, labels = FALSE), 
                .names = "quartile_{col}"))

dietaryscorewithquartile <- dietaryscorewithquartile[, c(1:2, 72:77, 9:71)]


# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
dietaryscorewithquartile$Survival.Obj <- with(dietaryscorewithquartile, Surv(survival.time, event))

# Loop through dietary score columns (quartile variables)
for (i in c(72:79)) {
  
  # Run Cox model with quartiles treated as factors (Q1 as reference)
  model.uni <- coxph(Survival.Obj ~ factor(dietaryscorewithquartile[, i]) + Sex + Age.group + multimorbidity + Qualification, data = dietaryscorewithquartile)
  
  # Extract variable name
  nutrient <- names(dietaryscorewithquartile)[i]
  response <- "Survival.Obj"
  
  # Loop through Q2–Q4 (compared to Q1)
  for (j in 2:4) {
    
    # Extract statistics
    coef_j <- coef(model.uni)[j - 1]
    se_j <- summary(model.uni)$coefficients[j - 1, "se(coef)"]
    z_j <- summary(model.uni)$coefficients[j - 1, "z"]
    p_j <- summary(model.uni)$coefficients[j - 1, "Pr(>|z|)"]
    ci_lower <- exp(confint(model.uni)[j - 1, 1])
    ci_upper <- exp(confint(model.uni)[j - 1, 2])
    hr <- exp(coef_j)
    
    # Format HR and CI
    formatted.hr <- sprintf("%.2f", hr)
    formatted.ci.lower <- sprintf("%.2f", ci_lower)
    formatted.ci.upper <- sprintf("%.2f", ci_upper)
    hr_ci <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
    
    # Subset data for current quartile level
    quartile_data <- dietaryscorewithquartile %>%
      filter(as.integer(dietaryscorewithquartile[, i]) == j)
    
    # Case count and person-years
    n <- sum(!is.na(quartile_data[, i]))
    person.years <- sum(quartile_data$survival.time, na.rm = TRUE)
    cases <- sum(quartile_data$group == "psoriasis")
    incidence <- cases / n * 100
    cases.total <- paste0(cases, "/", n, " (", sprintf("%.3f", incidence), ")")
    
    # Store the result row
    result_row <- data.frame(
      Nutrient = nutrient,
      Quartile = paste0("Q", j),
      case = cases,
      total = n,
      incidence = incidence,
      Case.total = cases.total,
      person.years = person.years,
      HR.95.CI = hr_ci,
      P.value = signif(p_j, 3)
    )
    
    # Append to results
    dietary.score.results <- rbind(dietary.score.results, result_row)
  }
}

# Print results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid quartile model 2.2.csv", row.names = FALSE)

# ---- Rmd chunk ----
# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
dietaryscorewithquartile$Survival.Obj <- with(dietaryscorewithquartile, Surv(survival.time, event))

# Loop through columns 3 to 89 (assumed dietary scores)
for (i in c(72:79)) {
  
  # Run univariate Cox model
  model.uni <- coxph(Survival.Obj ~ as.numeric(dietaryscorewithquartile[, i]) + Sex + Age.group + multimorbidity + IMD.quartile + Qualification, data = dietaryscorewithquartile)
  
  # Extract variable name
  x <- names(dietaryscorewithquartile)[i]
  y <- "Survival.Obj"
  
  # Extract model statistics
  coefficients <- coef(model.uni)[1]
  se <- summary(model.uni)$coefficients[1, "se(coef)"]
  z.value <- summary(model.uni)$coefficients[1, "z"]
  p.value <- summary(model.uni)$coefficients[1, "Pr(>|z|)"]
  ci.lower <- exp(confint(model.uni)[1, 1])
  ci.upper <- exp(confint(model.uni)[1, 2])
  hr <- exp(coefficients)
  
  # Format HR and CI
  formatted.hr <- sprintf("%.3f", hr)
  formatted.ci.lower <- sprintf("%.3f", ci.lower)
  formatted.ci.upper <- sprintf("%.3f", ci.upper)
  CI <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
  
  # Number of non-missing observations
  n <- sum(!is.na(dietaryscorewithquartile[, i]))
  
  # Calculate person-years
  person.years <- sum(dietaryscorewithquartile$survival.time[!is.na(dietaryscorewithquartile[, i])], na.rm = TRUE)
  
  # Calculate number of cases and incidence
  cases <- sum(dietaryscorewithquartile$group == "psoriasis" & !is.na(dietaryscorewithquartile[, i]))
  incidence <- cases / n * 100
  cases.total <- paste0(cases, "/", n, " (", sprintf("%.3f", incidence), ")")

  # Store the results
  results <- data.frame(
    Nutrient = x,
    Case.Total = cases.total,
    case = cases,
    total = n,
    incidence = incidence,
    Person.Years = round(person.years, 1),
    HR.95.CI = CI,
    P.value = signif(p.value, 3)
  )
  
  # Append to the main results data frame
  dietary.score.results <- rbind(dietary.score.results, results)
}

# Print the results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid trend model 2.2.csv", row.names = FALSE)

# ---- Rmd chunk ----
# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
polyphenolwithhealth$Survival.Obj <- with(polyphenolwithhealth, Surv(survival.time, event))

# Loop through columns 3 to 89 (assumed dietary scores)
for (i in 3:8) {
  
  # Run univariate Cox model
  model.uni <- coxph(Survival.Obj ~ polyphenolwithhealth[,i] + Sex + Age.group + Ethnic.background + MET.Quintiles + Smoking.status...Instance.0 + alcohol.category + multimorbidity + IMD.quartile + Qualification, data = polyphenolwithhealth)
  
  # Extract variable name
  x <- names(polyphenolwithhealth)[i]
  y <- "Survival.Obj"
  
  # Extract model statistics
  coefficients <- coef(model.uni)[1]
  se <- summary(model.uni)$coefficients[1, "se(coef)"]
  z.value <- summary(model.uni)$coefficients[1, "z"]
  p.value <- summary(model.uni)$coefficients[1, "Pr(>|z|)"]
  ci.lower <- exp(confint(model.uni)[1, 1])
  ci.upper <- exp(confint(model.uni)[1, 2])
  hr <- exp(coefficients)
  
  # Format HR and CI
  formatted.hr <- sprintf("%.3f", hr)
  formatted.ci.lower <- sprintf("%.3f", ci.lower)
  formatted.ci.upper <- sprintf("%.3f", ci.upper)
  CI <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
  
  # Number of non-missing observations
  n <- sum(!is.na(polyphenolwithhealth[, i]))
  
  # Calculate person-years
  person.years <- sum(polyphenolwithhealth$survival.time[!is.na(polyphenolwithhealth[, i])], na.rm = TRUE)
  
  # Calculate number of cases and incidence
  cases <- sum(polyphenolwithhealth$group == "psoriasis" & !is.na(polyphenolwithhealth[, i]))
  incidence <- cases / n * 100
  cases.total <- paste0(cases, "/", n, " (", sprintf("%.3f", incidence), ")")

  # Store the results
  results <- data.frame(
    Nutrient = x,
    case = cases,
    total = n,
    incidence = incidence,
    Case.Total = cases.total,
    Person.Years = round(person.years, 1),
    HR.95.CI = CI,
    P.value = signif(p.value, 3)
  )
  
  # Append to the main results data frame
  dietary.score.results <- rbind(dietary.score.results, results)
}

# Print the results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid dose model 3.2.csv")

# ---- Rmd chunk ----
dietaryscorewithquartile <- polyphenolwithhealth %>%
  mutate(across(3:8, ~ cut(., breaks = quantile(., probs = 0:4 / 4, na.rm = TRUE), 
                             include.lowest = TRUE, labels = FALSE), 
                .names = "quartile_{col}"))

dietaryscorewithquartile <- dietaryscorewithquartile[, c(1:2, 72:77, 9:71)]


# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
dietaryscorewithquartile$Survival.Obj <- with(dietaryscorewithquartile, Surv(survival.time, event))

# Loop through dietary score columns (quartile variables)
for (i in c(72:79)) {
  
  # Run Cox model with quartiles treated as factors (Q1 as reference)
  model.uni <- coxph(Survival.Obj ~ factor(dietaryscorewithquartile[, i]) + Sex + Age.group + MET.Quintiles + Smoking.status...Instance.0 + alcohol.category + multimorbidity + IMD.quartile + Qualification, data = dietaryscorewithquartile)
  
  # Extract variable name
  nutrient <- names(dietaryscorewithquartile)[i]
  response <- "Survival.Obj"
  
  # Loop through Q2–Q4 (compared to Q1)
  for (j in 2:4) {
    
    # Extract statistics
    coef_j <- coef(model.uni)[j - 1]
    se_j <- summary(model.uni)$coefficients[j - 1, "se(coef)"]
    z_j <- summary(model.uni)$coefficients[j - 1, "z"]
    p_j <- summary(model.uni)$coefficients[j - 1, "Pr(>|z|)"]
    ci_lower <- exp(confint(model.uni)[j - 1, 1])
    ci_upper <- exp(confint(model.uni)[j - 1, 2])
    hr <- exp(coef_j)
    
    # Format HR and CI
    formatted.hr <- sprintf("%.2f", hr)
    formatted.ci.lower <- sprintf("%.2f", ci_lower)
    formatted.ci.upper <- sprintf("%.2f", ci_upper)
    hr_ci <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
    
    # Subset data for current quartile level
    quartile_data <- dietaryscorewithquartile %>%
      filter(as.integer(dietaryscorewithquartile[, i]) == j)
    
    # Case count and person-years
    n <- sum(!is.na(quartile_data[, i]))
    person.years <- sum(quartile_data$survival.time, na.rm = TRUE)
    cases <- sum(quartile_data$group == "psoriasis")
    incidence <- cases / n * 100
    cases.total <- paste0(cases, "/", n, " (", sprintf("%.2f", incidence), ")")
    
    # Store the result row
    result_row <- data.frame(
      Nutrient = nutrient,
      Quartile = paste0("Q", j),
      case = cases,
      total = n,
      incidence = incidence,
      Case.total = cases.total,
      person.years = person.years,
      HR.95.CI = hr_ci,
      P.value = signif(p_j, 3)
    )
    
    # Append to results
    dietary.score.results <- rbind(dietary.score.results, result_row)
  }
}

# Print results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid quartile model 3.2.csv", row.names = FALSE)

# ---- Rmd chunk ----
# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
dietaryscorewithquartile$Survival.Obj <- with(dietaryscorewithquartile, Surv(survival.time, event))

# Loop through columns 3 to 89 (assumed dietary scores)
for (i in c(72:79)) {
  
  # Run univariate Cox model
  model.uni <- coxph(Survival.Obj ~ as.numeric(dietaryscorewithquartile[, i]) + Sex + Age.group + MET.Quintiles + Smoking.status...Instance.0 + alcohol.category + multimorbidity + IMD.quartile + Qualification, data = dietaryscorewithquartile)
  
  # Extract variable name
  x <- names(dietaryscorewithquartile)[i]
  y <- "Survival.Obj"
  
  # Extract model statistics
  coefficients <- coef(model.uni)[1]
  se <- summary(model.uni)$coefficients[1, "se(coef)"]
  z.value <- summary(model.uni)$coefficients[1, "z"]
  p.value <- summary(model.uni)$coefficients[1, "Pr(>|z|)"]
  ci.lower <- exp(confint(model.uni)[1, 1])
  ci.upper <- exp(confint(model.uni)[1, 2])
  hr <- exp(coefficients)
  
  # Format HR and CI
  formatted.hr <- sprintf("%.3f", hr)
  formatted.ci.lower <- sprintf("%.3f", ci.lower)
  formatted.ci.upper <- sprintf("%.3f", ci.upper)
  CI <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
  
  # Number of non-missing observations
  n <- sum(!is.na(dietaryscorewithquartile[, i]))
  
  # Calculate person-years
  person.years <- sum(dietaryscorewithquartile$survival.time[!is.na(dietaryscorewithquartile[, i])], na.rm = TRUE)
  
  # Calculate number of cases and incidence
  cases <- sum(dietaryscorewithquartile$group == "psoriasis" & !is.na(dietaryscorewithquartile[, i]))
  incidence <- cases / n * 100
  cases.total <- paste0(cases, "/", n, " (", sprintf("%.3f", incidence), ")")

  # Store the results
  results <- data.frame(
    Nutrient = x,
    Case.Total = cases.total,
    case = cases,
    total = n,
    incidence = incidence,
    Person.Years = round(person.years, 1),
    HR.95.CI = CI,
    P.value = signif(p.value, 3)
  )
  
  # Append to the main results data frame
  dietary.score.results <- rbind(dietary.score.results, results)
}

# Print the results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid trend model 3.2.csv", row.names = FALSE)

# ---- Rmd chunk ----
# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
polyphenolwithhealth$Survival.Obj <- with(polyphenolwithhealth, Surv(survival.time, event))

# Loop through columns 3 to 89 (assumed dietary scores)
for (i in c(3:8)) {

  # Run sunivariate Cox model
  model.uni <- coxph(Survival.Obj ~ polyphenolwithhealth[,i] + Sex + Age.group + MET.Quintiles + Smoking.status...Instance.0 + alcohol.category + multimorbidity + IMD.quartile + Qualification + Energy + Free.sugar + meat + wholegrain + BMI.Category, data = polyphenolwithhealth)
  
  # Extract variable name
  x <- names(polyphenolwithhealth)[i]
  y <- "Survival.Obj"
  
  # Extract model statistics
  coefficients <- coef(model.uni)[1]
  se <- summary(model.uni)$coefficients[1, "se(coef)"]
  z.value <- summary(model.uni)$coefficients[1, "z"]
  p.value <- summary(model.uni)$coefficients[1, "Pr(>|z|)"]
  ci.lower <- exp(confint(model.uni)[1, 1])
  ci.upper <- exp(confint(model.uni)[1, 2])
  hr <- exp(coefficients)
  
  # Format HR and CI
  formatted.hr <- sprintf("%.3f", hr)
  formatted.ci.lower <- sprintf("%.3f", ci.lower)
  formatted.ci.upper <- sprintf("%.3f", ci.upper)
  CI <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
  
  # Number of non-missing observations
  n <- sum(!is.na(polyphenolwithhealth[, i]))
  
  # Calculate person-years
  person.years <- sum(polyphenolwithhealth$survival.time[!is.na(polyphenolwithhealth[, i])], na.rm = TRUE)
  
  # Calculate number of cases and incidence
  cases <- sum(polyphenolwithhealth$group == "psoriasis" & !is.na(polyphenolwithhealth[, i]))
  incidence <- cases / n * 100
  cases.total <- paste0(cases, "/", n, " (", sprintf("%.3f", incidence), ")")

  # Store the results
  results <- data.frame(
    Nutrient = x,
    case = cases,
    total = n,
    incidence = incidence,
    Case.Total = cases.total,
    Person.Years = round(person.years, 1),
    HR.95.CI = CI,
    P.value = signif(p.value, 3)
  )
  
  # Append to the main results data frame
  dietary.score.results <- rbind(dietary.score.results, results)
}

# Print the results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid dose model PRS.csv")

# ---- Rmd chunk ----


dietaryscorewithquartile <- polyphenolwithhealth %>%
  mutate(across(3:8, ~ cut(., breaks = quantile(., probs = 0:4 / 4, na.rm = TRUE), 
                             include.lowest = TRUE, labels = FALSE), 
                .names = "quartile_{col}"))

dietaryscorewithquartile <- dietaryscorewithquartile[, c(1:2, 72:77, 9:71)]

# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
dietaryscorewithquartile$Survival.Obj <- with(dietaryscorewithquartile, Surv(survival.time, event))

# Loop through dietary score columns (quartile variables)
for (i in c(72:79)) {
  
  # Run Cox model with quartiles treated as factors (Q1 as reference)
  model.uni <- coxph(Survival.Obj ~ factor(dietaryscorewithquartile[, i]) + Sex + Age.group + MET.Quintiles + Smoking.status...Instance.0 + alcohol.category + multimorbidity + IMD.quartile + Qualification + Energy + Free.sugar + meat + wholegrain +BMI.Category, data = dietaryscorewithquartile)
  
  # Extract variable name
  nutrient <- names(dietaryscorewithquartile)[i]
  response <- "Survival.Obj"
  
  # Loop through Q2–Q4 (compared to Q1)
  for (j in 2:4) {
    
    # Extract statistics
    coef_j <- coef(model.uni)[j - 1]
    se_j <- summary(model.uni)$coefficients[j - 1, "se(coef)"]
    z_j <- summary(model.uni)$coefficients[j - 1, "z"]
    p_j <- summary(model.uni)$coefficients[j - 1, "Pr(>|z|)"]
    ci_lower <- exp(confint(model.uni)[j - 1, 1])
    ci_upper <- exp(confint(model.uni)[j - 1, 2])
    hr <- exp(coef_j)
    
    # Format HR and CI
    formatted.hr <- sprintf("%.2f", hr)
    formatted.ci.lower <- sprintf("%.2f", ci_lower)
    formatted.ci.upper <- sprintf("%.2f", ci_upper)
    hr_ci <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
    
    # Subset data for current quartile level
    quartile_data <- dietaryscorewithquartile %>%
      filter(as.integer(dietaryscorewithquartile[, i]) == j)
    
    # Case count and person-years
    n <- sum(!is.na(quartile_data[, i]))
    person.years <- sum(quartile_data$survival.time, na.rm = TRUE)
    cases <- sum(quartile_data$group == "psoriasis")
    incidence <- cases / n * 100
    cases.total <- paste0(cases, "/", n, " (", sprintf("%.2f", incidence), ")")
    
    # Store the result row
    result_row <- data.frame(
      Nutrient = nutrient,
      Quartile = paste0("Q", j),
      case = cases,
      total = n,
      incidence = incidence,
      Case.total = cases.total,
      person.years = person.years,
      HR.95.CI = hr_ci,
      P.value = signif(p_j, 3)
    )
    
    # Append to results
    dietary.score.results <- rbind(dietary.score.results, result_row)
  }
}

# Print results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid quartile model PRS.csv", row.names = FALSE)

# ---- Rmd chunk ----
# Initialize an empty data frame to store the results
dietary.score.results <- data.frame()

# Define the survival object
dietaryscorewithquartile$Survival.Obj <- with(dietaryscorewithquartile, Surv(survival.time, event))

# Loop through columns 3 to 89 (assumed dietary scores)
for (i in c(72:79)) {
  
  # Run univariate Cox model
  model.uni <- coxph(Survival.Obj ~ as.numeric(dietaryscorewithquartile[, i]) + Sex + Age.group + MET.Quintiles + Smoking.status...Instance.0 + alcohol.category + multimorbidity + IMD.quartile + Qualification + Energy + Free.sugar + meat + wholegrain + BMI.Category, data = dietaryscorewithquartile)
  
  # Extract variable name
  x <- names(dietaryscorewithquartile)[i]
  y <- "Survival.Obj"
  
  # Extract model statistics
  coefficients <- coef(model.uni)[1]
  se <- summary(model.uni)$coefficients[1, "se(coef)"]
  z.value <- summary(model.uni)$coefficients[1, "z"]
  p.value <- summary(model.uni)$coefficients[1, "Pr(>|z|)"]
  ci.lower <- exp(confint(model.uni)[1, 1])
  ci.upper <- exp(confint(model.uni)[1, 2])
  hr <- exp(coefficients)
  
  # Format HR and CI
  formatted.hr <- sprintf("%.3f", hr)
  formatted.ci.lower <- sprintf("%.3f", ci.lower)
  formatted.ci.upper <- sprintf("%.3f", ci.upper)
  CI <- paste0(formatted.hr, " (", formatted.ci.lower, "-", formatted.ci.upper, ")")
  
  # Number of non-missing observations
  n <- sum(!is.na(dietaryscorewithquartile[, i]))
  
  # Calculate person-years
  person.years <- sum(dietaryscorewithquartile$survival.time[!is.na(dietaryscorewithquartile[, i])], na.rm = TRUE)
  
  # Calculate number of cases and incidence
  cases <- sum(dietaryscorewithquartile$group == "psoriasis" & !is.na(dietaryscorewithquartile[, i]))
  incidence <- cases / n * 100
  cases.total <- paste0(cases, "/", n, " (", sprintf("%.3f", incidence), ")")

  # Store the results
  results <- data.frame(
    Nutrient = x,
    Case.Total = cases.total,
    case = cases,
    total = n,
    incidence = incidence,
    Person.Years = round(person.years, 1),
    HR.95.CI = CI,
    P.value = signif(p.value, 3)
  )
  
  # Append to the main results data frame
  dietary.score.results <- rbind(dietary.score.results, results)
}

# Print the results
print(dietary.score.results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(dietary.score.results, "diet quality output/flavonoid trend model PRS.csv", row.names = FALSE)

# ---- Rmd chunk ----
############## Show in table format ###################

# Define variables to be displayed as Median (IQR)
nonnormal_vars <- c(
  "Free.sugar", "meat", "wholegrain",
  "anthocyanins", "flavanones", "flavonols",
  "flavones", "flavan.3.ols", "total.flavonoid"
)

# Create the TableOne object
description.DQI <- CreateTableOne(
  data = polyphenolwithhealth,
  vars = c(
    "Ethnic.background", "Sex", "Age.at.recruitment",
    "Summed.MET.minutes.per.week.for.all.activity...Instance.0",
    "Smoking.status...Instance.0", "alcohol.category",
    "multimorbidity", "BMI", "PRS.tertile",
    "Free.sugar", "meat", "wholegrain",
    "anthocyanins", "flavanones", "flavonols",
    "flavones", "flavan.3.ols", "total.flavonoid",
    "Townsend.deprivation.index.at.recruitment",
    "Qualification", "Energy"
  ),
  
  factorVars = c(
    "Ethnic.background", "Sex", "Smoking.status...Instance.0",
    "Alcohol.drinker.status...Instance.0", "multimorbidity",
    "Qualification", "PRS.tertile"
  ),
  
  strata = "group"
)

# Print Table 1 with Median (IQR) and automatically generated P values
table1_output <- print(
  description.DQI,
  quote = FALSE,
  noSpaces = TRUE,
  nonnormal = nonnormal_vars,
  test = TRUE
)

# Display the output
table1_output

# Export the table as a CSV file
# Publication note: readr::write_csv is preferred for clean and reproducible CSV output
write.csv(table1_output, "Table1.csv")


# ---- Rmd chunk ----
# Exclude participants with missing PRS tertile information
dietaryscorewithquartilek <- dietaryscorewithquartile %>%
  filter(PRS.tertile != "Missing")

# Drop unused factor levels and set the desired PRS tertile order
dietaryscorewithquartilek$PRS.tertile <- factor(
  droplevels(dietaryscorewithquartilek$PRS.tertile),
  levels = c("Low", "Medium", "High")
)

# Fit the Cox model including the interaction between flavone quartiles and PRS tertile
model.1 <- coxph(
  Survival.Obj ~ factor(flavones_ord4) * PRS.tertile +
    Sex + Age.group + MET.Quintiles +
    Smoking.status...Instance.0 + alcohol.category + multimorbidity +
    IMD.quartile + Qualification + Energy + Free.sugar + meat +
    wholegrain + BMI.Category + Genotype.measurement.batch +
    Genetic.principal.components...Array.1 +
    Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 +
    Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 +
    Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 +
    Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 +
    Genetic.principal.components...Array.10,
  data = dietaryscorewithquartile
)

# Fit the Cox model without the interaction term
model.2 <- coxph(
  Survival.Obj ~ factor(flavones_ord4) + PRS.tertile +
    Sex + Age.group + MET.Quintiles +
    Smoking.status...Instance.0 + alcohol.category + multimorbidity +
    IMD.quartile + Qualification + Energy + Free.sugar + meat +
    wholegrain + BMI.Category + Genotype.measurement.batch +
    Genetic.principal.components...Array.1 +
    Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 +
    Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 +
    Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 +
    Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 +
    Genetic.principal.components...Array.10,
  data = dietaryscorewithquartile
)

# Compare models using a likelihood ratio test
anova(model.2, model.1, test = "LRT")

# Print the final interaction results table
print(interaction_results)

# Interaction P value: 0.5188


# ---- Rmd chunk ----
# Create the survival object
polyphenolwithhealth$Survival.Obj <- with(
  polyphenolwithhealth,
  Surv(survival.time, event)
)

# Fit the Cox model
model <- coxph(
  Survival.Obj ~ total.flavonoid * PRS.tertile +
    Sex + Age.group + Ethnic.background + MET.Quintiles +
    Smoking.status...Instance.0 + alcohol.category + multimorbidity +
    IMD.quartile + Energy + Free.sugar + meat + Sodium +
    AQI.Category + BMI.Category + Genotype.measurement.batch +
    Genetic.principal.components...Array.1 +
    Genetic.principal.components...Array.2 +
    Genetic.principal.components...Array.3 +
    Genetic.principal.components...Array.4 +
    Genetic.principal.components...Array.5 +
    Genetic.principal.components...Array.6 +
    Genetic.principal.components...Array.7 +
    Genetic.principal.components...Array.8 +
    Genetic.principal.components...Array.9 +
    Genetic.principal.components...Array.10,
  data = polyphenolwithhealth
)

# Extract HRs, confidence intervals, and P values for selected terms
term_ids <- c(1, 158, 159)

coef_values <- coef(model)[term_ids]
ci_values <- confint(model)[term_ids, ]

hr_values <- exp(coef_values)
ci_exp <- exp(ci_values)

p_values <- summary(model)$coefficients[term_ids, "Pr(>|z|)"]

# Format the results
results <- data.frame(
  Term = c(
    "total.flavonoid (Low PRS)",
    "total.flavonoid:PRS.tertileMedium",
    "total.flavonoid:PRS.tertileHigh"
  ),
  HR = sprintf("%.3f", hr_values),
  CI = paste0(
    "(",
    sprintf("%.3f", ci_exp[, 1]),
    "-",
    sprintf("%.3f", ci_exp[, 2]),
    ")"
  ),
  P.value = signif(p_values, 3)
)

# Combine HR and 95% CI into one column
results$HR.95.CI <- paste0(results$HR, " ", results$CI)

results <- results[, c("Term", "HR.95.CI", "P.value")]

# Print the results
print(results)

# Save the results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(results, "coxph_interaction_with_incidence total.flavonoid.csv")
