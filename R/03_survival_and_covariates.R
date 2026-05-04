# 03_survival_and_covariates.R
# Auto-refactored from the original analysis notebook: github.Rmd.
# Run via targets::tar_make() from the project root.
# NOTE: UK Biobank raw data are not included in this repository.

source(file.path("R", "00_setup.R"))
load_packages()

# ---- Rmd chunk ----
endtime <- read.csv("C:/Users/14999/OneDrive - King's College London/genetics/check the name of different sets.csv")
# Step 1: Convert ISO datetime strings to Date format (keep only the date)
endtime <- endtime %>%
  mutate(across(starts_with("When.diet.questionnaire.completed"),
                ~ as.Date(.x, format = "%Y-%m-%dT%H:%M:%OSZ")))

# Step 2: Calculate the latest (most recent) non-NA date across the 5 time points
endtime$Latest_Diet_Date <- apply(endtime[, c("When.diet.questionnaire.completed...Instance.0",
                                              "When.diet.questionnaire.completed...Instance.1",
                                              "When.diet.questionnaire.completed...Instance.2",
                                              "When.diet.questionnaire.completed...Instance.3",
                                              "When.diet.questionnaire.completed...Instance.4")],
                                   1, function(row) {
                                     if (all(is.na(row) | row == "")) {
                                       return(NA)
                                     } else {
                                       return(max(row, na.rm = TRUE))
                                     }
                                   })
endtime$Date.of.death...Instance.0 <- as.Date(endtime$Date.of.death...Instance.0, format = "%d/%m/%Y")

polyphenol <- read.csv("C:/Users/14999/OneDrive - King's College London/genetics/flavonoid/selected.polyphenol.csv")
#polyphenol <- polyphenol[, c(1,2,3,6,7,24,42,65,76,90,104,111,158,213)]

# Convert the 'group' column to a factor with specified levels
polyphenol$group <- factor(polyphenol$group, levels = c("non-psoriasis", "psoriasis"))

# Add a new column 'event' indicating the presence of psoriasis (1) or not (0)
polyphenol <- polyphenol %>%
  mutate(event = ifelse(group == "psoriasis", 1, 0))

# Load data related to psoriasis incidence and format date columns
# Load data related to psoriasis incidence and format date columns
filtered.incidence <- read.csv(input_files$incidence)
filtered.incidence <- merge(filtered.incidence, endtime, by = "Participant.ID")
filtered.incidence$Date.L40.first.reported..psoriasis. <- as.Date(filtered.incidence$Date.L40.first.reported..psoriasis., format = "%Y-%m-%d")

# Calculate survival time from the date of attending assessment centre to the first reported psoriasis date
filtered.incidence <- filtered.incidence %>%
  mutate(
    end_time = pmin(Date.L40.first.reported..psoriasis., Date.of.death...Instance.0, na.rm = TRUE),
    survival.time = interval(Latest_Diet_Date, end_time) / years(1)
  ) %>%
  filter(survival.time >  0)

# Load control group data and merge with non-psoriasis baseline characteristics
filtered.control <- read.csv(input_files$control)
filtered.control <- merge(filtered.control, endtime, by = "Participant.ID")


# Calculate survival time for control participants from the date of attending assessment centre
filtered.control <- filtered.control %>%
  mutate(
    end_time = pmin(as.Date("2023-10-12"), Date.of.death...Instance.0, na.rm = TRUE),
    survival.time = interval(Latest_Diet_Date, end_time) / years(1)
  ) %>%
  filter(survival.time > 0)

# Combine survival times from both incidence and control groups into one dataset
survivaltime <- rbind(filtered.incidence[,c(1,23)], filtered.control[,c(1,22)])

# Merge dietary scores dataset with survival times based on 'Participant.ID'
polyphenolwithst <- merge(polyphenol, survivaltime, by = "Participant.ID")

# Display the column names of the final dataset
names(polyphenolwithst)

#polyphenolwithst <- polyphenolwithst %>% filter(Participant.ID %in% participants.recall.over2$Participant.ID)

table(polyphenolwithst$group)
participants.recall.over2 <- read.csv(input_files$participants_recall_over2)

#polyphenolwithst <- polyphenolwithst %>% filter(Participant.ID %in% participants.recall.over2$Participant.ID)


# ---- Rmd chunk ----
# Load demographic data for participants with psoriasis
p.baseline <- read.csv(input_files$demographics) %>%
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)
p.baseline <- p.baseline[, c(1,2,3, 6, 16)]

white <- c("White", "British", "Irish", "Any other white background")
mixed <- c("White and Black Caribbean", "White and Black African", "White and Asian", "Mixed", "Any other mixed background")
asian <- c("Indian", "Chinese", "Pakistani", "Bangladeshi", "Asian or Asian British", "Any other Asian background")
black <- c("Caribbean", "African", "Black or Black British", "Any other Black background")
other.ethnic.group <- "Other ethnic group"
missing <- c("Do not know", "Prefer not to answer", "")

p.baseline <- p.baseline %>%
  mutate(
    Ethnic.background = case_when(
      p.baseline$Ethnic.background...Instance.0 %in% white ~ "White",
      p.baseline$Ethnic.background...Instance.0 %in% mixed ~ "Mixed",
      p.baseline$Ethnic.background...Instance.0 %in% asian ~ "Asian",
      p.baseline$Ethnic.background...Instance.0 %in% black ~ "Black",
      p.baseline$Ethnic.background...Instance.0 %in% other.ethnic.group ~ "Other ethnic group",
      p.baseline$Ethnic.background...Instance.0 %in% missing ~ "Missing"
  ))

education <- read.csv("basic files/education.csv") %>%
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)

qual_order <- c(
  "College or University degree",
  "Other professional qualifications eg: nursing, teaching",
  "A levels/AS levels or equivalent",
  "NVQ or HND or HNC or equivalent",
  "CSEs or equivalent",
  "O levels/GCSEs or equivalent",
  "Prefer not to answer",
  "None of the above",
  NA
)

# 处理数据：按优先级保留每个ID的最高学历
education_top <- education %>%
  mutate(priority = match(Qualifications, qual_order)) %>%
  group_by(Participant.ID) %>%
  slice_min(priority, with_ties = FALSE) %>%
  ungroup()


#Low: CSEs or equivalent, O levels/GCSEs or equivalent; 
low <- c("CSEs or equivalent", "O levels/GCSEs or equivalent")
#Medium: A levels/AS levels or equivalent, NVQ or HND or HNC or equivalent; 
medium <- c("A levels/AS levels or equivalent", "NVQ or HND or HNC or equivalent")
#High: College or University degree, other professional qualifications eg: nursing, teaching; 
high <- c("College or University degree", "Other professional qualifications eg: nursing, teaching")
#unknown/missing/prefer not to say
missing <- c(NA, "None of the above", "Prefer not to answer")
education_top <- education_top %>%
  mutate(
    Qualification = case_when(
      education_top$Qualifications %in% low ~ "low",
      education_top$Qualifications %in% medium ~ "medium",
      education_top$Qualifications %in% high ~ "high",
      education_top$Qualifications %in% missing ~ "missing"
  ))
p.baseline <- merge(p.baseline, education_top[, c(1,5)], by = "Participant.ID")

# Load lifestyle data for participants with and without psoriasis
p.lifestyle <- read.csv("basic files/healthandlifestyle.csv") %>%
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)
p.lifestyle <- p.lifestyle[, c(1, 12, 71:72, 85, 87, 88)]

# Rename the BMI column to a simpler name
names(p.lifestyle)[names(p.lifestyle) == "Body.mass.index..BMI....Instance.0.participant...p21001_i0."] <- "BMI"

covariatebasic <- merge(p.baseline, p.lifestyle, by = "Participant.ID")

p.comorbid <- read.csv("basic files/comorbidity.number.csv")%>%
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)
AQI <- read.csv("basic files/AQI.csv")%>%
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)
PRS <- read.csv("basic files/PRS.csv") %>%
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)

need <- merge(p.comorbid[, c(1, 155, 156)], AQI[, c(1,9)], by = "Participant.ID")
covariate <- merge(covariatebasic, need, by = "Participant.ID")
covariateothers <- merge(endtime[, c(3,9)], PRS[, c(1:15)], by = "Participant.ID")
covariate <- merge(covariate, covariateothers, by = "Participant.ID")

p.24h.foods.average <- read.csv(input_files$average_24h_foods) %>% 
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)
p.24h.nutrients.average <- read.csv("24hr/average 24hrecall nutrients.csv") %>% 
  filter(Participant.ID %in% polyphenolwithst$Participant.ID)
PDIscore <- read.csv("diet indices/PDIscore.csv")
nutrient <- merge(p.24h.nutrients.average[, c(1, 2, 14, 17, 20, 44)], PDIscore[,c(21, 18, 19)], by = "Participant.ID")

p.24h.foods.average$SSB <- rowSums(p.24h.foods.average[c(4, 5, 22, 24)], na.rm = TRUE)
p.24h.foods.average$cheese <- rowSums(p.24h.foods.average[c(106:116)], na.rm = TRUE)
p.24h.foods.average$egg <- rowSums(p.24h.foods.average[c(117:121)], na.rm = TRUE)
p.24h.foods.average$meat <- rowSums(p.24h.foods.average[c(122:125, 128:131)], na.rm = TRUE)
p.24h.foods.average$fish <- rowSums(p.24h.foods.average[c(132:140)], na.rm = TRUE)
p.24h.foods.average$wholegrain <- rowSums(p.24h.foods.average[c(33:35, 38:39, 48, 99, 101, 104)], na.rm = TRUE)
p.24h.foods.average$coffee <- rowSums(p.24h.foods.average[c(11:15)], na.rm = TRUE)
p.24h.foods.average$tea  <- rowSums(p.24h.foods.average[c(16:20)], na.rm = TRUE)
p.24h.foods.average$fruit.juice <- rowSums(p.24h.foods.average[c(6:9)], na.rm = TRUE)
p.24h.foods.average$fruit  <- rowSums(p.24h.foods.average[c(183:201)], na.rm = TRUE)
p.24h.foods.average$chocolate  <- rowSums(p.24h.foods.average[, c("Chocolate.bar.intake",
                                                 "White.chocolate.intake",
                                                 "Milk.chocolate.intake",
                                                 "Dark.chocolate.intake",
                                                 "Chocolate.covered.raisin.intake",
                                                 "Chocolate.sweet.intake",
                                                 "Chocolate.covered.biscuits.intake",
                                                 "Chocolate.biscuits.intake")], na.rm = TRUE)
p.24h.foods.average$vegetable  <- rowSums(p.24h.foods.average[c(150:154, 156:171, 173:182)], na.rm = TRUE)
p.24h.foods.average$wine  <- rowSums(p.24h.foods.average[c(26:31)], na.rm = TRUE)
p.24h.foods.average$nut  <- rowSums(p.24h.foods.average[c(85:89)], na.rm = TRUE)
p.24h.foods.average$legume  <- rowSums(p.24h.foods.average[c(145:146, 155, 172)], na.rm = TRUE)

#p.24h.foods.average$snk <- rowSums(p.24h.foods.average[c(50:71, 78:79, 82:84, 90:92, 94)], na.rm = TRUE)

covariate <- merge(covariate, nutrient, by = "Participant.ID")
covariate <- merge(covariate, p.24h.foods.average[, c(1, 202:216)], by = "Participant.ID")
polyphenolwithhealth <- merge(polyphenolwithst, covariate, by = "Participant.ID")

# Convert IPAQ activity group to a factor with specified levels, excluding blank entries
polyphenolwithhealth$IPAQ.activity.group...Instance.0 <- factor(polyphenolwithhealth$IPAQ.activity.group...Instance.0, levels = c("low", "moderate", "high", ""))

# Step 1: Replace "" and "Prefer not to answer" with "Missing"
polyphenolwithhealth$Smoking.status...Instance.0[
  polyphenolwithhealth$Smoking.status...Instance.0 %in% c("", "Prefer not to answer")
] <- "Missing"

# Step 2: Set factor levels, including "Missing"
polyphenolwithhealth$Smoking.status...Instance.0 <- factor(
  polyphenolwithhealth$Smoking.status...Instance.0,
  levels = c("Never", "Previous", "Current", "Missing")
)

# Step 1: Replace "" and "Prefer not to answer" with "Missing"
polyphenolwithhealth$Alcohol.drinker.status...Instance.0[
  polyphenolwithhealth$Alcohol.drinker.status...Instance.0 %in% c("", "Prefer not to answer")
] <- "Missing"

# Step 2: Set factor levels, including "Missing"
polyphenolwithhealth$Alcohol.drinker.status...Instance.0 <- factor(
  polyphenolwithhealth$Alcohol.drinker.status...Instance.0,
  levels = c("Never", "Previous", "Current", "Missing")
)

# Step 1: Replace "" and "Prefer not to answer" with "Missing"
polyphenolwithhealth$Alcohol.intake.frequency....Instance.0[
  polyphenolwithhealth$Alcohol.intake.frequency....Instance.0 %in% c("", "Prefer not to answer")
] <- "Missing"

# Step 2: Set factor levels, including "Missing"
polyphenolwithhealth$Alcohol.drinker.status...Instance.0 <- factor(polyphenolwithhealth$Alcohol.intake.frequency....Instance.0, levels = c("Never", "Special occasions only", "One to three times a month", "Once or twice a week", "Three or four times a week", "Daily or almost daily", "Missing"))


# Convert MET minutes
polyphenolwithhealth$Summed.MET.minutes.per.week.for.all.activity...Instance.0 <- 
  as.numeric(polyphenolwithhealth$Summed.MET.minutes.per.week.for.all.activity...Instance.0)

# Calculate unique quantile cut-points
valid_MET <- polyphenolwithhealth$Summed.MET.minutes.per.week.for.all.activity...Instance.0
quantile_breaks <- unique(quantile(valid_MET, probs = seq(0, 1, 0.2), na.rm = TRUE))

# Use cut() to create quintile groups for non-missing values
MET_quintiles <- cut(
  valid_MET,
  breaks = quantile_breaks,
  include.lowest = TRUE,
  labels = paste0("Q", 1:5)
)

# Replace missing values with "Missing"
MET_quintiles <- as.character(MET_quintiles)
MET_quintiles[is.na(MET_quintiles)] <- "Missing"

polyphenolwithhealth$MET.Quintiles <- factor(
  MET_quintiles,
  levels = c(paste0("Q", 1:5), "Missing")
)


# Step 1: Create BMI category for non-missing values
polyphenolwithhealth$BMI.Category <- cut(
  polyphenolwithhealth$BMI,
  breaks = c(-Inf, 24.9, 29.9, Inf),
  labels = c("<24.9", "25.0–29.9", "≥30"),
  include.lowest = TRUE
)

# Step 2: Convert to character to allow adding "Missing"
polyphenolwithhealth$BMI.Category <- as.character(polyphenolwithhealth$BMI.Category)

# Step 3: Replace missing values with "Missing"
polyphenolwithhealth$BMI.Category[is.na(polyphenolwithhealth$BMI.Category)] <- "Missing"

# Step 4: Convert back to factor and set the desired level order
polyphenolwithhealth$BMI.Category <- factor(
  polyphenolwithhealth$BMI.Category,
  levels = c("<24.9", "25.0–29.9", "≥30", "Missing")
)


# Create PRS tertiles
prs_cut <- cut(
  polyphenolwithhealth$Standard.PRS.for.psoriasis..PSO.,
  breaks = quantile(
    polyphenolwithhealth$Standard.PRS.for.psoriasis..PSO.,
    probs = c(0, 1/3, 2/3, 1),
    na.rm = TRUE
  ),
  include.lowest = TRUE,
  labels = c("Low", "Medium", "High")
)

prs_cut <- as.character(prs_cut)
prs_cut[is.na(prs_cut)] <- "Missing"

polyphenolwithhealth$PRS.tertile <- factor(
  prs_cut,
  levels = c("Low", "Medium", "High", "Missing")
)


# Create age groups
age_cut <- cut(
  polyphenolwithhealth$Age.at.recruitment,
  breaks = c(40, 50, 60, 70),
  include.lowest = TRUE,
  right = FALSE,
  labels = c("40–49", "50–59", "60–69")
)

age_cut <- as.character(age_cut)

polyphenolwithhealth$Age.group <- factor(
  age_cut,
  levels = c("40–49", "50–59", "60–69")
)


# Create IMD quartiles
imd_cut <- cut(
  as.numeric(polyphenolwithhealth$Townsend.deprivation.index.at.recruitment),
  breaks = quantile(
    polyphenolwithhealth$Townsend.deprivation.index.at.recruitment,
    probs = seq(0, 1, 0.25),
    na.rm = TRUE
  ),
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

imd_cut <- as.character(imd_cut)
imd_cut[is.na(imd_cut)] <- "Missing"

polyphenolwithhealth$IMD.quartile <- factor(
  imd_cut,
  levels = c("Q1", "Q2", "Q3", "Q4", "Missing")
)


# Convert AQI to a factor with the desired level order
polyphenolwithhealth$AQI.Category <- factor(
  polyphenolwithhealth$AQI,
  levels = c("1", "2", "Missing")
)


# Categorise alcohol intake and convert to an ordered factor
polyphenolwithhealth <- polyphenolwithhealth %>%
  mutate(
    alcohol.category = case_when(
      Alcohol < 1 ~ "<1 g/day",
      Alcohol >= 1 & Alcohol < 8 ~ "1–7 g/day",
      Alcohol >= 8 & Alcohol < 16 ~ "8–15 g/day",
      Alcohol >= 16 ~ "16+ g/day"
    ),
    alcohol.category = factor(
      alcohol.category,
      levels = c("<1 g/day", "1–7 g/day", "8–15 g/day", "16+ g/day")
    )
  )


# Create questionnaire completion level and multimorbidity categories
polyphenolwithhealth <- polyphenolwithhealth %>%
  mutate(
    level.questionnaires = factor(
      if_else(
        Number.of.diet.questionnaires.completed == 1,
        "level1",
        "level2"
      ),
      levels = c("level1", "level2")
    ),
    
    multimorbidity = factor(
      case_when(
        multimorbidity == 0 ~ "level0",
        multimorbidity == 1 ~ "level1",
        multimorbidity >= 2 ~ "level2"
      ),
      levels = c("level0", "level1", "level2")
    )
  )


# Replace empty strings with missing values in selected columns
polyphenolwithhealth[, 11:70][polyphenolwithhealth[, 11:70] == ""] <- NA


# Restrict the dataset to participants with more than two dietary recalls
polyphenolwithhealth <- polyphenolwithhealth %>%
  filter(Participant.ID %in% participants.recall.over2$Participant.ID)


# ---- Rmd chunk ----
df <- polyphenolwithhealth


# Define cut-points for each flavonoid variable in mg/day
specs4 <- list(
  anthocyanins = list(
    breaks = c(-Inf, 0.5, 25, 70, Inf),
    labels = c("<0.5", "0.5–<25", "25–<70", "≥70")
  ),
  
  `flavan.3.ols` = list(
    breaks = c(-Inf, 200, 500, 900, Inf),
    labels = c("<200", "200–<500", "500–<900", "≥900")
  ),
  
  flavanones = list(
    breaks = c(-Inf, 2, 20, 45, Inf),
    labels = c("<2", "2–<20", "20–<45", "≥45")
  ),
  
  flavonols = list(
    breaks = c(-Inf, 25, 40, 55, Inf),
    labels = c("<25", "25–<40", "40–<55", "≥55")
  ),
  
  flavones = list(
    breaks = c(-Inf, 0.5, 1.5, 3, Inf),
    labels = c("<0.5", "0.5–<1.5", "1.5–<3", "≥3")
  ),
  
  total.flavonoid = list(
    breaks = c(-Inf, 400, 800, 1100, Inf),
    labels = c("<450", "450–<800", "800–<1100", "≥1100")
  )
)


# Batch-create categorical variables (_cat4) and ordinal variables (_ord4)
for (nm in names(specs4)) {
  br <- specs4[[nm]]$breaks
  lab <- specs4[[nm]]$labels
  
  cat_var <- paste0(nm, "_cat4")
  ord_var <- paste0(nm, "_ord4")
  
  # Create categorical variable
  df[[cat_var]] <- cut(
    df[[nm]],
    breaks = br,
    right = FALSE,
    labels = lab
  )
  
  # Create ordinal variable
  df[[ord_var]] <- as.numeric(df[[cat_var]])
}


# Identify all newly created categorical and ordinal variables
cat_vars <- grep("_cat4$", names(df), value = TRUE)
ord_vars <- grep("_ord4$", names(df), value = TRUE)


# View all categorical variables
cat_table <- df %>%
  dplyr::select(all_of(cat_vars))


# View all ordinal variables
ord_table <- df %>%
  dplyr::select(all_of(ord_vars))


# Print the first few rows of categorical labels
cat_table %>% head()


# Print the first few rows of ordinal variables coded from 1 to 4
ord_table %>% head()


# Save the updated dataset
dietaryscorewithquartile <- df
