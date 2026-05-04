# 04_food_source_and_spline_analysis.R
# Auto-refactored from the original analysis notebook: github.Rmd.
# Run via targets::tar_make() from the project root.
# NOTE: UK Biobank raw data are not included in this repository.

source(file.path("R", "00_setup.R"))
load_packages()

# ---- Rmd chunk ----
df <- p.24h.foods.average %>%
  # Step 1: Extract the original variables
  dplyr::select(
    Participant.ID,
    Pure.fruit.vegetable.juice.intake,
    Olives.intake,
    Broccoli.intake,
    Celery.intake,
    Side.salad.intake,
    Red.wine.intake
  ) %>% 
  # Step 2: Calculate foodA and foodB from the original intake variables
  # Missing values are treated as 0
  mutate(
    foodA_raw = rowSums(across(
      c(
        Pure.fruit.vegetable.juice.intake,
        Olives.intake,
        Broccoli.intake,
        Celery.intake,
        Side.salad.intake,
        Red.wine.intake
      ), ~ replace_na(.x, 0)
    )),
    
    foodB_raw = rowSums(across(
      c(
        Pure.fruit.vegetable.juice.intake,
        Olives.intake,
        Broccoli.intake,
        Celery.intake,
        Side.salad.intake
      ), ~ replace_na(.x, 0)
    ))
  ) %>%
  # Step 3: Categorise each original intake variable
  # Missing values are coded as 0
  mutate(across(
    .cols = c(
      Pure.fruit.vegetable.juice.intake,
      Olives.intake,
      Broccoli.intake,
      Celery.intake,
      Side.salad.intake,
      Red.wine.intake
    ),
    .fns = ~ case_when(
      is.na(.x) | .x == 0 ~ 0,
      .x > 0 & .x <= 1 ~ 1,
      .x > 1 ~ 2
    ),
    .names = "{.col}_cat"
  )) %>%
  # Step 4: Categorise foodA_raw and foodB_raw
  mutate(
    foodA_cat = case_when(
      foodA_raw == 0 ~ 0,
      foodA_raw > 0 & foodA_raw <= 1 ~ 1,
      foodA_raw > 1 ~ 2
    ),
    foodB_cat = case_when(
      foodB_raw == 0 ~ 0,
      foodB_raw > 0 & foodB_raw <= 1 ~ 1,
      foodB_raw > 1 ~ 2
    )
  ) %>%
  # Step 5: Keep participants included in the main analysis dataset
  filter(Participant.ID %in% polyphenolwithhealth$Participant.ID)


dietaryscorewithquartile <- merge(
  polyphenolwithhealth,
  df[, c(1, 10:17)],
  by = "Participant.ID"
)

# ---- Rmd chunk ----
library(rms)

#==== 0) Data ====
df <- polyphenolwithhealth

# Complete-case analysis, or replace this with the prespecified imputation pipeline
df <- df %>% drop_na(survival.time, event, flavones)

# Optional transformations
# df$flavones <- log1p(df$flavones)
# df$flavanones <- log1p(df$flavanones)

#==== 1) rms setup ====
dd <- datadist(df)
options(datadist = "dd")

# Choose knots. Harrell commonly recommends 4–5 knots.
# Here, 5 knots are placed at selected percentiles.
kn_5 <- quantile(
  df$flavones,
  probs = c(0.05, 0.275, 0.50, 0.725, 0.95),
  na.rm = TRUE
)

#==== 2) Fit Cox model with restricted cubic splines ====
fit_rcs <- cph(
  Surv(survival.time, event) ~ 
    rcs(flavones, knots = kn_5) +
    Sex + Age.group + MET.Quintiles +
    Smoking.status...Instance.0 + alcohol.category +
    multimorbidity + IMD.quartile + Qualification +
    Energy + Free.sugar + meat + wholegrain + BMI.Category,
  data = df,
  x = TRUE,
  y = TRUE
)

#==== 3) Test overall and nonlinear associations ====
# In the anova output, the "Nonlinear" row provides the P value for nonlinearity.
anova(fit_rcs)

#==== 4) Plot the HR curve using a chosen reference value ====
# The 20th percentile is used as the reference value.
ref_val <- as.numeric(quantile(df$flavones, 0.20, na.rm = TRUE))

# Obtain spline predictions across the full range on the log-HR scale
pred_lin <- Predict(fit_rcs, flavones, conf.int = 0.95)

# Obtain the predicted value at the reference point
ref_lin <- Predict(fit_rcs, flavones = ref_val, conf.int = 0.95)

# Centre the curve at the reference point and convert estimates to HRs
pred_df <- as.data.frame(pred_lin)
ref_y <- as.numeric(ref_lin$yhat)

pred_df <- pred_df %>%
  mutate(
    HR = exp(yhat - ref_y),
    HRlo = exp(lower - ref_y),
    HRhi = exp(upper - ref_y)
  )

# Create a knot data frame for marking knot locations on the curve
knot_df <- data.frame(flavones = kn_5) %>%
  rowwise() %>%
  mutate(
    flavones_nearest = pred_df$flavones[which.min(abs(pred_df$flavones - flavones))],
    HR = pred_df$HR[which.min(abs(pred_df$flavones - flavones))],
    HRlo = pred_df$HRlo[which.min(abs(pred_df$flavones - flavones))],
    HRhi = pred_df$HRhi[which.min(abs(pred_df$flavones - flavones))]
  ) %>%
  ungroup()

pA_overall <- 0.012
pA_nonlinear <- 0.006

# Helper function to format P values
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

label_A <- bquote(
  italic(P)[overall] == .(fmt_p(pA_overall)) ~ "," ~
    italic(P)[nonlinear] == .(fmt_p(pA_nonlinear))
)
label_A_str <- as.character(as.expression(label_A))

# Plot the spline curve and mark knot locations
res_A <- ggplot(pred_df, aes(x = flavones, y = HR)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = HRlo, ymax = HRhi), alpha = 0.2) +
  geom_vline(xintercept = ref_val, linetype = 2, color = "grey50") +
  geom_point(
    data = knot_df,
    aes(x = flavones_nearest, y = HR),
    color = "red",
    size = 2.5
  ) +
  geom_text(
    data = knot_df,
    aes(
      x = flavones_nearest,
      y = HR,
      label = paste0("K", 1:nrow(knot_df))
    ),
    vjust = -1,
    color = "red",
    size = 3.2
  ) +
  labs(
    x = "Flavone intake (mg/day)",
    y = "Hazard Ratio"
  ) +
  theme_classic(base_size = 13) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_A_str,
    hjust = 1.05,
    vjust = 1.5,
    size = 4,
    parse = TRUE
  )


library(dplyr)
library(rms)
library(survival)
library(ggplot2)

#===========================
# 0. Data preparation
#===========================
df2 <- polyphenolwithhealth %>% 
  drop_na(survival.time, event, flavones, PRS.tertile)

#===========================
# 1. rms setup
#===========================
dd <- datadist(df)
options(datadist = "dd")

# Use common knot locations across PRS groups
kn <- quantile(
  df2$flavones,
  probs = c(0.05, 0.275, 0.50, 0.725, 0.95)
)

# Define the reference point
ref_val <- quantile(df2$flavones, 0.20, na.rm = TRUE)

#===========================
# 2. Function to fit a spline model within each PRS group
#===========================
fit_spline_group <- function(data_sub) {
  fit <- cph(
    Surv(survival.time, event) ~ 
      rcs(flavones, knots = kn) +
      Sex + Age.group + MET.Quintiles +
      Smoking.status...Instance.0 + alcohol.category +
      multimorbidity + IMD.quartile + Qualification +
      Energy + Free.sugar + meat + wholegrain + BMI.Category,
    data = data_sub,
    x = TRUE,
    y = TRUE
  )
  
  # Predict the HR curve
  pred_lin <- Predict(fit, flavones, conf.int = 0.95)
  pred_df <- as.data.frame(pred_lin)
  
  # Centre predictions at the reference value
  ref_hat <- Predict(fit, flavones = ref_val)$yhat
  
  pred_df <- pred_df %>%
    mutate(
      HR = exp(yhat - ref_hat),
      HRlo = exp(lower - ref_hat),
      HRhi = exp(upper - ref_hat)
    )
  
  pred_df
}

#===========================
# 3. Run spline models within each PRS group
#===========================
spl_low <- fit_spline_group(df2 %>% filter(PRS.tertile == "Low"))
spl_mid <- fit_spline_group(df2 %>% filter(PRS.tertile == "Medium"))
spl_high <- fit_spline_group(df2 %>% filter(PRS.tertile == "High"))

#===========================
# 4. Combine prediction datasets
#===========================
df_plot <- bind_rows(
  spl_low %>% mutate(PRS = "Low PRS"),
  spl_mid %>% mutate(PRS = "Medium PRS"),
  spl_high %>% mutate(PRS = "High PRS")
)

#===========================
# 5. Create P-value labels using bquote
#===========================
make_label <- function(overall, nonlinear) {
  as.character(
    as.expression(
      bquote(
        italic(P)[overall] == .(overall) ~ "," ~
          italic(P)[nonlinear] == .(nonlinear)
      )
    )
  )
}

label_low <- make_label(0.86, 0.73)
label_medium <- make_label(0.27, 0.16)
label_high <- make_label(0.003, 0.002)

#===========================
# 6. Plot PRS-stratified spline curves
#===========================
df_plot$PRS <- factor(
  df_plot$PRS,
  levels = c("High PRS", "Medium PRS", "Low PRS")
)

ggplot(df_plot, aes(x = flavones, y = HR, color = PRS)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = HRlo, ymax = HRhi, fill = PRS),
    alpha = 0.1,
    color = NA,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = ref_val, linetype = 2, color = "grey50") +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c("#0072B2", "#009E73", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#009E73", "#D55E00")) +
  
  # Set the y-axis display range
  coord_cartesian(ylim = c(0, 5)) +
  
  labs(
    x = "Flavone intake (mg/day)",
    y = "Hazard Ratio",
    color = "PRS tertile"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = as.character(as.expression(
      bquote(
        "Low PRS:" ~
          italic(P)[overall] == 0.88 ~ "," ~
          italic(P)[nonlinear] == 0.75
      )
    )),
    hjust = 1.05,
    vjust = 4.5,
    parse = TRUE,
    size = 3
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = as.character(as.expression(
      bquote(
        "Medium PRS:" ~
          italic(P)[overall] == 0.24 ~ "," ~
          italic(P)[nonlinear] == 0.14
      )
    )),
    hjust = 1.05,
    vjust = 3.0,
    parse = TRUE,
    size = 3
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = as.character(as.expression(
      bquote(
        "High PRS:" ~
          italic(P)[overall] == 0.004 ~ "," ~
          italic(P)[nonlinear] == 0.002
      )
    )),
    hjust = 1.05,
    vjust = 1.5,
    parse = TRUE,
    size = 3
  )

# ---- Rmd chunk ----
#==== 0) Data ====
df1 <- merge(df[, c(1, 8, 9)], polyphenolwithhealth, by = "Participant.ID")

# Complete-case analysis, or replace this with the prespecified imputation pipeline
df1 <- df1 %>% drop_na(survival.time, event, foodB_raw)

# Optional transformations
# df$foodB_raw <- log1p(df$foodB_raw)
# df$flavanones <- log1p(df$flavanones)

#==== 1) rms setup ====
dd <- datadist(df1)
options(datadist = "dd")

# Choose knots. Harrell commonly recommends 4–5 knots.
kn_5 <- quantile(
  df$foodB_raw,
  probs = c(0.05, 0.275, 0.50, 0.725, 0.95),
  na.rm = TRUE
)

#==== 2) Fit Cox model with restricted cubic splines ====
fit_rcs <- cph(
  Surv(survival.time, event) ~ 
    rcs(foodB_raw, knots = kn_5) +
    Sex + Age.group + MET.Quintiles +
    Smoking.status...Instance.0 + alcohol.category +
    multimorbidity + IMD.quartile + Qualification +
    Energy + Free.sugar + meat + wholegrain + BMI.Category,
  data = df1,
  x = TRUE,
  y = TRUE
)

#==== 3) Test overall and nonlinear associations ====
# In the anova output, the "Nonlinear" row provides the P value for nonlinearity.
anova(fit_rcs)

#==== 4) Plot the HR curve using a chosen reference value ====
# The 20th percentile is used as the reference value.
ref_val <- as.numeric(quantile(df$foodB_raw, 0.20, na.rm = TRUE))

# Obtain spline predictions across the full range on the log-HR scale
pred_lin <- Predict(fit_rcs, foodB_raw, conf.int = 0.95)

# Obtain the predicted value at the reference point
ref_lin <- Predict(fit_rcs, foodB_raw = ref_val, conf.int = 0.95)

# Centre the curve at the reference point and convert estimates to HRs
pred_df <- as.data.frame(pred_lin)
ref_y <- as.numeric(ref_lin$yhat)

pred_df <- pred_df %>%
  mutate(
    HR = exp(yhat - ref_y),
    HRlo = exp(lower - ref_y),
    HRhi = exp(upper - ref_y)
  )

# Create a knot data frame for marking knot locations on the curve
knot_df <- data.frame(foodB_raw = kn_5) %>%
  rowwise() %>%
  mutate(
    foodB_raw_nearest = pred_df$foodB_raw[which.min(abs(pred_df$foodB_raw - foodB_raw))],
    HR = pred_df$HR[which.min(abs(pred_df$foodB_raw - foodB_raw))],
    HRlo = pred_df$HRlo[which.min(abs(pred_df$foodB_raw - foodB_raw))],
    HRhi = pred_df$HRhi[which.min(abs(pred_df$foodB_raw - foodB_raw))]
  ) %>%
  ungroup()

pA_overall <- 0.023
pA_nonlinear <- 0.012

# Helper function to format P values
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

label_A <- bquote(
  italic(P)[overall] == .(fmt_p(pA_overall)) ~ "," ~
    italic(P)[nonlinear] == .(fmt_p(pA_nonlinear))
)
label_A_str <- as.character(as.expression(label_A))

# Plot the spline curve and mark knot locations
res_A <- ggplot(pred_df, aes(x = foodB_raw, y = HR)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = HRlo, ymax = HRhi), alpha = 0.2) +
  geom_vline(xintercept = ref_val, linetype = 2, color = "grey50") +
  geom_point(
    data = knot_df,
    aes(x = foodB_raw_nearest, y = HR),
    color = "red",
    size = 2.5
  ) +
  geom_text(
    data = knot_df,
    aes(
      x = foodB_raw_nearest,
      y = HR,
      label = paste0("K", 1:nrow(knot_df))
    ),
    vjust = -1,
    color = "red",
    size = 3.2
  ) +
  labs(
    x = "Food intake (servings/day)",
    y = "Hazard Ratio",
    title = "B. Flavone-rich food"
  ) +
  theme_classic(base_size = 13) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_A_str,
    hjust = 1.05,
    vjust = 1.5,
    size = 4,
    parse = TRUE
  )


df1 <- merge(df[, c(1, 8, 9)], polyphenolwithhealth, by = "Participant.ID")

# Complete-case analysis, or replace this with the prespecified imputation pipeline
df1 <- df1 %>% drop_na(survival.time, event, foodA_raw)

# Optional transformations
# df$foodA_raw <- log1p(df$foodA_raw)
# df$flavanones <- log1p(df$flavanones)

#==== 1) rms setup ====
dd <- datadist(df1)
options(datadist = "dd")

# Choose knots. Harrell commonly recommends 4–5 knots.
kn_5 <- quantile(
  df$foodA_raw,
  probs = c(0.05, 0.275, 0.50, 0.725, 0.95),
  na.rm = TRUE
)

#==== 2) Fit Cox model with restricted cubic splines ====
fit_rcs <- cph(
  Surv(survival.time, event) ~ 
    rcs(foodA_raw, knots = kn_5) +
    Sex + Age.group + MET.Quintiles +
    Smoking.status...Instance.0 + alcohol.category +
    multimorbidity + IMD.quartile + Qualification +
    Energy + Free.sugar + meat + wholegrain + BMI.Category,
  data = df1,
  x = TRUE,
  y = TRUE
)

#==== 3) Test overall and nonlinear associations ====
# In the anova output, the "Nonlinear" row provides the P value for nonlinearity.
anova(fit_rcs)

#==== 4) Plot the HR curve using a chosen reference value ====
# The 20th percentile is used as the reference value.
ref_val <- as.numeric(quantile(df$foodA_raw, 0.20, na.rm = TRUE))

# Obtain spline predictions across the full range on the log-HR scale
pred_lin <- Predict(fit_rcs, foodA_raw, conf.int = 0.95)

# Obtain the predicted value at the reference point
ref_lin <- Predict(fit_rcs, foodA_raw = ref_val, conf.int = 0.95)

# Centre the curve at the reference point and convert estimates to HRs
pred_df <- as.data.frame(pred_lin)
ref_y <- as.numeric(ref_lin$yhat)

pred_df <- pred_df %>%
  mutate(
    HR = exp(yhat - ref_y),
    HRlo = exp(lower - ref_y),
    HRhi = exp(upper - ref_y)
  )

# Create a knot data frame for marking knot locations on the curve
knot_df <- data.frame(foodA_raw = kn_5) %>%
  rowwise() %>%
  mutate(
    foodA_raw_nearest = pred_df$foodA_raw[which.min(abs(pred_df$foodA_raw - foodA_raw))],
    HR = pred_df$HR[which.min(abs(pred_df$foodA_raw - foodA_raw))],
    HRlo = pred_df$HRlo[which.min(abs(pred_df$foodA_raw - foodA_raw))],
    HRhi = pred_df$HRhi[which.min(abs(pred_df$foodA_raw - foodA_raw))]
  ) %>%
  ungroup()

pA_overall <- 0.14
pA_nonlinear <- 0.08

# Helper function to format P values
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

label_A <- bquote(
  italic(P)[overall] == .(fmt_p(pA_overall)) ~ "," ~
    italic(P)[nonlinear] == .(fmt_p(pA_nonlinear))
)
label_A_str <- as.character(as.expression(label_A))

# Plot the spline curve and mark knot locations
res_B <- ggplot(pred_df, aes(x = foodA_raw, y = HR)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = HRlo, ymax = HRhi), alpha = 0.2) +
  geom_vline(xintercept = ref_val, linetype = 2, color = "grey50") +
  geom_point(
    data = knot_df,
    aes(x = foodA_raw_nearest, y = HR),
    color = "red",
    size = 2.5
  ) +
  geom_text(
    data = knot_df,
    aes(
      x = foodA_raw_nearest,
      y = HR,
      label = paste0("K", 1:nrow(knot_df))
    ),
    vjust = -1,
    color = "red",
    size = 3.2
  ) +
  labs(
    x = "Food intake (servings/day)",
    y = "Hazard Ratio",
    title = "A. Flavone-rich food and wine"
  ) +
  theme_classic(base_size = 13) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_A_str,
    hjust = 1.05,
    vjust = 1.5,
    size = 4,
    parse = TRUE
  )

library(cowplot)

combined_plot <- plot_grid(
  res_B,
  res_A,
  ncol = 2,
  rel_widths = c(2, 2)
)

combined_plot


# ---- Rmd chunk ----
# Load the data
table3 <- read.csv("diet quality output/table2.csv")
table3$Person.years <- round(table3$Person.years, 0)
table3$Person.years <- ifelse(is.na(table3$Person.years), "", table3$Person.years)

# Prepare columns for the forest plot
table3 <- table3 %>%
  mutate(
    "Flavonoid subclass (mg/day)" = Flavonoid.subclass..mg.d.,
    "Cases/Total (Incidence %)" = Case.total..incidence...,
    "Person years" = Person.years,
    `Model 1;HR (95% CI)` = Model.1..HR..95..CI.,
    `Model 2;HR (95% CI)` = Model.2..HR..95..CI.,
    `Model 3;HR (95% CI)` = Model.3..HR..95..CI.
  )

# Parse HRs and 95% CIs from formatted text columns
parse_ci <- function(ci_vec) {
  strcapture(
    pattern = "^\\s*([0-9.+-eE]+)\\s*\\(\\s*([0-9.+-eE]+)\\s*-\\s*([0-9.+-eE]+)\\s*\\)\\s*$",
    x = ci_vec,
    proto = list(HR = double(), lower = double(), upper = double())
  )
}

table3 <- cbind(
  table3,
  parse_ci(table3$`Model 1;HR (95% CI)`),
  parse_ci(table3$`Model 2;HR (95% CI)`),
  parse_ci(table3$`Model 3;HR (95% CI)`)
)

names(table3) <- make.unique(names(table3))
names(table3)

table3$Model1 <- paste(rep(" ", 22), collapse = " ")
table3$Model2 <- paste(rep(" ", 22), collapse = " ")
table3$Model3 <- paste(rep(" ", 22), collapse = " ")

table3 <- table3 %>%
  mutate(
    HR = as.numeric(HR),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    HR.1 = as.numeric(HR.1),
    lower.1 = as.numeric(lower.1),
    upper.1 = as.numeric(upper.1),
    HR.2 = as.numeric(HR.2),
    lower.2 = as.numeric(lower.2),
    upper.2 = as.numeric(upper.2)
  )

library(knitr)
library(forestploter)
library(grid)

# Create the initial forest plot
plot3 <- forest(
  table3[, c(7, 8, 9, 10, 22, 11, 23, 12, 24)],
  est = list(table3$HR, table3$HR.1, table3$HR.2),
  lower = list(table3$lower, table3$lower.1, table3$lower.2),
  upper = list(table3$upper, table3$upper.1, table3$upper.2),
  ci_column = c(5, 7, 9),
  ref.line = 1
)

# Customise plot appearance
tm <- forest_theme(
  base_size = 10,
  ci_pch = 15,
  ci_col = "#762a83",
  ci_fill = "blue",
  ci_alpha = 0.8,
  ci_lty = 1,
  ci_lwd = 1.5,
  ci_Theight = 0.2,
  refline_gp = gpar(col = "grey20", lwd = 1, lty = "dashed")
)

# Apply theme and set axis limits
plot3 <- forest(
  table3[, c(7, 8, 9, 10, 22, 11, 23, 12, 24)],
  est = list(table3$HR, table3$HR.1, table3$HR.2),
  lower = list(table3$lower, table3$lower.1, table3$lower.2),
  upper = list(table3$upper, table3$upper.1, table3$upper.2),
  ci_column = c(5, 7, 9),
  ref_line = 1,
  xlim = c(0.50, 1.50),
  ticks_at = c(0.55, 0.70, 0.85, 1, 1.15, 1.30, 1.45),
  theme = tm
)

plot3

ggsave(
  filename = "figure2.tiff",
  plot = plot3,
  dpi = 900,
  width = 15,
  height = 9
)


# Load the data
table3 <- read.csv("diet quality output/table3.csv")
table3$Person.years <- round(table3$Person.years, 0)
table3$Person.years <- ifelse(is.na(table3$Person.years), "", table3$Person.years)

# Prepare columns for the forest plot
table3 <- table3 %>%
  mutate(
    "Flavone food" = Food,
    "Cases/Total (Incidence %)" = Case.total..incidence...,
    "Person years" = Person.years,
    `Model 1;HR (95% CI)` = Model.1..HR..95..CI.,
    `Model 2;HR (95% CI)` = Model.2..HR..95..CI.,
    `Model 3;HR (95% CI)` = Model.3..HR..95..CI.
  )

# Parse HRs and 95% CIs from formatted text columns
parse_ci <- function(ci_vec) {
  strcapture(
    pattern = "^\\s*([0-9.+-eE]+)\\s*\\(\\s*([0-9.+-eE]+)\\s*-\\s*([0-9.+-eE]+)\\s*\\)\\s*$",
    x = ci_vec,
    proto = list(HR = double(), lower = double(), upper = double())
  )
}

table3 <- cbind(
  table3,
  parse_ci(table3$`Model 1;HR (95% CI)`),
  parse_ci(table3$`Model 2;HR (95% CI)`),
  parse_ci(table3$`Model 3;HR (95% CI)`)
)

names(table3) <- make.unique(names(table3))
names(table3)

table3$Model1 <- paste(rep(" ", 22), collapse = " ")
table3$Model2 <- paste(rep(" ", 22), collapse = " ")
table3$Model3 <- paste(rep(" ", 22), collapse = " ")

table3 <- table3 %>%
  mutate(
    HR = as.numeric(HR),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    HR.1 = as.numeric(HR.1),
    lower.1 = as.numeric(lower.1),
    upper.1 = as.numeric(upper.1),
    HR.2 = as.numeric(HR.2),
    lower.2 = as.numeric(lower.2),
    upper.2 = as.numeric(upper.2)
  )

library(knitr)
library(forestploter)
library(grid)

# Create the initial forest plot
plot3 <- forest(
  table3[, c(7, 8, 9, 10, 22, 11, 23, 12, 24)],
  est = list(table3$HR, table3$HR.1, table3$HR.2),
  lower = list(table3$lower, table3$lower.1, table3$lower.2),
  upper = list(table3$upper, table3$upper.1, table3$upper.2),
  ci_column = c(5, 7, 9),
  ref.line = 1
)

# Customise plot appearance
tm <- forest_theme(
  base_size = 10,
  ci_pch = 15,
  ci_col = "#762a83",
  ci_fill = "blue",
  ci_alpha = 0.8,
  ci_lty = 1,
  ci_lwd = 1.5,
  ci_Theight = 0.2,
  refline_gp = gpar(col = "grey20", lwd = 1, lty = "dashed")
)

# Apply theme and set axis limits
plot3 <- forest(
  table3[, c(7, 8, 9, 10, 22, 11, 23, 12, 24)],
  est = list(table3$HR, table3$HR.1, table3$HR.2),
  lower = list(table3$lower, table3$lower.1, table3$lower.2),
  upper = list(table3$upper, table3$upper.1, table3$upper.2),
  ci_column = c(5, 7, 9),
  ref_line = 1,
  xlim = c(0.50, 1.50),
  ticks_at = c(0.55, 0.70, 0.85, 1, 1.15, 1.30, 1.45),
  theme = tm
)

plot3

ggsave(
  filename = "figure3.tiff",
  plot = plot3,
  dpi = 900,
  width = 15,
  height = 10
)

# Suggested figure size: 1500 × 900 pixels


# ---- Rmd chunk ----
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(cowplot)

# ============================================================
# 1. Create the flavone factor variable using quartiles 1–4
# ============================================================
dietaryscorewithquartile$flavone <- factor(
  dietaryscorewithquartile$flavones_ord4,
  levels = c(1, 2, 3, 4),
  labels = c("Q1", "Q2", "Q3", "Q4")
)

# ============================================================
# 2. Fit the Cox model
# x = TRUE and model = TRUE are required for adjusted survival prediction.
# ============================================================
model.uni <- coxph(
  Survival.Obj ~ flavone + Sex + Age.group + MET.Quintiles +
    Smoking.status...Instance.0 + alcohol.category + multimorbidity +
    IMD.quartile + Qualification + Energy + Free.sugar + meat +
    wholegrain + BMI.Category,
  data = dietaryscorewithquartile,
  ties = "efron",
  x = TRUE,
  model = TRUE
)

# ============================================================
# 3. Function to identify the modal category for factor covariates
# ============================================================
mode_ <- function(x) {
  names(which.max(table(x)))
}

# ============================================================
# 4. Build the reference covariate row
# ============================================================
ref_row <- data.frame(
  Sex = mode_(dietaryscorewithquartile$Sex),
  Age.group = mode_(dietaryscorewithquartile$Age.group),
  MET.Quintiles = mode_(dietaryscorewithquartile$MET.Quintiles),
  Smoking.status...Instance.0 = mode_(dietaryscorewithquartile$Smoking.status...Instance.0),
  alcohol.category = mode_(dietaryscorewithquartile$alcohol.category),
  multimorbidity = mode_(dietaryscorewithquartile$multimorbidity),
  IMD.quartile = mode_(dietaryscorewithquartile$IMD.quartile),
  Qualification = mode_(dietaryscorewithquartile$Qualification),
  
  # Use medians for numeric covariates
  Energy = median(dietaryscorewithquartile$Energy, na.rm = TRUE),
  Free.sugar = median(dietaryscorewithquartile$Free.sugar, na.rm = TRUE),
  meat = median(dietaryscorewithquartile$meat, na.rm = TRUE),
  wholegrain = median(dietaryscorewithquartile$wholegrain, na.rm = TRUE),
  
  BMI.Category = mode_(dietaryscorewithquartile$BMI.Category)
)

# ============================================================
# 5. Create four copies of the reference row, one for each flavone quartile
# ============================================================
nd_list <- replicate(4, ref_row, simplify = FALSE)

# Combine into a single data frame
newdata <- do.call(rbind, nd_list)

# Add the flavone quartile variable
newdata$flavone <- factor(
  c("Q1", "Q2", "Q3", "Q4"),
  levels = c("Q1", "Q2", "Q3", "Q4")
)

# ----------------------------
# 6. Match factor levels with the fitted model
# ----------------------------
for (v in names(model.uni$xlevels)) {
  if (v %in% names(newdata)) {
    newdata[[v]] <- factor(
      newdata[[v]],
      levels = model.uni$xlevels[[v]]
    )
  }
}

print(newdata)

# ----------------------------
# 7. Predict adjusted survival curves
# ----------------------------
max_time <- max(dietaryscorewithquartile$survival.time, na.rm = TRUE)

sf_adj <- survfit(
  model.uni,
  newdata = newdata,
  se.fit = TRUE,
  conf.type = "log",
  # Specify the full time sequence for prediction
  time = seq(0, max_time, by = 0.1)
)

# Inspect the survival object
print(sf_adj)

# Plot adjusted survival curves with confidence interval ribbons
p_adj <- ggsurvplot(
  sf_adj,
  data = newdata,
  conf.int = TRUE,
  # legend.title = "flavone",
  # legend = "none",
  palette = "Dark2",
  ggtheme = theme_classic(base_size = 12),
  xlab = "Time (years)",
  ylab = "Survival probability",
  ylim = c(0.976, 1.00)
)

p_adj$plot <- p_adj$plot + 
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

p_adj$plot <- p_adj$plot + labs(title = "C. High genetic risk")

# Fit Kaplan-Meier curves only for the risk table
fit_km <- survfit(Survival.Obj ~ flavone, data = dietaryscorewithquartile)

km_plot <- ggsurvplot(
  fit_km,
  data = dietaryscorewithquartile,
  risk.table = TRUE,
  risk.table.type = "nrisk_cumevents",
  risk.table.title = "Participants at risk",
  risk.table.height = 0.30,
  risk.table.fontsize = 3.6,
  risk.table.y.text = TRUE,
  risk.table.y.text.col = FALSE,
  break.time.by = 2.5,
  conf.int = FALSE,
  legend = "none",
  xlab = "Time (years)",
  ylab = NULL,
  ggtheme = theme_classic(base_size = 12),
  tables.theme = theme_cleantable()
)

final_plot3 <- plot_grid(
  p_adj$plot,
  km_plot$table,
  ncol = 1,
  rel_heights = c(0.72, 0.28)
)

print(final_plot3)

library(cowplot)

combined_plot <- plot_grid(
  final_plot1,
  final_plot2,
  final_plot3,
  nrow = 3,
  rel_widths = c(3, 3, 3)
)

combined_plot

# Suggested figure size: 1000 × 800 pixels
