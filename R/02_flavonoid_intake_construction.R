# 02_flavonoid_intake_construction.R
# Auto-refactored from the original analysis notebook: github.Rmd.
# Run via targets::tar_make() from the project root.
# NOTE: UK Biobank raw data are not included in this repository.

source(file.path("R", "00_setup.R"))
load_packages()

# ---- Rmd chunk ----
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(stringr)

# --- 1) Define the original column names included in each of the five flavonoid subclasses ---
# Columns that are not present in the dataset will be ignored automatically.
cols <- list(
  anthocyanins = c("Cyanidin", "Delphinidin", "Malvidin", "Pelargonidin", "Petunidin", "Peonidin"),
  
  flavanones = c("Eriodictyol", "Hesperetin", "Naringenin"),
  
  flavonols = c(
    "Quercetin", "Kaempferol", "Myricetin", "Isorhamnetin",
    "Jaceidin", "Patuletin", "Spinacetin"
  ),
  
  flavones = c("Luteolin", "Apigenin"),
  
  flavan3ols = c(
    "X....Catechin", "X....Epicatechin", "X....Epigallocatechin", "X....Gallocatechin",
    "X....Catechin.3.O.gallate", "X....Epicatechin.3.O.gallate",
    "X....Epigallocatechin.3.O.gallate", "X....Gallocatechin.3.O.gallate",
    "Theaflavin", "Thearubigins", "X02.mers", "X03.mers", "X04.06.mers",
    "X07.10.mers", "Polymers...10.mers."
  )
)

# --- 2) Helper function: convert individual flavonoid columns into subclass totals ---
# Input:
#   flavonoid.checkout: rows represent foods; columns represent individual flavonoids or polymers.
#
# Output:
#   flavonoid_subclass_per100g: rows represent foods; columns represent the five flavonoid subclasses.
#   Values are expressed as mg/100 g.
build_subclass_matrix <- function(df, cols_map) {
  out <- lapply(cols_map, function(cand_cols) {
    keep <- intersect(cand_cols, colnames(df))
    
    if (length(keep) == 0) {
      rep(0, nrow(df))
    } else {
      rowSums(df[, keep, drop = FALSE], na.rm = TRUE)
    }
  })
  
  out <- as.data.frame(out, check.names = FALSE)
  rownames(out) <- rownames(df)
  
  out
}


row_map <- tribble(
  ~row, ~food, ~components,
  14L, "Black tea", list(list(col="Standard.tea.intake", portion=240, weight=1)),
  15L, "Green tea", list(list(col="Green.tea.intake", portion=240, weight=1)),
  16L, "Other tea A", list(list(col="Other.tea.intake", portion=240, weight=0.5)),
  17L, "Other tea B", list(list(col="Other.tea.intake", portion=240, weight=0.5)),
  
  32L, "Dark chocolate", list(list(col="Dark.chocolate.intake", portion=50, weight=1)),
  33L, "Milk chocolate", list(list(col="Milk.chocolate.intake", portion=50, weight=1)),
  40L, "Hot chocolate (all)", list(list(col="Hot.chocolate.intake", portion=240, weight=1),
                                   list(col="Low.calorie.hot.chocolate.intake", portion=240, weight=1)),
  50L, "Chocolate biscuits", list(list(col="Chocolate.biscuits.intake", portion=13.5, weight=1)),
  51L, "Mixed chocolate snacks", list(
    list(col="Chocolate.covered.biscuits.intake", portion=13.5, weight=1),
    list(col="Chocolate.covered.raisin.intake",   portion=30,   weight=1),
    list(col="Chocolate.bar.intake",              portion=50,   weight=1),
    list(col="Chocolate.sweet.intake",            portion=30,   weight=1)
  ),
  
  58L, "Red wine",        list(list(col="Red.wine.intake",        portion=165.76, weight=1)),
  59L, "Rose wine",       list(list(col="Rose.wine.intake",       portion=163.5,  weight=1)),
  60L, "White wine",      list(list(col="White.wine.intake",      portion=162.96, weight=1)),
  61L, "Fortified wine",  list(list(col="Fortified.wine.intake",  portion=50,     weight=1)),
  
  85L, "Orange juice",           list(list(col="Orange.juice.intake",          portion=250, weight=1)),
  78L, "Grapefruit juice",       list(list(col="Grapefruit.juice.intake",      portion=250, weight=1)),
  90L, "Pure fruit/veg juice",   list(list(col="Pure.fruit.vegetable.juice.intake", portion=250, weight=1)),
  89L, "Fruit smoothie",         list(list(col="Fruit.smoothie.intake",        portion=250, weight=1)),
  
  348L, "Stewed fruit (apples)", list(list(col="Stewed.fruit.intake", portion=80, weight=1)),
  298L, "Prunes",                list(list(col="Prune.intake",        portion=30, weight=1)),
  366L, "Dried fruit",           list(list(col="Dried.fruit.intake",  portion=30, weight=1)),
  304L, "Mixed fruit",           list(list(col="Mixed.fruit.intake",  portion=80, weight=1)),
  252L, "Apples",                list(list(col="Apple.intake",        portion=80, weight=1)),
  300L, "Bananas",               list(list(col="Banana.intake",       portion=80, weight=1)),
  354L, "Berries",               list(list(col="Berry.intake",        portion=80, weight=1)),
  346L, "Cherries",              list(list(col="Cherry.intake",       portion=80, weight=1)),
  369L, "Grapefruit",            list(list(col="Grapefruit.intake",   portion=80, weight=1)),
  
  # grape two type 
  377L, "Grapes - black", list(list(col="Grape.intake", portion=80, weight=0.32)),
  378L, "Grapes - green", list(list(col="Grape.intake", portion=80, weight=0.68)),
  
  307L, "Mango",      list(list(col="Mango.intake",      portion=80, weight=1)),
  311L, "Melon",      list(list(col="Melon.intake",      portion=80, weight=1)),
  372L, "Orange",     list(list(col="Orange.intake",     portion=80, weight=1)),
  373L, "Tangerine",  list(list(col="Satsuma.intake",    portion=80, weight=1)),
  340L, "Peach/Nectarine", list(list(col="Peach.nectarine.intake", portion=80, weight=1)),
  255L, "Pear",       list(list(col="Pear.intake",       portion=80, weight=1)),
  319L, "Pineapple",  list(list(col="Pineapple.intake",  portion=80, weight=1)),
  344L, "Plum",       list(list(col="Plum.intake",       portion=40, weight=1)),
  321L, "Other fruit - pomegranate", list(list(col="Other.fruit.intake", portion=80, weight=0.5)),
  305L, "Other fruit - kiwi",        list(list(col="Other.fruit.intake", portion=80, weight=0.5)),
  
  578L, "Fried potatoes",          list(list(col="Fried.potatoes.intake",          portion=180, weight=1)),
  572L, "Boiled/Baked potatoes",   list(list(col="Boiled.baked.potatoes.intake",   portion=180, weight=1)),
  576L, "Mashed potatoes",         list(list(col="Mashed.potato.intake",           portion=180, weight=1)),
  
  614L, "Mixed vegetables",  list(list(col="Mixed.vegetable.intake", portion=80, weight=1)),
  615L, "Vegetable pieces",  list(list(col="Vegetable.pieces.intake", portion=80, weight=1)),
  623L, "Coleslaw",          list(list(col="Coleslaw.intake",         portion=80, weight=1)),
  667L, "Side salad",        list(list(col="Side.salad.intake",       portion=80, weight=1)),
  379L, "Avocado",           list(list(col="Avocado.intake",          portion=160, weight=1)),
  610L, "Green beans",       list(list(col="Green.bean.intake",       portion=80, weight=0.93)),
  587L, "Beetroot",          list(list(col="Beetroot.intake",         portion=80, weight=1)),
  505L, "Broccoli",          list(list(col="Broccoli.intake",         portion=80, weight=1)),
  565L, "Butternut squash",  list(list(col="Butternut.squash.intake", portion=80, weight=1)),
  503L, "Cabbage",           list(list(col="Cabbage.kale.intake",     portion=80, weight=0.5)),
  500L, "Kale",              list(list(col="Cabbage.kale.intake",     portion=80, weight=0.5)),
  592L, "Carrots",           list(list(col="Carrot.intake",           portion=80, weight=1)),
  499L, "Cauliflower",       list(list(col="Cauliflower.intake",      portion=80, weight=1)),
  465L, "Celery",            list(list(col="Celery.intake",           portion=26.7, weight=1)),
  559L, "Courgette",         list(list(col="Courgette.intake",        portion=80, weight=1)),
  467L, "Cucumber",          list(list(col="Cucumber.intake",         portion=32, weight=1)),
  419L, "Garlic",            list(list(col="Garlic.intake",           portion=5,  weight=1)),
  522L, "Leek",              list(list(col="Leek.intake",             portion=80, weight=1)),
  702L, "Lettuce",           list(list(col="Lettuce.intake",          portion=80, weight=1)),
  712L, "Mushroom",          list(list(col="Mushroom.intake",         portion=80, weight=1)),
  521L, "Onion",             list(list(col="Onion.intake",            portion=80, weight=1)),
  586L, "Parsnip",           list(list(col="Parsnip.intake",          portion=80, weight=1)),
  397L, "Sweet pepper",      list(list(col="Sweet.pepper.intake",     portion=80, weight=1)),
  # spinach cooked and raw
  513L, "Spinach (cooked)",  list(list(col="Spinach.intake",          portion=80, weight=0.5)),
  706L, "Spinach (raw)",     list(list(col="Spinach.intake",          portion=80, weight=0.5)),
  507L, "Sprouts",           list(list(col="Sprouts.intake",          portion=80, weight=1)),
  139L, "Sweetcorn",         list(list(col="Sweetcorn.intake",        portion=80, weight=1)),
  571L, "Sweet potato",      list(list(col="Sweet.potato.intake",     portion=130, weight=1)),
  392L, "Tomato (fresh)",    list(list(col="Fresh.tomato.intake",     portion=80, weight=1)),
  390L, "Tomato (tinned)",   list(list(col="Tinned.tomato.intake",    portion=200, weight=1)),
  597L, "Turnip",            list(list(col="Turnip.swede.intake",     portion=80, weight=0.5)),
  589L, "Swede",             list(list(col="Turnip.swede.intake",     portion=80, weight=0.5)),
  695L, "Watercress",        list(list(col="Watercress.intake",       portion=80, weight=1)),
  564L, "Asparagus",         list(list(col="Other.vegetables.intake",        portion=80, weight=0.4)),
  394L, "Aubergine",         list(list(col="Other.vegetables.intake",        portion=80, weight=0.3)),
  566L, "Squash (other)",    list(list(col="Other.vegetables.intake",           portion=80, weight=0.3)),
  
  847L, "Olives - black",    list(list(col="Olives.intake",           portion=50, weight=0.5)),
  848L, "Olives - green",    list(list(col="Olives.intake",           portion=50, weight=0.5))
)

# 4) Build the five flavonoid subclass matrix (mg/100 g)
flavonoid_subclass_per100g <- build_subclass_matrix(flavonoid.checkout, cols)

# 5) Check for missing columns (optional)
missing_cols <- lapply(cols, function(v) setdiff(v, colnames(flavonoid.checkout)))

# 6) Functions for mean-based intake contribution (new)

mean_intake_col <- function(df, col, na_to_zero = TRUE) {
  if (!col %in% names(df)) return(NA_real_)
  v <- df[[col]]
  if (na_to_zero) v[is.na(v)] <- 0
  mean(v, na.rm = TRUE)
}

compute_row_mean_contrib <- function(row, na_to_zero = TRUE) {
  stopifnot(row >= 1, row <= nrow(flavonoid_subclass_per100g))
  
  vec <- as.numeric(flavonoid_subclass_per100g[row, ])
  names(vec) <- colnames(flavonoid_subclass_per100g)
  
  cfg <- dplyr::filter(row_map, row == !!row)
  if (nrow(cfg) == 0) return(setNames(rep(0, length(vec)), names(vec)))
  
  out <- setNames(rep(0, length(vec)), names(vec))
  
  for (comp in cfg$components[[1]]) {
    colname <- comp$col
    portion <- comp$portion
    weight  <- comp$weight
    
    m_intake <- mean_intake_col(p.24h.foods.average, colname, na_to_zero = na_to_zero)
    
    if (is.na(m_intake)) {
      warning(sprintf("Intake column not found: %s (row = %s)", colname, row))
      next
    }
    
    grams_mean <- as.numeric(m_intake) * portion
    
    out <- out + vec * (grams_mean / 100) * weight
  }
  
  out
}

compute_mean_contrib_table <- function(rows, na_to_zero = TRUE) {
  rows <- as.integer(rows)
  
  res <- purrr::map(rows, ~{
    v <- compute_row_mean_contrib(.x, na_to_zero = na_to_zero)
    
    tibble(
      row = .x,
      subclass = names(v),
      mean_mg_day = as.numeric(v)
    )
    
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::left_join(
      dplyr::select(row_map, row, food),
      by = "row"
    ) %>%
    dplyr::relocate(food, .after = row) %>%
    dplyr::group_by(row, food) %>%
    dplyr::mutate(total_mean_mg_day = sum(mean_mg_day, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  list(
    per_food_by_subclass_mean =
      res %>%
      dplyr::arrange(
        dplyr::desc(total_mean_mg_day),
        food,
        dplyr::desc(mean_mg_day)
      ),
    
    per_food_total_mean =
      res %>%
      dplyr::distinct(row, food, total_mean_mg_day) %>%
      dplyr::arrange(dplyr::desc(total_mean_mg_day)),
    
    grand_total_mean =
      res %>%
      dplyr::distinct(row, total_mean_mg_day) %>%
      dplyr::summarise(sum(total_mean_mg_day), .groups = "drop") %>%
      dplyr::pull(1)
  )
}

# 7) Rows of interest (retain as specified)
rows_of_interest <- c(
  14, 15, 16, 17, 32, 33, 40, 50, 51, 58, 59, 60, 61,
  85, 78, 90, 89, 348, 298, 366, 304, 252, 300, 354,
  346, 369, 377, 378, 307, 311, 372, 373, 340, 255,
  319, 344, 321, 305, 578, 572, 576, 614, 615, 623,
  667, 379, 610, 587, 505, 565, 503, 500, 592, 499,
  465, 559, 467, 419, 522, 702, 712, 521, 586, 397,
  513, 706, 507, 139, 571, 392, 390, 597, 589, 695,
  564, 394, 566, 847, 848
)

# 8) Run the computation
out_mean <- compute_mean_contrib_table(rows_of_interest, na_to_zero = TRUE)

per_food_total_mean        <- out_mean$per_food_total_mean
per_food_by_subclass_mean  <- out_mean$per_food_by_subclass_mean
grand_total_mean           <- out_mean$grand_total_mean

# ---- Rmd chunk ----

### Select the variables into a separate dataset for flavonoid evaluation ###
diet.tea <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.tea[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.tea) <- c("Participant.ID", "group")


sum.tea.values <- function(tea.rows, flavonoid.data, diet.tea) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 223
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.tea (each participant)
  for (participant in seq_len(nrow(diet.tea))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified tea row and add its values to the summed.values vector
    for (i in seq_along(tea.rows)) {
      tea.row <- tea.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (tea.row == 252) {
        # Apple (row 252)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ])  * p.24h.foods.average$Apple.intake[participant]) * 80 / 100
      } else if (tea.row == 69) {
        # Red wine (row 69)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Red.wine.intake[participant]) * 165.76 / 100
      } else if (tea.row == 354) {
        # Berry (row 354)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ])  * p.24h.foods.average$Berry.intake[participant]) * 80 / 100
      } else if (tea.row == 369) {
        # Grapefruit (row 369)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Grapefruit.intake[participant]) * 80 / 100
      } else if (tea.row == 377) {
        # Grape purple (row = 377)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Grape.intake[participant]) * 80 / 100 *0.32
      } else if (tea.row == 378) {
        # Grape green (row = 378)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Grape.intake[participant]) * 80 / 100 *0.68
      } else if (tea.row == 372) {
        # Orange intake (row 372)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Orange.intake[participant]) * 80 / 100
      } else if (tea.row == 373) {
        # Tangerine intake (row 373)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Satsuma.intake[participant]) * 80 / 100
      } else if (tea.row == 521) {
        # Onion intake (row 521)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Onion.intake[participant]) * 80 / 100
      } else if (tea.row == 397) {
        # Sweet pepper intake (row 397)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Sweet.pepper.intake[participant]) * 80 / 100
      } else if (tea.row == 33) {
        # Dark.chocolate.intake (row 33)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Dark.chocolate.intake[participant]) * 50 / 100
      } else if (tea.row == 3) {
        # Black tea (row 3)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Standard.tea.intake[participant]) * 240 / 100
      } else if (tea.row == 5) {
        # Green tea (row 5)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Green.tea.intake[participant]) * 240 / 100
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.tea dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.tea[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.tea)
}

# Define the tea rows (adjust based on the actual rows you want to sum)
tea.rows <- c(252, 69, 354, 369, 377, 378, 372, 373, 521, 397, 33, 3, 5)

# Call the function to calculate and update the diet.tea dataframe, incorporating p.24h.foods.average adjustments
diet.tea <- sum.tea.values(tea.rows, flavonoid.checkout, diet.tea)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.tea, "flavonoid/test.csv", row.names = FALSE)

# ---- Rmd chunk ----
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.tea <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.tea[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.tea) <- c("Participant.ID", "group")

# Define a function to sum the tea values, adjust with p.24h.foods.average, and add to diet.tea dataframe
sum.tea.values <- function(tea.rows, flavonoid.data, diet.tea) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.tea (each participant)
  for (participant in seq_len(nrow(diet.tea))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified tea row and add its values to the summed.values vector
    for (i in seq_along(tea.rows)) {
      tea.row <- tea.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (tea.row == 14) {
        # Black tea (row 14)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Standard.tea.intake[participant]) * 240 / 100
      } else if (tea.row == 881) {
        # Rooibos tea (row 881)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Rooibos.tea.intake[participant]) * 240 / 100
      } else if (tea.row == 5) {
        # Green tea (row 5)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Green.tea.intake[participant]) * 240 / 100
      } else if (tea.row == 16) {
        # Other tea (row 16)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Other.tea.intake[participant]) * 240 / 100 /2
      } else if (tea.row == 17) {
        # Other tea (row 17)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[tea.row, ]) * p.24h.foods.average$Other.tea.intake[participant]) * 240 / 100 /2
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.tea dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.tea[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.tea)
}

# Define the tea rows (adjust based on the actual rows you want to sum)
tea.rows <- c(14, 881, 5, 16, 17)

# Call the function to calculate and update the diet.tea dataframe, incorporating p.24h.foods.average adjustments
diet.tea <- sum.tea.values(tea.rows, flavonoid.checkout, diet.tea)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.tea, "flavonoid/diet.tea.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
diet.tea$anthocyanins.tea <- diet.tea$Cyanidin + diet.tea$Delphinidin + diet.tea$Malvidin + diet.tea$Pelargonidin + diet.tea$Petunidin + diet.tea$Peonidin

# flavanones 
diet.tea$flavanones.tea <- diet.tea$Eriodictyol + diet.tea$Hesperetin + diet.tea$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.tea$flavan.3.ols.tea <- diet.tea$X....Catechin + diet.tea$X....Epicatechin + diet.tea$X....Epigallocatechin + diet.tea$X....Gallocatechin + diet.tea$X....Catechin.3.O.gallate + diet.tea$X....Epicatechin.3.O.gallate + diet.tea$X....Gallocatechin.3.O.gallate +  diet.tea$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.tea$flavonols.tea <- diet.tea$Quercetin + diet.tea$Kaempferol + diet.tea$Myricetin + diet.tea$Isorhamnetin

# flavones
diet.tea$flavones.tea <- diet.tea$Luteolin + diet.tea$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.tea$polymers.tea <-  diet.tea$X02.mers + diet.tea$X03.mers + diet.tea$X04.06.mers + diet.tea$X07.10.mers + diet.tea$Polymers...10.mers.

# Proanthocyanidins
diet.tea$proanthocyanidins.tea <- diet.tea$Theaflavin + diet.tea$Thearubigins

# total flavonoid
diet.tea$total.flavonoid.tea <- diet.tea$anthocyanins.tea + diet.tea$flavanones.tea + diet.tea$flavan.3.ols.tea + diet.tea$flavonols.tea + diet.tea$flavones + diet.tea$polymers.tea + diet.tea$proanthocyanidins.tea

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.tea$group <- factor(diet.tea$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.tea)[96:103]  # Selecting the variable names from column 3 to 217
description.tea <- CreateTableOne(data = diet.tea,
                                  vars = vars,  # Variable names extracted above
                                  strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.tea, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.tea

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.tea), "flavonoid/tearesult.csv")

# ---- Rmd chunk ----
rm(diet.coffee)
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.chocolate <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.chocolate[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.chocolate) <- c("Participant.ID", "group")

# Define a function to sum the chocolate values, adjust with p.24h.foods.average, and add to diet.chocolate dataframe
sum.chocolate.values <- function(chocolate.rows, flavonoid.data, diet.chocolate) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.chocolate (each participant)
  for (participant in seq_len(nrow(diet.chocolate))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified chocolate row and add its values to the summed.values vector
    for (i in seq_along(chocolate.rows)) {
      chocolate.row <- chocolate.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (chocolate.row == 41) {
        # Milk.chocolate.intake (row 41)
        summed.values <- summed.values +
          (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Milk.chocolate.intake[participant]) * 50 / 100
      } else if (chocolate.row == 42) {
        # Dark.chocolate.intake (row 42)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Dark.chocolate.intake[participant]) * 50 / 100
      } else if (chocolate.row == 40) {
        # (Low) Hot.chocolate.intake (row 40)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Hot.chocolate.intake[participant]) * 240 / 100 + (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Low.calorie.hot.chocolate.intake[participant]) * 240 / 100
      } else if (chocolate.row == 50) {
        # Chocolate.biscuits.intake(row 50)
        summed.values <- summed.values + 
          (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Chocolate.biscuits.intake[participant]) * 13.5 / 100
      } else if (chocolate.row == 51) {
        # Chocolate.covered.raisin.intake and Chocolate.covered biscuits intake, chocolate bar, Chocolate.sweet.intake (row 51)
        summed.values <- summed.values + 
          (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Chocolate.covered.biscuits.intake[participant]) *13.5 / 100 + (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Chocolate.covered.raisin.intake[participant]) * 30 / 100 + (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Chocolate.bar.intake[participant]) * 50 / 100 + (as.numeric(flavonoid.columns[chocolate.row, ]) * p.24h.foods.average$Chocolate.sweet.intake[participant]) * 30 / 100
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.chocolate dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.chocolate[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.chocolate)
}

# Define the chocolate rows (adjust based on the actual rows you want to sum)
chocolate.rows <- c(41, 42, 40, 50, 51)

# Call the function to calculate and update the diet.chocolate dataframe, incorporating p.24h.foods.average adjustments
diet.chocolate <- sum.chocolate.values(chocolate.rows, flavonoid.checkout, diet.chocolate)
diet.chocolate$X....Epicatechin[is.na(diet.chocolate$X....Epicatechin)] <- 0
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.chocolate, "flavonoid/diet.chocolate.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
diet.chocolate$anthocyanins.chocolate <- diet.chocolate$Cyanidin + diet.chocolate$Delphinidin + diet.chocolate$Malvidin + diet.chocolate$Pelargonidin + diet.chocolate$Petunidin + diet.chocolate$Peonidin

# flavanones 
diet.chocolate$flavanones.chocolate <- diet.chocolate$Eriodictyol + diet.chocolate$Hesperetin + diet.chocolate$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.chocolate$flavan.3.ols.chocolate <- diet.chocolate$X....Catechin + diet.chocolate$X....Epicatechin + diet.chocolate$X....Epigallocatechin + diet.chocolate$X....Gallocatechin + diet.chocolate$X....Catechin.3.O.gallate + diet.chocolate$X....Epicatechin.3.O.gallate + diet.chocolate$X....Gallocatechin.3.O.gallate +  diet.chocolate$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.chocolate$flavonols.chocolate <- diet.chocolate$Quercetin + diet.chocolate$Kaempferol + diet.chocolate$Myricetin + diet.chocolate$Isorhamnetin

# flavones
diet.chocolate$flavones.chocolate <- diet.chocolate$Luteolin + diet.chocolate$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.chocolate$polymers.chocolate <-  diet.chocolate$X02.mers + diet.chocolate$X03.mers + diet.chocolate$X04.06.mers + diet.chocolate$X07.10.mers + diet.chocolate$Polymers...10.mers.

# Proanthocyanidins
diet.chocolate$proanthocyanidins.chocolate <- diet.chocolate$Theaflavin + diet.chocolate$Thearubigins

# total flavonoid
diet.chocolate$total.flavonoid.chocolate <- diet.chocolate$anthocyanins.chocolate + diet.chocolate$flavanones.chocolate + diet.chocolate$flavan.3.ols.chocolate + diet.chocolate$flavonols.chocolate + diet.chocolate$flavones + diet.chocolate$polymers.chocolate + diet.chocolate$proanthocyanidins.chocolate

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.chocolate$group <- factor(diet.chocolate$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.chocolate)[96:103]  # Selecting the variable names from column 3 to 217
description.chocolate <- CreateTableOne(data = diet.chocolate,
                                        vars = vars,  # Variable names extracted above
                                        strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.chocolate, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.chocolate

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.chocolate), "flavonoid/chocolateresult.csv")

# ---- Rmd chunk ----
rm(diet.chocolate)
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.alcohol <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.alcohol[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.alcohol) <- c("Participant.ID", "group")

# Define a function to sum the alcohol values, adjust with p.24h.foods.average, and add to diet.alcohol dataframe
sum.alcohol.values <- function(alcohol.rows, flavonoid.data, diet.alcohol) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.alcohol (each participant)
  for (participant in seq_len(nrow(diet.alcohol))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified alcohol row and add its values to the summed.values vector
    for (i in seq_along(alcohol.rows)) {
      alcohol.row <- alcohol.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (alcohol.row == 58) {
        # Red.wine.intake (row 58)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[alcohol.row, ]) * p.24h.foods.average$Red.wine.intake[participant]) * 165.76 / 100
      } else if (alcohol.row == 59) {
        # Rose.wine.intake (row 59)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[alcohol.row, ]) * p.24h.foods.average$Rose.wine.intake[participant]) * 163.50 / 100
      } else if (alcohol.row == 60) {
        # White.wine.intake (row 60)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[alcohol.row, ]) * p.24h.foods.average$White.wine.intake[participant]) * 162.96 / 100
      } else if (alcohol.row == 61) {
        # Fortified.wine.intake (row 61)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[alcohol.row, ]) * p.24h.foods.average$Fortified.wine.intake[participant]) * 50 / 100
      } else if (alcohol.row == 68) {
        # Beer.cider.intake (row 68)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[alcohol.row, ]) * p.24h.foods.average$Beer.cider.intake[participant]) * 568 / 100
      } 
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.alcohol dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.alcohol[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.alcohol)
}

# Define the alcohol rows (adjust based on the actual rows you want to sum)
alcohol.rows <- c(58, 59, 60, 61, 68)

# Call the function to calculate and update the diet.alcohol dataframe, incorporating p.24h.foods.average adjustments
diet.alcohol <- sum.alcohol.values(alcohol.rows, flavonoid.checkout, diet.alcohol)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.alcohol, "flavonoid/diet.alcohol.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
diet.alcohol$anthocyanins.alcohol <- diet.alcohol$Cyanidin + diet.alcohol$Delphinidin + diet.alcohol$Malvidin + diet.alcohol$Pelargonidin + diet.alcohol$Petunidin + diet.alcohol$Peonidin

# flavanones 
diet.alcohol$flavanones.alcohol <- diet.alcohol$Eriodictyol + diet.alcohol$Hesperetin + diet.alcohol$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.alcohol$flavan.3.ols.alcohol <- diet.alcohol$X....Catechin + diet.alcohol$X....Epicatechin + diet.alcohol$X....Epigallocatechin + diet.alcohol$X....Gallocatechin + diet.alcohol$X....Catechin.3.O.gallate + diet.alcohol$X....Epicatechin.3.O.gallate + diet.alcohol$X....Gallocatechin.3.O.gallate +  diet.alcohol$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.alcohol$flavonols.alcohol <- diet.alcohol$Quercetin + diet.alcohol$Kaempferol + diet.alcohol$Myricetin + diet.alcohol$Isorhamnetin

# flavones
diet.alcohol$flavones.alcohol <- diet.alcohol$Luteolin + diet.alcohol$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.alcohol$polymers.alcohol <-  diet.alcohol$X02.mers + diet.alcohol$X03.mers + diet.alcohol$X04.06.mers + diet.alcohol$X07.10.mers + diet.alcohol$Polymers...10.mers.

# Proanthocyanidins
diet.alcohol$proanthocyanidins.alcohol <- diet.alcohol$Theaflavin + diet.alcohol$Thearubigins

# total flavonoid
diet.alcohol$total.flavonoid.alcohol <- diet.alcohol$anthocyanins.alcohol + diet.alcohol$flavanones.alcohol + diet.alcohol$flavan.3.ols.alcohol + diet.alcohol$flavonols.alcohol + diet.alcohol$flavones + diet.alcohol$polymers.alcohol + diet.alcohol$proanthocyanidins.alcohol

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.alcohol$group <- factor(diet.alcohol$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.alcohol)[96:103]  # Selecting the variable names from column 3 to 217
description.alcohol <- CreateTableOne(data = diet.alcohol,
                                      vars = vars,  # Variable names extracted above
                                      strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.alcohol, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.alcohol

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.alcohol), "flavonoid/alcoholresult.ap.csv")

# ---- Rmd chunk ----
rm(diet.alcohol)
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.juice <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.juice[, 1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.juice) <- c("Participant.ID", "group")

# Define a function to sum the juice values, adjust with p.24h.foods.average, and add to diet.juice dataframe
sum.juice.values <- function(juice.rows, flavonoid.data, diet.juice) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.juice (each participant)
  for (participant in seq_len(nrow(diet.juice))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified juice row and add its values to the summed.values vector
    for (i in seq_along(juice.rows)) {
      juice.row <- juice.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (juice.row == 85) {
        # Orange.juice.intake (row 85)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[juice.row, ]) * p.24h.foods.average$Orange.juice.intake[participant]) * 250 / 100
      } else if (juice.row == 78) {
        # Grapefruit.juice.intake (row 78)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[juice.row, ]) * p.24h.foods.average$Grapefruit.juice.intake[participant]) * 250 / 100
      } else if (juice.row == 90) {
        # Pure.fruit.vegetable.juice.intake (row 90)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[juice.row, ]) * p.24h.foods.average$Pure.fruit.vegetable.juice.intake[participant]) * 250 / 100
      } else if (juice.row == 88) {
        # Fruit.smoothie.intake (row 88)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[juice.row, ]) * p.24h.foods.average$Fruit.smoothie.intake[participant]) * 250 / 100
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.juice dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.juice[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.juice)
}

# Define the juice rows (adjust based on the actual rows you want to sum)
juice.rows <- c(85, 78, 90, 88)

# Call the function to calculate and update the diet.juice dataframe, incorporating p.24h.foods.average adjustments
diet.juice <- sum.juice.values(juice.rows, flavonoid.checkout, diet.juice)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.juice, "flavonoid/diet.juice.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
diet.juice$anthocyanins.juice <- diet.juice$Cyanidin + diet.juice$Delphinidin + diet.juice$Malvidin + diet.juice$Pelargonidin + diet.juice$Petunidin + diet.juice$Peonidin

# flavanones 
diet.juice$flavanones.juice <- diet.juice$Eriodictyol + diet.juice$Hesperetin + diet.juice$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.juice$flavan.3.ols.juice <- diet.juice$X....Catechin + diet.juice$X....Epicatechin + diet.juice$X....Epigallocatechin + diet.juice$X....Gallocatechin + diet.juice$X....Catechin.3.O.gallate + diet.juice$X....Epicatechin.3.O.gallate + diet.juice$X....Gallocatechin.3.O.gallate +  diet.juice$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.juice$flavonols.juice <- diet.juice$Quercetin + diet.juice$Kaempferol + diet.juice$Myricetin + diet.juice$Isorhamnetin

# flavones
diet.juice$flavones.juice <- diet.juice$Luteolin + diet.juice$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.juice$polymers.juice <-  diet.juice$X02.mers + diet.juice$X03.mers + diet.juice$X04.06.mers + diet.juice$X07.10.mers + diet.juice$Polymers...10.mers.

# Proanthocyanidins
diet.juice$proanthocyanidins.juice <- diet.juice$Theaflavin + diet.juice$Thearubigins

# total flavonoid
diet.juice$total.flavonoid.juice <- diet.juice$anthocyanins.juice + diet.juice$flavanones.juice + diet.juice$flavan.3.ols.juice + diet.juice$flavonols.juice + diet.juice$flavones + diet.juice$polymers.juice + diet.juice$proanthocyanidins.juice

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.juice$group <- factor(diet.juice$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.juice)[96:103]  # Selecting the variable names from column 3 to 217
description.juice <- CreateTableOne(data = diet.juice,
                                    vars = vars,  # Variable names extracted above
                                    strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.juice, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.juice

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.juice), "flavonoid/juiceresult.csv")

# ---- Rmd chunk ----
rm(diet.juice)
appendixfood <- read.csv("24hr/appendixtofoods.csv") %>%
  filter(Participant.ID %in% demographics$Participant.ID)
# Replace "half" with 0.5, "6+" with 6, and convert to numeric, replacing NAs with 0
appendixfood <- appendixfood %>%
  mutate(across(contains("intake"), ~ gsub("half", "0.5", .))) %>%
  mutate(across(contains("intake"), ~ gsub("6\\+", "6", .))) %>%
  mutate(across(contains("intake"), ~ ifelse(is.na(.), "0", .))) %>%   # Replace NA with "0"
  mutate(across(contains("intake"), as.numeric))  # Convert to numeric

# Get all column names for intake related to bread types
sliced_intake_columns <- paste0("Sliced.bread.intake...Instance.", 0:4)
baguette_intake_columns <- paste0("Baguette.intake...Instance.", 0:4)
bap_intake_columns <- paste0("Bap.intake...Instance.", 0:4)
bread_roll_intake_columns <- paste0("Bread.roll.intake...Instance.", 0:4)

# Define the correct type columns for each bread type
sliced_type_columns <- paste0("Type.of.sliced.bread.eaten...Instance.", 0:4)
baguette_type_columns <- paste0("Type.of.baguette.eaten...Instance.", 0:4)  # Corrected this
bap_type_columns <- paste0("Type.of.large.bap.eaten...Instance.", 0:4)
bread_roll_type_columns <- paste0("Type.of.bread.roll.eaten...Instance.", 0:4)


# Replace all NA values in intake columns with 0
appendixfood <- appendixfood %>%
  mutate(across(contains("intake"), ~ ifelse(is.na(.), 0, .)))

# Function to calculate mean intake by type (wholemeal, white, mixed, seeded, others)
calc_mean_by_type <- function(data_row, type, type_columns, intake_columns) {
  # For each column, check if the type is present and return the intake value, otherwise return 0
  intake_values <- mapply(function(type_col, intake_col) {
    ifelse(grepl(type, data_row[[type_col]]), data_row[[intake_col]], 0)
  }, type_columns, intake_columns)
  
  # Calculate the mean of the intake values
  mean(intake_values)
}

appendixfood <- appendixfood %>%
  rowwise() %>%
  
  # Sliced Bread
  mutate(sliced_wholemeal.mean = calc_mean_by_type(cur_data(), "wholemeal", sliced_type_columns, sliced_intake_columns),
         sliced_white.mean = calc_mean_by_type(cur_data(), "white", sliced_type_columns, sliced_intake_columns),
         sliced_mixed.mean = calc_mean_by_type(cur_data(), "mixed", sliced_type_columns, sliced_intake_columns),
         sliced_seeded.mean = calc_mean_by_type(cur_data(), "seeded", sliced_type_columns, sliced_intake_columns),
         sliced_others.mean = calc_mean_by_type(cur_data(), "others", sliced_type_columns, sliced_intake_columns)) %>%
  
  # Baguette
  mutate(baguette_wholemeal.mean = calc_mean_by_type(cur_data(), "wholemeal", baguette_type_columns, baguette_intake_columns),
         baguette_white.mean = calc_mean_by_type(cur_data(), "white", baguette_type_columns, baguette_intake_columns),
         baguette_mixed.mean = calc_mean_by_type(cur_data(), "mixed", baguette_type_columns, baguette_intake_columns),
         baguette_seeded.mean = calc_mean_by_type(cur_data(), "seeded", baguette_type_columns, baguette_intake_columns),
         baguette_others.mean = calc_mean_by_type(cur_data(), "others", baguette_type_columns, baguette_intake_columns)) %>%
  
  # Bap
  mutate(bap_wholemeal.mean = calc_mean_by_type(cur_data(), "wholemeal", bap_type_columns, bap_intake_columns),
         bap_white.mean = calc_mean_by_type(cur_data(), "white", bap_type_columns, bap_intake_columns),
         bap_mixed.mean = calc_mean_by_type(cur_data(), "mixed", bap_type_columns, bap_intake_columns),
         bap_seeded.mean = calc_mean_by_type(cur_data(), "seeded", bap_type_columns, bap_intake_columns),
         bap_others.mean = calc_mean_by_type(cur_data(), "others", bap_type_columns, bap_intake_columns)) %>%
  
  # Bread Roll
  mutate(bread_roll_wholemeal.mean = calc_mean_by_type(cur_data(), "wholemeal", bread_roll_type_columns, bread_roll_intake_columns),
         bread_roll_white.mean = calc_mean_by_type(cur_data(), "white", bread_roll_type_columns, bread_roll_intake_columns),
         bread_roll_mixed.mean = calc_mean_by_type(cur_data(), "mixed", bread_roll_type_columns, bread_roll_intake_columns),
         bread_roll_seeded.mean = calc_mean_by_type(cur_data(), "seeded", bread_roll_type_columns, bread_roll_intake_columns),
         bread_roll_others.mean = calc_mean_by_type(cur_data(), "others", bread_roll_type_columns, bread_roll_intake_columns)) %>%
  
  ungroup()

# Save the bread data
appendix <- appendixfood[,c(1, 87:106)]
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(appendix, "flavonoid/processed bread data.csv", row.names = FALSE)

# ---- Rmd chunk ----
appendix <- read.csv("flavonoid/processed bread data.csv")
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.cereal <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.cereal[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.cereal) <- c("Participant.ID", "group")

# Define a function to sum the cereal values, adjust with p.24h.foods.average, and add to diet.cereal dataframe
sum.cereal.values <- function(cereal.rows, flavonoid.data, diet.cereal) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.cereal (each participant)
  for (participant in seq_len(nrow(diet.cereal))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified cereal row and add its values to the summed.values vector
    for (i in seq_along(cereal.rows)) {
      cereal.row <- cereal.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (cereal.row == 106) {
        # Porridge intake and Oat crunch intake (row 106)
        summed.values <- summed.values + 
          (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Porridge.intake[participant]) *40 /100 +
          (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Oat.crunch.intake[participant]) * 30 / 100
      } else if (cereal.row == 92) {
        # Muesli.intake (row 92)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Muesli.intake[participant]) * 30 / 100
      } else if (cereal.row == 93) {
        # Bran.cereal.intake (row 93)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Bran.cereal.intake[participant]) * 30 / 100
      } else if (cereal.row == 104) {
        # Whole.wheat.cereal.intake (row 104)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Whole.wheat.cereal.intake[participant]) * 30 / 100
      } else if (cereal.row == 120) {
        # Sliced.bread, baguette, bap, bread roll -- wholemeal (row 120)
        summed.values <- summed.values + ((as.numeric(flavonoid.columns[cereal.row, ]) * appendix$sliced_wholemeal.mean[participant]) * 35 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$baguette_wholemeal.mean[participant]) * 140 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bap_wholemeal.mean[participant]) * 45 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bread_roll_wholemeal.mean[participant]) * 60 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$sliced_mixed.mean[participant]) * 35 / 100 *0.5 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$baguette_mixed.mean[participant]) * 140 / 100 *0.5 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bap_mixed.mean[participant]) * 45 / 100 *0.5 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bread_roll_mixed.mean[participant]) * 60 / 100 *0.5)* 48/64*0.9 
      } else if (cereal.row == 202) {
        # Sliced.bread, baguette, bap, bread roll -- seeded (row 202)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$sliced_seeded.mean[participant]) * 35 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$baguette_seeded.mean[participant]) * 140 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bap_seeded.mean[participant]) * 45 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bread_roll_seeded.mean[participant]) * 60 / 100
      } else if (cereal.row == 167) {
        # Sliced.bread, baguette, bap, bread roll -- others (row 167)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$sliced_others.mean[participant]) * 35 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$baguette_others.mean[participant]) * 140 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bap_others.mean[participant]) * 45 / 100 + (as.numeric(flavonoid.columns[cereal.row, ]) * appendix$bread_roll_others.mean[participant]) * 60 / 100
      } else if (cereal.row == 884) {
        # Naan bread intake (row 884)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Naan.bread.intake[participant]) * 140 / 100 
      } else if (cereal.row == 164) {
        # Garlic bread intake (row 164)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Garlic.bread.intake[participant]) * 35 / 100 
      } else if (cereal.row == 166) {
        # Crispbread intake (row 166)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Crispbread.intake[participant]) * 8 / 100 
      } else if (cereal.row == 134) {
        # Oatcakes intake (row 134)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Oatcakes.intake[participant]) * 12 / 100 
      } else if (cereal.row == 153) {
        # Other bread intake -- tortilla wraps,corn (row 153)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Other.bread.intake[participant]) * 60 / 100
      } else if (cereal.row == 217) {
        # White pasta intake (row 217)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$White.pasta.intake[participant]) * 180 / 100
      } else if (cereal.row == 135) {
        # Wholemeal pasta intake (row 135)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Wholemeal.pasta.intake[participant]) * 180 / 100
      } else if (cereal.row == 196) {
        # Brown rice intake (row 196)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Brown.rice.intake[participant]) * 180 / 100
      } else if (cereal.row == 180) {
        # Snackpot intake (row 180)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Snackpot.intake[participant]) * 65 / 100
      } else if (cereal.row == 128) {
        # Couscous intake (row 128)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Couscous.intake[participant]) * 150 / 100
      } else if (cereal.row == 203) {
        # Other grain intake (row 203)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[cereal.row, ]) * p.24h.foods.average$Other.grain.intake[participant]) * 180 / 100
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.cereal dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.cereal[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.cereal)
}

# Define the cereal rows (adjust based on the actual rows you want to sum)
cereal.rows <- c(106, 92, 93, 104, 120, 202, 167, 884, 164, 166, 133, 153, 217, 135, 196, 180, 128)

# Call the function to calculate and update the diet.cereal dataframe, incorporating p.24h.foods.average adjustments
diet.cereal <- sum.cereal.values(cereal.rows, flavonoid.checkout, diet.cereal)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.cereal, "flavonoid/diet.cereal.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
diet.cereal$anthocyanins.cereal <- diet.cereal$Cyanidin + diet.cereal$Delphinidin + diet.cereal$Malvidin + diet.cereal$Pelargonidin + diet.cereal$Petunidin + diet.cereal$Peonidin

# flavanones 
diet.cereal$flavanones.cereal <- diet.cereal$Eriodictyol + diet.cereal$Hesperetin + diet.cereal$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.cereal$flavan.3.ols.cereal <- diet.cereal$X....Catechin + diet.cereal$X....Epicatechin + diet.cereal$X....Epigallocatechin + diet.cereal$X....Gallocatechin + diet.cereal$X....Catechin.3.O.gallate + diet.cereal$X....Epicatechin.3.O.gallate + diet.cereal$X....Gallocatechin.3.O.gallate +  diet.cereal$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.cereal$flavonols.cereal <- diet.cereal$Quercetin + diet.cereal$Kaempferol + diet.cereal$Myricetin + diet.cereal$Isorhamnetin

# flavones
diet.cereal$flavones.cereal <- diet.cereal$Luteolin + diet.cereal$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.cereal$polymers.cereal <-  diet.cereal$X02.mers + diet.cereal$X03.mers + diet.cereal$X04.06.mers + diet.cereal$X07.10.mers + diet.cereal$Polymers...10.mers.

# Proanthocyanidins
diet.cereal$proanthocyanidins.cereal <- diet.cereal$Theaflavin + diet.cereal$Thearubigins

# total flavonoid
diet.cereal$total.flavonoid.cereal <- diet.cereal$anthocyanins.cereal + diet.cereal$flavanones.cereal + diet.cereal$flavan.3.ols.cereal + diet.cereal$flavonols.cereal + diet.cereal$flavones + diet.cereal$polymers.cereal + diet.cereal$proanthocyanidins.cereal

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.cereal$group <- factor(diet.cereal$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.cereal)[96:103]  # Selecting the variable names from column 3 to 217
description.cereal <- CreateTableOne(data = diet.cereal,
                                     vars = vars,  # Variable names extracted above
                                     strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.cereal, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.cereal


# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.cereal), "flavonoid/cerealresult.csv")

# ---- Rmd chunk ----
rm(diet.vegan)
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.fruit <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.fruit[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.fruit) <- c("Participant.ID", "group")

# Define a function to sum the fruit values, adjust with p.24h.foods.average, and add to diet.fruit dataframe
sum.fruit.values <- function(fruit.rows, flavonoid.data, diet.fruit) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.fruit (each participant)
  for (participant in seq_len(nrow(diet.fruit))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified fruit row and add its values to the summed.values vector
    for (i in seq_along(fruit.rows)) {
      fruit.row <- fruit.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (fruit.row == 348) {
        # Stewed fruit intake - Apples (row 348)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Stewed.fruit.intake[participant]) * 80 / 100
      } else if (fruit.row == 298) {
        # Prune intake (row 298)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Prune.intake[participant]) * 30 / 100
      } else if (fruit.row == 366) {
        # Dried fruit intake (row 366)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Dried.fruit.intake[participant]) * 30 / 100
      } else if (fruit.row == 304) {
        # Mixed fruit intake (row 304)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Mixed.fruit.intake[participant]) * 80 / 100
      } else if (fruit.row == 252) {
        # Apple intake (row 252)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Apple.intake[participant]) * 80 / 100
      } else if (fruit.row == 300) {
        # Banana intake (row 300)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Banana.intake[participant]) * 80 / 100
      } else if (fruit.row == 354) {
        # berry intake (row 354)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Berry.intake[participant]) * 80 / 100
      } else if (fruit.row == 346) {
        # Cherry intake (row 346)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Cherry.intake[participant]) * 80 / 100
      } else if (fruit.row == 369) {
        # Grapefruit intake (row 369)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Grapefruit.intake[participant]) * 80 / 100
      } else if (fruit.row == 377) {
        # Grape intake - black (row 377)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Grape.intake[participant]) * 80 / 100 *0.32
      } else if (fruit.row == 378) {
        # Grape intake - green (row 378)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Grape.intake[participant]) * 80 / 100 *0.68
      } else if (fruit.row == 307) {
        # Mango intake (row 307)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Mango.intake[participant]) * 80 / 100
      } else if (fruit.row == 311) {
        # Melon intake (row 311)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Melon.intake[participant]) * 80 / 100
      } else if (fruit.row == 372) {
        # Orange intake (row 372)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Orange.intake[participant]) * 80 / 100
      } else if (fruit.row == 373) {
        # Tangerine intake (row 373)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Satsuma.intake[participant]) * 80 / 100
      } else if (fruit.row == 340) {
        # Peach intake (row 340)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Peach.nectarine.intake[participant]) * 80 / 100
      } else if (fruit.row == 255) {
        # Pear intake (row 255)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Pear.intake[participant]) * 80 / 100
      } else if (fruit.row == 319) {
        # Pineapple intake (row 319)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Pineapple.intake[participant]) * 80 / 100
      } else if (fruit.row == 344) {
        # Plum intake (row 344)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Plum.intake[participant]) * 40 / 100
      } else if (fruit.row == 321) {
        # Other fruit intake - pomegranate (row 321)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Other.fruit.intake[participant]) * 80 / 100 *0.5
      } else if (fruit.row == 305) {
        # Other fruit intake - kiwi (row 305)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[fruit.row, ]) * p.24h.foods.average$Other.fruit.intake[participant]) * 80 / 100 *0.5
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.fruit dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.fruit[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.fruit)
}

# Define the fruit rows (adjust based on the actual rows you want to sum)
fruit.rows <- c(348, 298, 366, 304, 252, 300, 354, 346, 369, 377, 378, 307, 311, 372, 373, 340, 255, 319, 344, 321, 305)

# Call the function to calculate and update the diet.fruit dataframe, incorporating p.24h.foods.average adjustments
diet.fruit <- sum.fruit.values(fruit.rows, flavonoid.checkout, diet.fruit)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.fruit, "flavonoid/diet.fruit.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
diet.fruit$anthocyanins.fruit <- diet.fruit$Cyanidin + diet.fruit$Delphinidin + diet.fruit$Malvidin + diet.fruit$Pelargonidin + diet.fruit$Petunidin + diet.fruit$Peonidin

# flavanones 
diet.fruit$flavanones.fruit <- diet.fruit$Eriodictyol + diet.fruit$Hesperetin + diet.fruit$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.fruit$flavan.3.ols.fruit <- diet.fruit$X....Catechin + diet.fruit$X....Epicatechin + diet.fruit$X....Epigallocatechin + diet.fruit$X....Gallocatechin + diet.fruit$X....Catechin.3.O.gallate + diet.fruit$X....Epicatechin.3.O.gallate + diet.fruit$X....Gallocatechin.3.O.gallate +  diet.fruit$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.fruit$flavonols.fruit <- diet.fruit$Quercetin + diet.fruit$Kaempferol + diet.fruit$Myricetin + diet.fruit$Isorhamnetin

# flavones
diet.fruit$flavones.fruit <- diet.fruit$Luteolin + diet.fruit$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.fruit$polymers.fruit <-  diet.fruit$X02.mers + diet.fruit$X03.mers + diet.fruit$X04.06.mers + diet.fruit$X07.10.mers + diet.fruit$Polymers...10.mers.

# Proanthocyanidins
diet.fruit$proanthocyanidins.fruit <- diet.fruit$Theaflavin + diet.fruit$Thearubigins

# total flavonoid
diet.fruit$total.flavonoid.fruit <- diet.fruit$anthocyanins.fruit + diet.fruit$flavanones.fruit + diet.fruit$flavan.3.ols.fruit + diet.fruit$flavonols.fruit + diet.fruit$flavones + diet.fruit$polymers.fruit + diet.fruit$proanthocyanidins.fruit

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.fruit$group <- factor(diet.fruit$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.fruit)[96:103]  # Selecting the variable names from column 3 to 217
description.fruit <- CreateTableOne(data = diet.fruit,
                                    vars = vars,  # Variable names extracted above
                                    strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.fruit, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.fruit

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.fruit), "flavonoid/fruitresult.csv")

# ---- Rmd chunk ----
rm(diet.fruit)
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.vegetable <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.vegetable[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.vegetable) <- c("Participant.ID", "group")

# Define a function to sum the vegetable values, adjust with p.24h.foods.average, and add to diet.vegetable dataframe
sum.vegetable.values <- function(vegetable.rows, flavonoid.data, diet.vegetable) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.vegetable (each participant)
  for (participant in seq_len(nrow(diet.vegetable))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified vegetable row and add its values to the summed.values vector
    for (i in seq_along(vegetable.rows)) {
      vegetable.row <- vegetable.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (vegetable.row == 578) {
        # Fried potatoes intake (row 578)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Fried.potatoes.intake[participant]) * 180 / 100
      } else if (vegetable.row == 572) {
        # Boiled/baked potatoes intake (row 572)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Boiled.baked.potatoes.intake[participant]) * 180 / 100
      } else if (vegetable.row == 576) {
        # Mashed potatoes intake (row 576)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Mashed.potato.intake[participant]) * 180 / 100
      } else if (vegetable.row == 614) {
        # Mixed vegetable intake (row 614)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Mixed.vegetable.intake[participant]) * 80 / 100
      } else if (vegetable.row == 615) {
        # Vegetable pieces intake (row 615)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Vegetable.pieces.intake[participant]) * 80 / 100
      } else if (vegetable.row == 623) {
        # Coleslaw intake (row 623)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Coleslaw.intake[participant]) * 80 / 100
      } else if (vegetable.row == 667) {
        # Side salad intake (row 667)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Side.salad.intake[participant]) * 80 / 100
      } else if (vegetable.row == 379) {
        # Avocado intake (row 379)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Avocado.intake[participant]) * 160 / 100
      } else if (vegetable.row == 610) {
        # Green bean intake (row 610)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Green.bean.intake[participant]) * 80 / 100 * 0.93
      } else if (vegetable.row == 587) {
        # Beetroot intake (row 587)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Beetroot.intake[participant]) * 80 / 100
      } else if (vegetable.row == 505) {
        # Broccoli intake (row 505)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Broccoli.intake[participant]) * 80 / 100
      } else if (vegetable.row == 565) {
        # Butternut squash intake (row 565)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Butternut.squash.intake[participant]) * 80 / 100
      } else if (vegetable.row == 503) {
        # Cabbage intake (row 503)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Cabbage.kale.intake[participant]) * 80 / 100 *0.5
      } else if (vegetable.row == 500) {
        # Kale intake (row 500)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Cabbage.kale.intake[participant]) * 80 / 100 *0.5
      } else if (vegetable.row == 592) {
        # Carrot intake (row 592)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Carrot.intake[participant]) * 80 / 100
      } else if (vegetable.row == 499) {
        # Cauliflower intake (row 499)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Cauliflower.intake[participant]) * 80 / 100
      } else if (vegetable.row == 465) {
        # Celery intake (row 465)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Celery.intake[participant]) * 26.7 / 100
      } else if (vegetable.row == 559) {
        # Courgette intake (row 559)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Courgette.intake[participant]) * 80 / 100
      } else if (vegetable.row == 467) {
        # Cucumber intake (row 467)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Cucumber.intake[participant]) * 32 / 100
      } else if (vegetable.row == 419) {
        # Garlic intake (row 419)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Garlic.intake[participant]) * 5 / 100
      } else if (vegetable.row == 522) {
        # Leek intake (row 522)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Leek.intake[participant]) * 80 / 100
      } else if (vegetable.row == 702) {
        # Lettuce intake (row 702)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Lettuce.intake[participant]) * 80 / 100
      } else if (vegetable.row == 712) {
        # Mushroom intake (row 712)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Mushroom.intake[participant]) * 80 / 100
      } else if (vegetable.row == 521) {
        # Onion intake (row 521)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Onion.intake[participant]) * 80 / 100
      } else if (vegetable.row == 586) {
        # Parsnip intake (row 586)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Parsnip.intake[participant]) * 80 / 100
      } else if (vegetable.row == 397) {
        # Sweet pepper intake (row 397)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Sweet.pepper.intake[participant]) * 80 / 100
      } else if (vegetable.row == 513) {
        # Spinach intake - cooked (row 513)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Spinach.intake[participant]) * 80 / 100 *0.5
      } else if (vegetable.row == 706) {
        # Spinach intake - raw (row 706)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Spinach.intake[participant]) * 80 / 100 *0.5
      } else if (vegetable.row == 507) {
        # Sprouts intake (row 507)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Sprouts.intake[participant]) * 80 / 100
      } else if (vegetable.row == 139) {
        # Sweetcorn intake (row 139)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Sweetcorn.intake[participant]) * 80 / 100
      } else if (vegetable.row == 571) {
        # Sweet potato intake (row 571)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Sweet.potato.intake[participant]) * 130 / 100
      } else if (vegetable.row == 392) {
        # Fresh tomato intake (row 392)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Fresh.tomato.intake[participant]) * 80 / 100
      } else if (vegetable.row == 390) {
        # Tinned tomato intake (row 390)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Tinned.tomato.intake[participant]) * 200 / 100
      } else if (vegetable.row == 597) {
        # Turnip intake (row 597)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Turnip.swede.intake[participant]) * 80 / 100 *0.5
      } else if (vegetable.row == 589) {
        # Swede intake (row 589)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Turnip.swede.intake[participant]) * 80 / 100 *0.5
      } else if (vegetable.row == 695) {
        # Watercress intake (row 695)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Turnip.swede.intake[participant]) * 80 / 100
      } else if (vegetable.row == 564) {
        # Other vegetables intake - asparagus (row 564)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Other.vegetables.intake[participant]) * 80 / 100 *0.4
      } else if (vegetable.row == 394) {
        # Other vegetables intake - aubergine (row 394)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Other.vegetables.intake[participant]) * 80 / 100 *0.3
      } else if (vegetable.row == 566) {
        # Other vegetables intake - squash (row 565)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Other.vegetables.intake[participant]) * 80 / 100 *0.3
      } else if (vegetable.row == 847) {
        # Olives Intake - Black (row 847)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Olives.intake[participant]) * 50 / 100 *0.5
      } else if (vegetable.row == 848) {
        # Olives Intake - Green (row 848)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[vegetable.row, ]) * p.24h.foods.average$Olives.intake[participant]) * 50 / 100 *0.5
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.vegetable dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.vegetable[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.vegetable)
}

# Define the vegetable rows (adjust based on the actual rows you want to sum)
vegetable.rows <- c(578, 572, 576, 614, 615, 623, 667, 379, 610, 587, 505, 565, 503, 500, 592, 499, 465, 559, 467, 419, 522, 702, 712, 521, 586, 397, 513, 706, 507, 139, 571, 392, 390, 597, 589, 695, 564, 394, 566, 847, 848)

# Call the function to calculate and update the diet.vegetable dataframe, incorporating p.24h.foods.average adjustments
diet.vegetable <- sum.vegetable.values(vegetable.rows, flavonoid.checkout, diet.vegetable)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.vegetable, "flavonoid/diet.vegetable.csv", row.names = FALSE)


# ---- Rmd chunk ----
# anthocyanins
diet.vegetable$anthocyanins.vegetable <- diet.vegetable$Cyanidin + diet.vegetable$Delphinidin + diet.vegetable$Malvidin + diet.vegetable$Pelargonidin + diet.vegetable$Petunidin + diet.vegetable$Peonidin

# flavanones 
diet.vegetable$flavanones.vegetable <- diet.vegetable$Eriodictyol + diet.vegetable$Hesperetin + diet.vegetable$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.vegetable$flavan.3.ols.vegetable <- diet.vegetable$X....Catechin + diet.vegetable$X....Epicatechin + diet.vegetable$X....Epigallocatechin + diet.vegetable$X....Gallocatechin + diet.vegetable$X....Catechin.3.O.gallate + diet.vegetable$X....Epicatechin.3.O.gallate + diet.vegetable$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.vegetable$flavonols.vegetable <- diet.vegetable$Quercetin + diet.vegetable$Kaempferol + diet.vegetable$Myricetin + diet.vegetable$Isorhamnetin

# flavones
diet.vegetable$flavones.vegetable <- diet.vegetable$Luteolin + diet.vegetable$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.vegetable$polymers.vegetable <- diet.vegetable$Theaflavin + diet.vegetable$Thearubigins + diet.vegetable$X02.mers + diet.vegetable$X03.mers + diet.vegetable$X04.06.mers + diet.vegetable$X07.10.mers + diet.vegetable$Polymers...10.mers.

# Proanthocyanidins
diet.vegetable$proanthocyanidins.vegetable <- diet.vegetable$X....Catechin + diet.vegetable$X....Epicatechin + diet.vegetable$X....Epigallocatechin + diet.vegetable$X....Gallocatechin + diet.vegetable$X....Catechin.3.O.gallate + diet.vegetable$X....Epicatechin.3.O.gallate + diet.vegetable$X....Epigallocatechin.3.O.gallate + diet.vegetable$X02.mers + diet.vegetable$X03.mers + diet.vegetable$X04.06.mers + diet.vegetable$X07.10.mers + diet.vegetable$Polymers...10.mers.

diet.vegetable$total.flavonoid.vegetable <- diet.vegetable$anthocyanins + diet.vegetable$flavanones + diet.vegetable$flavan.3.ols + diet.vegetable$flavonols + diet.vegetable$flavones + diet.vegetable$polymers

diet.vegetable <- merge(diet.vegetable,dietaryscorewithquartile[, c(1,10)], by = "Participant.ID")# Set the 'group'
# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.vegetable$group <- factor(diet.vegetable$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the vegetable variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.vegetable)[96:103]  # Selecting the variable names from column 3 to 217
description.vegetable <- CreateTableOne(data = diet.vegetable,
                                        vars = vars,  # Variable names extracted above
                                        strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.vegetable, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.vegetable

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.vegetable), "flavonoid/vegetableresult.csv")

# ---- Rmd chunk ----
rm(diet.snack)
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.legume <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.legume[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.legume) <- c("Participant.ID", "group")

# Define a function to sum the legume values, adjust with p.24h.foods.average, and add to diet.legume dataframe
sum.legume.values <- function(legume.rows, flavonoid.data, diet.legume) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.legume (each participant)
  for (participant in seq_len(nrow(diet.legume))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified legume row and add its values to the summed.values vector
    for (i in seq_along(legume.rows)) {
      legume.row <- legume.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (legume.row == 433) {
        # Pulse intake -Beans, kidney, all types, mature seeds, cooked, boiled, without salt (row 433)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Pulses.intake[participant]) * 150 / 100 *0.20
      } else if (legume.row == 439) {
        # Pulse intake -Chickpeas (garbanzo beans, bengal gram), mature seeds, cooked, boiled, without salt (row 439)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Pulses.intake[participant]) * 150 / 100 *0.30
      } else if (legume.row == 458) {
        # Pulse intake -Lima beans, Butter beans, mature seeds, cooked, boiled, drained, without salt (row 458)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Pulses.intake[participant]) * 150 / 100 *0.10
      } else if (legume.row == 527) {
        # Pulse intake -Lentils, red, split, dried, boiled in unsalted water (row 527)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Pulses.intake[participant]) * 150 / 100 *0.25
      } else if (legume.row == 544) {
        # Pulse intake -Cannellini beans, dried, cooked in unsalted water (row 544)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Pulses.intake[participant]) * 150 / 100 *0.10
      } else if (legume.row == 541) {
        # Pulse intake -Beans, pinto, dried, boiled in unsalted water (row 541)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Pulses.intake[participant]) * 150 / 100 *0.05
      } else if (legume.row == 626) {
        # Broad bean intake (row 626)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Broad.bean.intake[participant]) * 150 / 100
      } else if (legume.row == 529) {
        # Pea intake (row 529)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[legume.row, ]) * p.24h.foods.average$Pea.intake[participant]) * 80 / 100
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.legume dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.legume[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.legume)
}

# Define the legume rows (adjust based on the actual rows you want to sum)
legume.rows <- c(433, 439, 458, 527, 544, 541, 626, 529)

# Call the function to calculate and update the diet.legume dataframe, incorporating p.24h.foods.average adjustments
diet.legume <- sum.legume.values(legume.rows, flavonoid.checkout, diet.legume)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.legume, "flavonoid/diet.legume.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
diet.legume$anthocyanins.legume <- diet.legume$Cyanidin + diet.legume$Delphinidin + diet.legume$Malvidin + diet.legume$Pelargonidin + diet.legume$Petunidin + diet.legume$Peonidin

# flavanones 
diet.legume$flavanones.legume <- diet.legume$Eriodictyol + diet.legume$Hesperetin + diet.legume$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.legume$flavan.3.ols.legume <- diet.legume$X....Catechin + diet.legume$X....Epicatechin + diet.legume$X....Epigallocatechin + diet.legume$X....Gallocatechin + diet.legume$X....Catechin.3.O.gallate + diet.legume$X....Epicatechin.3.O.gallate + diet.legume$X....Gallocatechin.3.O.gallate +  diet.legume$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.legume$flavonols.legume <- diet.legume$Quercetin + diet.legume$Kaempferol + diet.legume$Myricetin + diet.legume$Isorhamnetin

# flavones
diet.legume$flavones.legume <- diet.legume$Luteolin + diet.legume$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.legume$polymers.legume <-  diet.legume$X02.mers + diet.legume$X03.mers + diet.legume$X04.06.mers + diet.legume$X07.10.mers + diet.legume$Polymers...10.mers.

# Proanthocyanidins
diet.legume$proanthocyanidins.legume <- diet.legume$Theaflavin + diet.legume$Thearubigins

# total flavonoid
diet.legume$total.flavonoid.legume <- diet.legume$anthocyanins.legume + diet.legume$flavanones.legume + diet.legume$flavan.3.ols.legume + diet.legume$flavonols.legume + diet.legume$flavones + diet.legume$polymers.legume + diet.legume$proanthocyanidins.legume

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.legume$group <- factor(diet.legume$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.legume)[96:103]  # Selecting the variable names from column 3 to 217
description.legume <- CreateTableOne(data = diet.legume,
                                     vars = vars,  # Variable names extracted above
                                     strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.legume, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.legume

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.legume), "flavonoid/legumeresult.csv")


# ---- Rmd chunk ----
rm(diet.legume)
### Select the variables into a separate dataset for flavonoid evaluation ###
diet.nut <- data.frame(matrix(ncol = 2, nrow = nrow(demographics)))

# Select Participant ID and group (column 20) from the baseline dataset
diet.nut[,1:2] <- demographics[, c(1,27)]

# Rename the columns to match the original variable names
names(diet.nut) <- c("Participant.ID", "group")

# Define a function to sum the nut values, adjust with p.24h.foods.average, and add to diet.nut dataframe
sum.nut.values <- function(nut.rows, flavonoid.data, diet.nut) {
  # Specify the column range from column 9 to column 223
  start.col <- 9
  end.col <- 101
  
  # Extract the data from the specified column range and replace NAs with 0
  flavonoid.columns <- flavonoid.data[, start.col:end.col]
  
  # Loop through each row in diet.nut (each participant)
  for (participant in seq_len(nrow(diet.nut))) {
    
    # Initialize an empty vector to store the sum of each column for the current participant
    summed.values <- numeric(ncol(flavonoid.columns))
    
    # Loop through each specified nut row and add its values to the summed.values vector
    for (i in seq_along(nut.rows)) {
      nut.row <- nut.rows[i]
      
      # Multiply by the corresponding p.24h.foods.average value and divide by 240 for each participant
      if (nut.row == 791) {
        # (UN)Salted Nuts Intake -- Almonds (row 791)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Salted.nuts.intake[participant]) * 30 / 100 *0.2659 + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Unsalted.nuts.intake[participant]) * 30 / 100 *0.2659
      } else if (nut.row == 794) {
        # (UN)Salted Nuts Intake -- Cashew (row 794)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Salted.nuts.intake[participant]) * 30 / 100 *0.3017 + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Unsalted.nuts.intake[participant]) * 30 / 100 *0.3017
      } else if (nut.row == 796) {
        # (UN)Salted Nuts Intake -- Hazelnuts (row 796)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Salted.nuts.intake[participant]) * 30 / 100 *0.1240 + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Unsalted.nuts.intake[participant]) * 30 / 100 *0.1240
      } else if (nut.row == 811) {
        # (UN)Salted Nuts Intake --	Macadamia (row 811)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Salted.nuts.intake[participant]) * 30 / 100 *0.005 + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Unsalted.nuts.intake[participant]) * 30 / 100 *0.005
      } else if (nut.row == 805) {
        # (UN)Salted Nuts Intake --	Pecan (row 805)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Salted.nuts.intake[participant]) * 30 / 100 *0.0358 + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Unsalted.nuts.intake[participant]) * 30 / 100 *0.0358
      } else if (nut.row == 799) {
        # (UN)Salted Nuts Intake --	Pistachios (row 799)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Salted.nuts.intake[participant]) * 30 / 100 *0.0605 + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Unsalted.nuts.intake[participant]) * 30 / 100 *0.0605
      } else if (nut.row == 809) {
        # (UN)Salted Nuts Intake --	Walnuts (row 809)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Salted.nuts.intake[participant]) * 30 / 100 *0.1444 + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Unsalted.nuts.intake[participant]) * 30 / 100 *0.1444
      } else if (nut.row == 840) {
        # Seeds Intake (row 840)
        summed.values <- summed.values + (as.numeric(flavonoid.columns[nut.row, ]) * p.24h.foods.average$Seeds.intake[participant]) * 20 / 100
      }
    }
    
    # Add the resulting summed values for this participant as a new row in the diet.nut dataframe
    # Column names will correspond to the column names in flavonoid_data
    for (i in seq_along(colnames(flavonoid.columns))) {
      col.name <- colnames(flavonoid.columns)[i]
      diet.nut[participant, col.name] <- summed.values[i]
    }
  }
  
  return(diet.nut)
}

# Define the nut rows (adjust based on the actual rows you want to sum)
nut.rows <- c(791, 794, 796, 811, 805, 799, 809, 840)

# Call the function to calculate and update the diet.nut dataframe, incorporating p.24h.foods.average adjustments
diet.nut <- sum.nut.values(nut.rows, flavonoid.checkout, diet.nut)

# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(diet.nut, "flavonoid/diet.nut.csv", row.names = FALSE)


# ---- Rmd chunk ----
# anthocyanins
diet.nut$anthocyanins.nut <- diet.nut$Cyanidin + diet.nut$Delphinidin + diet.nut$Malvidin + diet.nut$Pelargonidin + diet.nut$Petunidin + diet.nut$Peonidin

# flavanones 
diet.nut$flavanones.nut <- diet.nut$Eriodictyol + diet.nut$Hesperetin + diet.nut$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
diet.nut$flavan.3.ols.nut <- diet.nut$X....Catechin + diet.nut$X....Epicatechin + diet.nut$X....Epigallocatechin + diet.nut$X....Gallocatechin + diet.nut$X....Catechin.3.O.gallate + diet.nut$X....Epicatechin.3.O.gallate + diet.nut$X....Gallocatechin.3.O.gallate +  diet.nut$X....Epigallocatechin.3.O.gallate

# flavonols # add for jaceidin, patuletin, spinacetin
diet.nut$flavonols.nut <- diet.nut$Quercetin + diet.nut$Kaempferol + diet.nut$Myricetin + diet.nut$Isorhamnetin

# flavones
diet.nut$flavones.nut <- diet.nut$Luteolin + diet.nut$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
diet.nut$polymers.nut <-  diet.nut$X02.mers + diet.nut$X03.mers + diet.nut$X04.06.mers + diet.nut$X07.10.mers + diet.nut$Polymers...10.mers.

# Proanthocyanidins
diet.nut$proanthocyanidins.nut <- diet.nut$Theaflavin + diet.nut$Thearubigins

# total flavonoid
diet.nut$total.flavonoid.nut <- diet.nut$anthocyanins.nut + diet.nut$flavanones.nut + diet.nut$flavan.3.ols.nut + diet.nut$flavonols.nut + diet.nut$flavones + diet.nut$polymers.nut + diet.nut$proanthocyanidins.nut

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
diet.nut$group <- factor(diet.nut$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(diet.nut)[96:103]  # Selecting the variable names from column 3 to 217
description.nut <- CreateTableOne(data = diet.nut,
                                  vars = vars,  # Variable names extracted above
                                  strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.nut, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.nut

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.nut), "flavonoid/nutresult.csv")

# ---- Rmd chunk ----
diet.tea <- read.csv("flavonoid/diet.tea.csv")
#diet.coffee <- read.csv("flavonoid/diet.coffee.csv")
diet.chocolate <- read.csv("flavonoid/diet.chocolate.csv")
diet.alcohol <- read.csv("flavonoid/diet.alcohol.csv")
diet.juice <- read.csv("flavonoid/diet.juice.csv")
diet.cereal <- read.csv("flavonoid/diet.cereal.csv")
#diet.vegan <- read.csv("flavonoid/diet.vegan.csv")
diet.fruit <- read.csv("flavonoid/diet.fruit.csv")
diet.vegetable <- read.csv("flavonoid/diet.vegetable.csv")
diet.legume <- read.csv("flavonoid/diet.legume.csv")
#diet.snack <- read.csv("flavonoid/diet.snack.csv")
diet.nut <- read.csv("flavonoid/diet.nut.csv")
#diet.k <- read.csv("flavonoid/necessary.csv")

# Create a new dataset by keeping the first two columns from one of the datasets
final.polyphenol <- diet.tea  # Copy the structure of dataset1

# Keep the first two columns as they are from 'dataset1'
final.polyphenol[, 1:2] <- diet.tea[, 1:2]

diet.tea[, 3:95] <- replace(diet.tea[, 3:95], is.na(diet.tea[, 3:95]), 0)
diet.chocolate[, 3:95] <- replace(diet.chocolate[, 3:95], is.na(diet.chocolate[, 3:95]), 0)
diet.alcohol[, 3:95] <- replace(diet.alcohol[, 3:95], is.na(diet.alcohol[, 3:95]), 0)
diet.juice[, 3:95] <- replace(diet.juice[, 3:95], is.na(diet.juice[, 3:95]), 0)
diet.fruit[, 3:95] <- replace(diet.fruit[, 3:95], is.na(diet.fruit[, 3:95]), 0)
diet.vegetable[, 3:95] <- replace(diet.vegetable[, 3:95], is.na(diet.vegetable[, 3:95]), 0)
#diet.snack[, 3:95] <- replace(diet.snack[, 3:95], is.na(diet.snack[, 3:95]), 0)
#diet.coffee[, 3:95] <- replace(diet.coffee[, 3:95], is.na(diet.coffee[, 3:95]), 0)
diet.cereal[, 3:95] <- replace(diet.cereal[, 3:95], is.na(diet.cereal[, 3:95]), 0)
#diet.vegan[, 3:95] <- replace(diet.vegan[, 3:95], is.na(diet.vegan[, 3:95]), 0)
diet.legume[, 3:95] <- replace(diet.legume[, 3:95], is.na(diet.legume[, 3:95]), 0)
diet.nut[, 3:95] <- replace(diet.nut[, 3:95], is.na(diet.nut[, 3:95]), 0)

# Add the values of corresponding cells in columns 3 to 217 of both datasets
final.polyphenol[, 3:95] <- diet.tea[, 3:95] + 
  diet.chocolate[, 3:95] + 
  diet.alcohol[, 3:95] + 
  diet.juice[, 3:95] + 
  diet.fruit[, 3:95] + 
  diet.vegetable[, 3:95] #+ 
#diet.snack[, 3:95] +
#diet.coffee[, 3:95] + 
diet.cereal[, 3:95] + 
  #diet.vegan[, 3:95] + 
  diet.legume[, 3:95] +
  diet.nut[, 3:95]

#final.polyphenol[, 3:95] <- diet.k[, 3:95]

# The result is stored in 'summed_dataset', with columns 3 to 217 containing the summed values
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(final.polyphenol, "flavonoid/final.polyphenol.csv", row.names = FALSE)

# ---- Rmd chunk ----
# anthocyanins
final.polyphenol$anthocyanins1 <- final.polyphenol$Cyanidin + final.polyphenol$Delphinidin + final.polyphenol$Malvidin + final.polyphenol$Pelargonidin + final.polyphenol$Petunidin + final.polyphenol$Peonidin

# flavanones 
final.polyphenol$flavanones1 <- final.polyphenol$Eriodictyol + final.polyphenol$Hesperetin + final.polyphenol$Naringenin

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
final.polyphenol$flavan.3.ols1 <- final.polyphenol$X....Catechin + final.polyphenol$X....Epicatechin 
#+ final.polyphenol$X....Epigallocatechin + final.polyphenol$X....Gallocatechin + final.polyphenol$X....Catechin.3.O.gallate + final.polyphenol$X....Epicatechin.3.O.gallate + final.polyphenol$X....Gallocatechin.3.O.gallate +  final.polyphenol$X....Epigallocatechin.3.O.gallate


# flavonols # add for jaceidin, patuletin, spinacetin
final.polyphenol$flavonols1 <- final.polyphenol$Quercetin + final.polyphenol$Kaempferol + final.polyphenol$Myricetin + final.polyphenol$Isorhamnetin

# flavones
final.polyphenol$flavones1 <- final.polyphenol$Luteolin + final.polyphenol$Apigenin

# polymers (including proanthocyanidins [excluding monomers], theaflavins, and thearubigins), proanthocyanidins (dimers, trimers, 4–6mers, 7–10mers, polymers, and monomers). 
final.polyphenol$polymers1 <-  final.polyphenol$X02.mers + final.polyphenol$X03.mers + final.polyphenol$X04.06.mers + final.polyphenol$X07.10.mers + final.polyphenol$Polymers...10.mers.

# Proanthocyanidins
final.polyphenol$proanthocyanidins1 <- final.polyphenol$Theaflavin + final.polyphenol$Thearubigins

# total flavonoid
final.polyphenol$total.flavonoid1 <- final.polyphenol$anthocyanins1 + final.polyphenol$flavanones1 + final.polyphenol$flavan.3.ols1 + final.polyphenol$flavonols1 + final.polyphenol$flavones1 + final.polyphenol$polymers1 + final.polyphenol$proanthocyanidins1

# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
final.polyphenol$group <- factor(final.polyphenol$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the tea variables (columns 3 to 217) stratified by the 'group'
vars <- names(final.polyphenol)[96:103]  # Selecting the variable names from column 3 to 217
description.polyphenol <- CreateTableOne(data = final.polyphenol,
                                         vars = vars,  # Variable names extracted above
                                         strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.polyphenol, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.polyphenol

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.polyphenol), "flavonoid/polyphenolresult.csv")

# ---- Rmd chunk ----
### Select the variables into a separate dataset for flavonoid evaluation ###
selected.polyphenol <- data.frame(matrix(ncol = 2, nrow = nrow(diet.tea)))

# Select Participant ID and group (column 20) from the baseline dataset
selected.polyphenol[,1:2] <- diet.tea[, c(1,2)]

# Rename the columns to match the original variable names
names(selected.polyphenol) <- c("Participant.ID", "group")


# anthocyanins
selected.polyphenol$anthocyanins <- final.polyphenol$anthocyanins1

# flavanones 
selected.polyphenol$flavanones <- final.polyphenol$flavanones1

# flavonols # add for jaceidin, patuletin, spinacetin
selected.polyphenol$flavonols <- final.polyphenol$flavonols1

# flavones
selected.polyphenol$flavones <- final.polyphenol$flavones1

# flavan-3-ols # add for Epigallocatechin Gallocatechin Catechin.3.O.gallate Epicatechin.3.O.gallate Epigallocatechin.3.O.gallate
selected.polyphenol$flavan.3.ols <- final.polyphenol$flavan.3.ols1 + final.polyphenol$polymers1 + final.polyphenol$proanthocyanidins1


# total flavonoid
selected.polyphenol$total.flavonoid <- selected.polyphenol$anthocyanins + selected.polyphenol$flavanones + selected.polyphenol$flavan.3.ols + selected.polyphenol$flavonols + selected.polyphenol$flavones


# The result is stored in 'summed_dataset', with columns 3 to 217 containing the summed values
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(selected.polyphenol, "flavonoid/selected.polyphenol.csv", row.names = FALSE)

# ---- Rmd chunk ----
selected.polyphenol <- read.csv("flavonoid/selected.polyphenol.csv")
#selected.polyphenol <- selected.polyphenol %>% filter(Participant.ID %in% polyphenolwithst$Participant.ID)
# Set the 'group' column as a factor with specified levels: 'psoriasis' and 'non-psoriasis'
selected.polyphenol$group <- factor(selected.polyphenol$group, levels = c("psoriasis", "non-psoriasis"))

# Create a descriptive table for the fruit variables (columns 3 to 217) stratified by the 'group'
vars <- names(selected.polyphenol)[3:8]  # Selecting the variable names from column 3 to 217
description.polyphenol <- CreateTableOne(data = selected.polyphenol,
                                         vars = vars,  # Variable names extracted above
                                         strata = "group")  # Stratified by 'group' column

# Print the description table without quotes, without spaces, and without printing to the console
print(description.polyphenol, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
description.polyphenol

# Save the descriptive table results to a CSV file
# Publication note: prefer readr::write_csv for clean reproducible CSV output
write.csv(print(description.polyphenol), "flavonoid/selectedpolyphenolresult.csv")

# ---- Rmd chunk ----
# Step 1: Specify flavonoid subclass column names
# Extract column names from any one dataset (e.g., diet.nut)
flavonoid_columns <- names(diet.nut)[96:103]

# Step 2: Reconstruct the source dataset list
# Keep only Participant.ID and flavonoid variables
source_list <- list(
  nut       = diet.nut[, c("Participant.ID", flavonoid_columns)],
  tea       = diet.tea[, c("Participant.ID", flavonoid_columns)],
  chocolate = diet.chocolate[, c("Participant.ID", flavonoid_columns)],
  alcohol   = diet.alcohol[, c("Participant.ID", flavonoid_columns)],
  juice     = diet.juice[, c("Participant.ID", flavonoid_columns)],
  cereal    = diet.cereal[, c("Participant.ID", flavonoid_columns)],
  fruit     = diet.fruit[, c("Participant.ID", flavonoid_columns)],
  vegetable = diet.vegetable[, c("Participant.ID", flavonoid_columns)],
  legume    = diet.legume[, c("Participant.ID", flavonoid_columns)]
)

# Step 3: Extract quartile information
quartile_info <- final.polyphenol[, c("Participant.ID", "quartile_total.flavonoid")]

# Step 4: Initialise an empty dataset
final_long_data <- data.frame()

# Step 5: Loop through each food source
# Merge quartile information and reshape to long format
for (source_name in names(source_list)) {
  
  df <- source_list[[source_name]]
  
  # Merge quartile information
  df <- merge(df, quartile_info, by = "Participant.ID")
  
  # Add food source label
  df$Source <- source_name
  
  # Convert to long format
  temp_long <- tidyr::pivot_longer(
    df,
    cols = all_of(flavonoid_columns),
    names_to = "SubType",
    values_to = "Amount"
  )
  
  # Append to the final dataset
  final_long_data <- rbind(final_long_data, temp_long)
}

# Calculate the difference between Q4 and Q3
q3_q4_diff <- final_long_data %>%
  group_by(Source, SubType, quartile_total.flavonoid) %>%
  summarise(mean = mean(Amount, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = quartile_total.flavonoid,
    values_from = mean,
    names_prefix = "Q"
  ) %>%
  mutate(diff = Q4 - Q3) %>%
  arrange(desc(diff))

# Export results
# Publication note: use readr::write_csv for cleaner and reproducible output
write.csv(q3_q4_diff, "q3_q4_diff.csv")