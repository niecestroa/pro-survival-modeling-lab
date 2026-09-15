###############################################################################
# Title: GBCS Survival Analysis + ASA RandomForestSRC Compatibility (Part II)
# Author:       Aaron Niecestro
# Created on:   September 15, 2026
# Last Edit:    September 15, 2026
#
# Description:
# This script prepares the GBCS breast‑cancer dataset for use with the ASA
# Traveling Course: Tree-Based Machine Learning Methods (Part II: Inference
# and Prediction). The dataset is cleaned, recoded, and transformed to ensure
# compatibility with randomForestSRC functions used throughout the ASA
# workshop, including:
#
#   • OOB classification inference (predicted, predicted.oob, class.oob)
#   • OOB survival inference (survival, survival.oob, chf, chf.oob)
#   • Prediction on new data (predict.rfsrc)
#   • Train/test splits with factor-level mismatches
#   • Restore mode (perf.type = "brier", forest.wt, custom estimators)
#   • Partial dependence plots (plot.variable, partial)
#
# The GBCS dataset is modified to include:
#   • Surv(time, status) columns for survival forests
#   • y_class for classification forests
#   • y_reg for regression forests
#   • Removal of unsupported columns (IDs, dates, raw encodings)
#   • Proper factor encoding for categorical variables
#   • No missing values (required by randomForestSRC)
#
# Citation:
# Portions of this workflow reference the ASA Traveling Course:
# “Tree-Based Machine Learning Methods,” American Statistical Association.
# Code examples adapted from the ASA student R-code companion.
###############################################################################

###############################################################################
# Prepare GBCS dataset for ASA Traveling Course Part II (Inference & Prediction)
###############################################################################

library(tidyverse)
library(survival)
library(randomForestSRC)

# Load original GBCS dataset
gbcs <- read_csv("~/BIST0627 Applied Survival Data Analysis/BIST0627 Project/project_data_files/gbcs.csv")

# --- STEP 1: Clean and recode for ASA compatibility --------------------------

gbcs_asa <- gbcs %>%
  mutate(
    # Survival forest requirements
    time   = as.numeric(survtime),              # rename survtime → time
    status = case_when(
      censdead == 1 ~ 1,                        # event occurred
      censdead == 0 ~ 0,                        # censored
      TRUE ~ 0
    ),
    
    # Classification target (Glioma-style)
    y_class = factor(if_else(censdead == 1, "Dead", "Alive")),
    
    # Regression target (SalePrice-style)
    y_reg = as.numeric(survtime),
    
    # Numeric predictors
    age        = as.numeric(age),
    size       = as.numeric(size),
    nodes      = as.numeric(nodes),
    prog_recp  = as.numeric(prog_recp),
    estrg_recp = as.numeric(estrg_recp),
    rectime    = as.numeric(rectime),
    
    # Factor predictors
    menopause_f = factor(if_else(menopause == "1", "Yes", "No")),
    hormone_f   = factor(if_else(hormone == "1", "Yes", "No")),
    grade_f     = factor(grade),
    censrec_f   = factor(if_else(censrec == "1", "Recurrence", "Censored"))
  ) %>%
  # Remove columns ASA forests cannot handle
  select(
    -id, -diagdateb, -recdate, -deathdate,
    -menopause, -hormone, -grade, -censrec
  ) %>%
  # Remove missing values (required by randomForestSRC)
  drop_na()

# --- STEP 2: Verify ASA compatibility ----------------------------------------

str(gbcs_asa)       # must show only numeric + factor columns
anyNA(gbcs_asa)     # must be FALSE

# --- STEP 3: Test ASA Part II functionality ----------------------------------

# 1. Classification forest (Glioma example)
fit_class <- rfsrc(y_class ~ ., data = gbcs_asa)
print(fit_class)

# 2. Survival forest (PBC example)
fit_surv <- rfsrc(Surv(time, status) ~ ., data = gbcs_asa)
print(fit_surv)

# 3. Regression forest (SalePrice example)
fit_reg <- rfsrc(y_reg ~ ., data = gbcs_asa)
print(fit_reg)

# 4. OOB quantities (predicted, predicted.oob, survival.oob, chf.oob)
head(fit_surv$survival.oob)
head(fit_class$predicted.oob)

# 5. Prediction on new data (veteran example)
test_idx <- sample(1:nrow(gbcs_asa), 10)
pred_surv <- predict(fit_surv, gbcs_asa[test_idx, ])
pred_class <- predict(fit_class, gbcs_asa[test_idx, ])
pred_reg <- predict(fit_reg, gbcs_asa[test_idx, ])

# 6. Restore mode (perf.type = "brier")
restore_surv <- predict(fit_surv, perf.type = "brier")
restore_class <- predict(fit_class, perf.type = "brier")

# 7. Partial plots (peakVO2 example)
plot.variable(fit_surv, surv.type = "surv", xvar.names = "age",
              time = 365, partial = TRUE)
