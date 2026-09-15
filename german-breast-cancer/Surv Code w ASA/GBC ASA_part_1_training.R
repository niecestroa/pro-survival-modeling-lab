###############################################################################
# Title: GBCS Survival Analysis + ASA RandomForestSRC Compatibility
# Author:       Aaron Niecestro
# Date Created: September 15, 2026
# Last Edit:    September 15, 2026
#
# Description:
# This script prepares the GBCS breast‑cancer dataset for use with both
# traditional survival modeling (Cox, Weibull, diagnostics) and the ASA
# Traveling Course tree‑based machine‑learning methods (randomForestSRC).
# The dataset is cleaned, recoded, and transformed to ensure compatibility
# with rfsrc() for regression, classification, and survival forests.
#
# Citation:
# Portions of this workflow reference the ASA Traveling Course:
# “Tree‑Based Machine Learning Methods,” American Statistical Association.
# Code examples adapted from the ASA student R‑code companion.

# Part I: Training
###############################################################################

###############################################################################
# GBCS Dataset Preparation for ASA Traveling Course (randomForestSRC)
###############################################################################

library(tidyverse)
library(survival)
library(randomForestSRC)

# Load your original cleaned dataset (from your modular pipeline)
gbcs <- read_csv("~/BIST0627 Applied Survival Data Analysis/BIST0627 Project/project_data_files/gbcs.csv")

# --- STEP 1: Clean and recode for ASA compatibility --------------------------

gbcs_asa <- gbcs %>%
  mutate(
    # Convert survival time to ASA-required name
    time = as.numeric(survtime),        # ASA forests require "time" for Surv(time, status)
    
    # Convert death indicator to ASA-required name
    status = case_when(
      censdead == 1 ~ 1,                # event occurred
      censdead == 0 ~ 0,                # censored
      TRUE ~ 0                          # fallback (should not occur)
    ),
    
    # Ensure all numeric predictors are numeric
    age = as.numeric(age),
    size = as.numeric(size),
    nodes = as.numeric(nodes),
    prog_recp = as.numeric(prog_recp),
    estrg_recp = as.numeric(estrg_recp),
    rectime = as.numeric(rectime),
    
    # Convert categorical variables to factors
    menopause_f = factor(if_else(menopause == "1", "Yes", "No")),
    hormone_f   = factor(if_else(hormone == "1", "Yes", "No")),
    grade_f     = factor(grade),
    censrec_f   = factor(if_else(censrec == "1", "Recurrence", "Censored")),
    
    # Create ASA-friendly classification target
    y_class = factor(if_else(censdead == 1, "Dead", "Alive")),
    
    # Create ASA-friendly regression target
    y_reg = as.numeric(survtime)
  ) %>%
  # Remove columns ASA forests cannot handle
  select(
    -id, -diagdateb, -recdate, -deathdate,
    -menopause, -hormone, -grade, -censrec
  ) %>%
  # Remove rows with missing values (ASA requirement)
  drop_na()

# --- STEP 2: Verify ASA compatibility ----------------------------------------

str(gbcs_asa)       # Should show only numeric + factor columns
summary(gbcs_asa)   # No missing values allowed
anyNA(gbcs_asa)     # Must be FALSE

# --- STEP 3: Test ASA forest calls -------------------------------------------

# Survival forest (ASA Slide 6 style)
fit_surv <- rfsrc(Surv(time, status) ~ ., data = gbcs_asa)
print(fit_surv)

# Regression forest (SalePrice example → use y_reg)
fit_reg <- rfsrc(y_reg ~ ., data = gbcs_asa)
print(fit_reg)

# Classification forest (Glioma example → use y_class)
fit_class <- rfsrc(y_class ~ ., data = gbcs_asa)
print(fit_class)

# Integrated ASA analysis
run.rfsrc(Surv(time, status) ~ ., data = gbcs_asa)

