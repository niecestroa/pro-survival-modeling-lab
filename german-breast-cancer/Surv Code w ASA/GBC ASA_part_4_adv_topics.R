###############################################################################
# Title: GBCS Survival Analysis + ASA RandomForestSRC Compatibility (Part IV)
# Author:       Aaron Niecestro
# Created on:   September 15, 2024
# Last Edit:    September 15, 2026
#
# Description:
# This script prepares the GBCS breast‑cancer dataset for use with the ASA
# Traveling Course: Tree-Based Machine Learning Methods (Part IV: Advanced
# Topics). The dataset is cleaned, recoded, and transformed to ensure
# compatibility with all advanced methods demonstrated in the ASA workshop,
# including:
#
#   • Imbalanced classification (RFQ, BRF, G‑mean VIMP)
#   • Subsampling confidence intervals for VIMP
#   • Missing‑data imputation (supervised, unsupervised)
#   • missForest and mForest grouped regressions
#   • impute.learn + predict for test‑time imputation
#   • Out‑of‑distribution (OOD) scoring
#   • Super Greedy Trees (SGT) and hcut tuning
#   • SGT model explainers (beta coefficients, partial effects)
#   • Random Hazard Forests (RHF)
#   • Counting‑process survival format conversion
#   • Time‑dependent AUC curves
#   • Smoothed hazard curves
#   • Time‑localized RHF variable importance
#
# The GBCS dataset is modified to include:
#   • Surv(time, status) columns for survival forests and RHF
#   • y_class for classification forests (including imbalanced examples)
#   • y_reg for regression forests
#   • Removal of unsupported columns (IDs, dates, raw encodings)
#   • Proper factor encoding for categorical variables
#   • No missing values (required by randomForestSRC, SGT, RHF, VarPro)
#
# Citation:
# Portions of this workflow reference the ASA Traveling Course:
# “Tree-Based Machine Learning Methods,” American Statistical Association.
# Code examples adapted from the ASA student R-code companion.
###############################################################################

###############################################################################
# Prepare GBCS dataset for ASA Traveling Course Part IV (Advanced Topics)
###############################################################################

library(tidyverse)
library(survival)
library(randomForestSRC)
library(varPro)
library(randomForestSGT)
library(randomForestRHF)

# Load original GBCS dataset
gbcs <- read_csv("~/BIST0627 Applied Survival Data Analysis/BIST0627 Project/project_data_files/gbcs.csv")

# --- STEP 1: Clean and recode for ASA compatibility --------------------------

gbcs_asa <- gbcs %>%
  mutate(
    # Survival forest + RHF requirements
    time   = as.numeric(survtime),              # rename survtime → time
    status = case_when(
      censdead == 1 ~ 1,                        # event occurred
      censdead == 0 ~ 0,                        # censored
      TRUE ~ 0
    ),
    
    # Classification target (imbalanced examples)
    y_class = factor(if_else(censdead == 1, "Dead", "Alive")),
    
    # Regression target (SGT, RFQ, missForest)
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
  # Remove missing values (required by randomForestSRC, SGT, RHF)
  drop_na()

# --- STEP 2: Verify ASA compatibility ----------------------------------------

str(gbcs_asa)       # must show only numeric + factor columns
anyNA(gbcs_asa)     # must be FALSE

# --- STEP 3: Test ASA Part IV functionality ----------------------------------

# 1. Imbalanced classification (RFQ, BRF)
o_rf <- imbalanced(y_class ~ ., data = gbcs_asa)
o_brf <- imbalanced(y_class ~ ., data = gbcs_asa, method = "brf")

# 2. G‑mean VIMP + subsampling
o_vimp <- imbalanced(y_class ~ ., data = gbcs_asa,
                     importance = "permute", block.size = 20)
oo_vimp <- subsample(o_vimp)

# 3. Missing‑data imputation (supervised + unsupervised)
imp_unsup <- impute(data = gbcs_asa)
imp_sup <- impute(Surv(time, status) ~ ., data = gbcs_asa)

# 4. missForest / mForest
imp_mf1 <- impute(data = gbcs_asa, mf.q = 1)
imp_mf40 <- impute(data = gbcs_asa, mf.q = 40)

# 5. impute.learn + predict
fit_imp <- impute.learn(
  data = gbcs_asa,
  mf.q = 1,
  max.iter = 5,
  full.sweep.options = list(ntree = 25, nsplit = 5),
  target.mode = "all"
)
pred_imp <- predict(fit_imp, gbcs_asa[1:20, ])

# 6. OOD scoring
fit_sup <- impute.learn(
  data = gbcs_asa,
  mf.q = 1,
  supervised.formula = y_reg ~ .,
  supervised.args = list(ntree = 50, nsplit = 5),
  full.sweep.options = list(ntree = 25, nsplit = 5),
  save.ood = TRUE
)
ood_scores <- impute.ood(fit_sup, gbcs_asa[1:20, ])

# 7. Super Greedy Trees (SGT)
sgt_filter <- tune.hcut(y_reg ~ ., data = gbcs_asa, hcut = 3)
sgt_fit <- rfsgt(y_reg ~ ., data = gbcs_asa, filter = sgt_filter)

# 8. SGT explainers
sgt_beta <- get.beta(sgt_fit, bag = "oob")

# 9. Random Hazard Forests (RHF)
# Convert to counting-process format
gbcs_count <- convert.counting(Surv(time, status) ~ ., data = gbcs_asa)

# Fit RHF
f_rhf <- "Surv(id, start, stop, event) ~ ."
rhf_fit <- rhf(f_rhf, data = gbcs_count)

# Time-dependent AUC
auc_chf <- auct.rhf(rhf_fit)
auc_haz <- auct.rhf(rhf_fit, marker = "haz")

# Smoothed hazards
shaz <- smoothed.hazard(rhf_fit)

# Time-localized RHF importance
imp_rhf <- importance.rhf(rhf_fit)

# All ASA Part IV examples now run successfully on GBCS.
