###############################################################################
# Title: GBCS Survival Analysis + ASA RandomForestSRC Compatibility (Part III)
# Author: Aaron Niecestro
# Created on: January 8, 2024
# Last Edited on: February 15, 2026
#
# Description:
# This script prepares the GBCS breast‑cancer dataset for use with the ASA
# Traveling Course: Tree-Based Machine Learning Methods (Part III: Variable
# Selection). The dataset is cleaned, recoded, and transformed to ensure
# compatibility with all variable‑selection methods demonstrated in the ASA
# workshop, including:
#
#   • Permutation VIMP (Breiman–Cutler)
#   • Anti-VIMP and random-VIMP
#   • Block-size VIMP
#   • Joint VIMP for variable pairs
#   • Subsampling inference for VIMP confidence intervals
#   • Minimal depth variable selection
#   • Guided trees using split weights
#   • VarPro (supervised variable priority)
#   • cv.varpro (cross‑validated variable priority)
#   • ivarpro (individual variable priority)
#   • shap.ivarpro (SHAP-style variable priority)
#   • uvarpro (unsupervised variable priority)
#   • High-dimensional VarPro examples
#
# The GBCS dataset is modified to include:
#   • Surv(time, status) columns for survival forests
#   • y_class for classification forests
#   • y_reg for regression forests
#   • Removal of unsupported columns (IDs, dates, raw encodings)
#   • Proper factor encoding for categorical variables
#   • No missing values (required by randomForestSRC and VarPro)
#
# Citation:
# Portions of this workflow reference the ASA Traveling Course:
# “Tree-Based Machine Learning Methods,” American Statistical Association.
# Code examples adapted from the ASA student R-code companion.
###############################################################################

###############################################################################
# Prepare GBCS dataset for ASA Traveling Course Part III (Variable Selection)
###############################################################################

library(tidyverse)
library(survival)
library(randomForestSRC)
library(varPro)

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
    
    # Classification target (Species / Glioma style)
    y_class = factor(if_else(censdead == 1, "Dead", "Alive")),
    
    # Regression target (mpg / SalePrice style)
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
  # Remove missing values (required by randomForestSRC + VarPro)
  drop_na()

# --- STEP 2: Verify ASA compatibility ----------------------------------------

str(gbcs_asa)       # must show only numeric + factor columns
anyNA(gbcs_asa)     # must be FALSE

# --- STEP 3: Test ASA Part III functionality ---------------------------------

# 1. Permutation VIMP (Slide 7)
vimp_test <- rfsrc(y_reg ~ ., data = gbcs_asa, importance = "permute")$importance

# 2. Joint VIMP (Slide 7)
obj <- rfsrc(y_reg ~ ., data = gbcs_asa)
joint_vimp <- vimp(obj, c("age", "size"), joint = TRUE, importance = "permute")

# 3. Subsampling inference (Slide 11)
o_surv <- rfsrc(Surv(time, status) ~ ., gbcs_asa, importance = "permute")
oo <- subsample(o_surv)

# 4. Minimal depth (Slide 18)
md <- max.subtree(o_surv)$order[, 1]

# 5. Guided trees using split weights (Slide 19)
xvar.used <- predict(o_surv, var.used = "all.trees")$var.used
os <- rfsrc(Surv(time, status) ~ ., gbcs_asa, xvar.wt = xvar.used)
md_guided <- max.subtree(os)$order[, 1]

# 6. VarPro (Slide 26)
vp <- varpro(Surv(time, status) ~ ., gbcs_asa)
vp_cv <- cv.varpro(Surv(time, status) ~ ., gbcs_asa)

# 7. ivarpro + shap.ivarpro (Slide 34)
ivp <- ivarpro(vp)
ivp_shap <- shap.ivarpro(ivp)

# 8. uvarpro (Slide 38)
uvp <- uvarpro(gbcs_asa)

# All ASA Part III examples now run successfully on GBCS.
