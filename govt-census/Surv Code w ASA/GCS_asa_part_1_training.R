#======================================================================
# Author:       Aaron Niecestro
# Created:      September 14, 2026
# Last Edit:    September 15, 2026
#
# Project: US Government Survival Analysis – PBC276
# 
# FILE DESCRIPTION
#
# This script is adapted from the **ASA Traveling Course: Tree-Based Machine
# Learning Methods — Part I: Training**. The original workshop code demonstrates
# how to train Random Forest models for regression, classification, and survival
# analysis using the randomForestSRC ecosystem.
#
# This version rewrites the survival-analysis portion of the ASA Part I training
# code so that it works with the **PBC276 dataset** used in your clinical
# survival modeling project. The Mayo Clinic PBC dataset from the workshop is
# replaced with your PBC276 dataset, and all syntax is updated accordingly.
#
# The workflow includes:
#
#   1. Data Import & Cleaning
#      - Loads PBC276 from local storage
#      - Cleans column names and converts categorical variables to factors
#
#   2. Survival Object Construction
#      - Builds a right-censored Surv(futime, status) object
#      - Uses status = 0/1 (censored/event)
#
#   3. Random Forest Survival Model Training
#      - Fits a baseline survival forest using logrank splitting
#      - Trains alternative forests using logrankscore and bs.gradient rules
#      - Prints model summaries including OOB prediction error
#
#   4. Prediction
#      - Generates survival predictions for the first 10 patients
#      - Extracts predicted survival curves
#
#   5. Integrated Analysis
#      - Runs run.rfsrc() for a full automated workflow
#
# This script provides a complete machine-learning survival pipeline parallel to
# your Cox proportional hazards workflow, enabling comparison of parametric vs.
# nonparametric survival modeling approaches.

###############################################################################

###############################################################################
# Survival example: PBC276 (your dataset)
###############################################################################

library(randomForestSRC)
library(survival)
library(tidyverse)
library(janitor)

#----------------------------------------------------------------------
# Load and prepare your dataset
#----------------------------------------------------------------------

PBC276 <- read_csv("~/UTH Survival Analysis/PH1831_Project/PH1831 Project Data/PBC276.csv") %>%
  clean_names() %>%
  mutate(
    stage = factor(stage),
    drug  = factor(drug),
    sex   = factor(sex)
  )

# Inspect structure
glimpse(PBC276)

#----------------------------------------------------------------------
# Create right‑censored survival object
#----------------------------------------------------------------------

# Your dataset uses:
#   futime = follow‑up time
#   status = 0/1 event indicator (1 = death)

# No competing risks here, so Surv() is straightforward:
surv_obj <- Surv(PBC276$futime, PBC276$status)

#----------------------------------------------------------------------
# Fit survival random forest
#----------------------------------------------------------------------

o <- rfsrc(
  Surv(futime, status) ~ .,     # use all predictors
  data = PBC276,
  ntree = 1000,                 # more stable survival forest
  splitrule = "logrank"         # default survival split rule
)

print(o)

#----------------------------------------------------------------------
# Alternative split rules (as in workshop slides)
#----------------------------------------------------------------------

# Log-rank score split rule
o_lrscore <- rfsrc(
  Surv(futime, status) ~ ., 
  data = PBC276,
  splitrule = "logrankscore"
)

print(o_lrscore)

# Gradient-based split rule
o_bsgrad <- rfsrc(
  Surv(futime, status) ~ ., 
  data = PBC276,
  splitrule = "bs.gradient"
)

print(o_bsgrad)

#----------------------------------------------------------------------
# Prediction example (like Slide 20)
#----------------------------------------------------------------------

# Predict survival for first 10 patients
pred10 <- predict(o, newdata = PBC276[1:10, ])

# Extract predicted survival curves
head(pred10$survival)

#----------------------------------------------------------------------
# Integrated analysis (Slide 47 equivalent)
#----------------------------------------------------------------------

library(randomForestSRC.run)

run.rfsrc(
  Surv(futime, status) ~ ., 
  data = PBC276,
  ntree = 1000
)
