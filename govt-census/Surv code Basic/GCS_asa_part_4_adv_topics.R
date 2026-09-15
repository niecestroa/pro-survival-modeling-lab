###############################################################################
# Author:       Aaron Niecestro
# DATE CREATED: February 15, 2026
# DATE EDITED:  February 15, 2026
#
# FILE DESCRIPTION:
#
# This script adapts the ASA Traveling Course:
# "Tree-Based Machine Learning Methods — Part IV: Advanced Topics"
# for use with the German PBC dataset (PBC276). The original workshop examples
# include imbalanced classification, missing-data imputation, Super Greedy Trees
# (SGT), and Random Hazard Forests (RHF). This version rewrites all sections
# that can logically operate on survival data (PBC276) and fully comments the
# remaining examples that require incompatible data structures.
#
# Adapted sections for PBC276:
#   • Survival imputation (supervised + unsupervised)
#   • Test-time imputation using impute.learn()
#   • Out-of-distribution (OOD) scoring
#
# Sections kept as ASA originals (with full comments):
#   • Imbalanced classification (PBC276 is not a classification dataset)
#   • Super Greedy Trees (SGT) — requires regression or high-dimensional data
#   • Random Hazard Forests (RHF) — requires counting-process longitudinal data
#
# Origin: ASA Traveling Course — Part IV: Advanced Topics
# Adaptation: Aaron Niecestro (German PBC dataset)
###############################################################################

#===============================================================================
# Load required packages
#===============================================================================

library(randomForestSRC)     # RSF, impute(), OOD scoring
library(survival)            # Surv() object
library(tidyverse)           # data manipulation
library(janitor)             # clean_names()
library(varPro)              # optional: subsample(), VarPro tools

#===============================================================================
# Load and prepare the German PBC dataset (PBC276)
#===============================================================================

PBC276 <- read_csv("~/UTH Survival Analysis/PH1831_Project/PH1831 Project Data/PBC276.csv") %>%
  clean_names() %>%                     # convert column names to snake_case
  mutate(
    stage = factor(stage),              # categorical variables → factors
    drug  = factor(drug),
    sex   = factor(sex)
  )

glimpse(PBC276)                         # inspect dataset structure

###############################################################################
# SECTION 1 — Imbalanced Classification (NOT APPLICABLE TO PBC276)
###############################################################################
# The ASA example uses glioma classification. PBC276 is a survival dataset
# (time-to-event), not a classification dataset. Therefore, this section is
# preserved as comments only.

# Example (kept for reference):
# library(varPro)
# data(glioma, package = "varPro")
# table(glioma$y)
# ... classification examples omitted ...

###############################################################################
# SECTION 2 — Missing-data imputation (ADAPTED FOR PBC276)
###############################################################################

# Supervised on-the-fly imputation for survival data
pbc_sup_imp <- impute(
  Surv(futime, status) ~ .,             # supervised survival formula
  data = PBC276                         # German PBC dataset
)

# Unsupervised on-the-fly imputation
pbc_unsup_imp <- impute(
  data = PBC276                         # no formula → unsupervised imputation
)

# Print imputed dataset summaries
summary(pbc_sup_imp)
summary(pbc_unsup_imp)

###############################################################################
# SECTION 3 — missForest and mForest imputation (ADAPTED FOR PBC276)
###############################################################################

# missForest-style imputation (one variable at a time)
pbc_missforest <- impute(
  data = PBC276,
  mf.q = 1                              # mf.q = 1 → missForest mode
)

# mForest-style grouped multivariate imputation
pbc_mforest <- impute(
  data = PBC276,
  mf.q = 0.5                            # mf.q between 0 and 1 → grouped regressions
)

# Print summaries
summary(pbc_missforest)
summary(pbc_mforest)

###############################################################################
# SECTION 4 — Test-time imputation using impute.learn() (ADAPTED FOR PBC276)
###############################################################################

# Create a train/test split
set.seed(1831)
id <- sample(seq_len(nrow(PBC276)), 100)
train <- PBC276[id, ]
test  <- PBC276[-id, ]

# Learn imputation system from training data
fit_imp <- impute.learn(
  data = train,
  mf.q = 1,                             # missForest-style
  max.iter = 5,                         # number of imputation iterations
  full.sweep.options = list(
    ntree = 25,                         # number of trees per sweep
    nsplit = 5                          # number of split candidates
  ),
  target.mode = "all"                   # impute all variables
)

# Apply learned imputation to test data
test_imp <- predict(
  fit_imp,
  newdata = test,
  max.predict.iter = 2                  # number of prediction iterations
)

# Print imputed test data
summary(test_imp)

###############################################################################
# SECTION 5 — Out-of-distribution (OOD) scoring (ADAPTED FOR PBC276)
###############################################################################

# Train supervised imputation model with OOD scoring enabled
fit_ood <- impute.learn(
  data = train,
  mf.q = 1,
  supervised.formula = Surv(futime, status) ~ .,   # supervised survival model
  supervised.args = list(
    ntree = 50,
    nsplit = 5
  ),
  full.sweep.options = list(
    ntree = 25,
    nsplit = 5
  ),
  save.ood = TRUE                                  # enable OOD scoring
)

# Score and impute test cases
ood_scores <- impute.ood(fit_ood, test)

# Print OOD scores and percentiles
head(ood_scores$score)
head(ood_scores$score.percentile)

###############################################################################
# SECTION 6 — Super Greedy Trees (SGT) (NOT APPLICABLE TO PBC276)
###############################################################################
# SGT requires regression or high-dimensional continuous data.
# PBC276 is survival data → cannot be used directly.
# Section preserved as comments only.

# Example (kept for reference):
# library(randomForestSGT)
# n <- 2500
# p <- 50
# noise <- matrix(runif(n * p), ncol = p)
# dta <- data.frame(mlbench:::mlbench.friedman1(n, sd = 0), noise = noise)
# filter <- tune.hcut(y ~ ., data = dta, hcut = 3)
# o.sgt <- rfsgt(y ~ ., data = dta, filter = filter)

###############################################################################
# SECTION 7 — Random Hazard Forests (RHF) (NOT APPLICABLE TO PBC276)
###############################################################################
# RHF requires longitudinal counting-process data:
#   id, start, stop, event
# PBC276 is one-row-per-subject baseline survival → cannot be converted.
# Section preserved as comments only.

# Example (kept for reference):
# library(randomForestRHF)
# d <- convert.counting(Surv(ttodead, died) ~ ., data = peakVO2)
# f <- "Surv(id, start, stop, event) ~ ."
# o <- rhf(f, data = d)

###############################################################################
# End of ASA Part IV — PBC276 Adaptation
###############################################################################
