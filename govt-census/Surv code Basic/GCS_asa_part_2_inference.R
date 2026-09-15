###############################################################################
# FILE NAME: GSC_asa_part2_inference.R
#
# Author:       Aaron Niecestro
# DATE CREATED: February 15, 2026
# DATE EDITED:  February 15, 2026
#
# FILE DESCRIPTION:
#
# This script adapts the ASA Traveling Course:
# "Tree-Based Machine Learning Methods — Part II: Inference and Prediction"
# for use with the German PBC dataset (PBC276). The original workshop examples
# used the Mayo Clinic PBC dataset; this version updates all syntax, data
# preparation, and survival modeling steps to match the structure of the
# PBC276 dataset used in your clinical survival analysis project.
#
# The script demonstrates:
#
#   • Out-of-bag (OOB) inference for survival forests
#   • Extracting key RSF quantities (survival, CHF, mortality)
#   • Plotting inbag vs OOB survival curves
#   • Computing prediction error (C-index, Brier score)
#   • Predicting on new data, including unseen factor levels
#   • Restore mode for performance estimation
#   • Partial dependence plots for survival variables
#
# This file provides a complete inference and prediction workflow for
# Random Forest Survival models parallel to your Cox PH modeling pipeline.
#
# Origin: ASA Traveling Course — Part II: Inference & Prediction
# Adaptation: Aaron Niecestro (German PBC dataset)
###############################################################################

#===============================================================================
# Load required packages
#===============================================================================

library(randomForestSRC)   # core Random Forest Survival package
library(survival)          # Surv() object and survival utilities
library(tidyverse)         # data manipulation and pipes
library(janitor)           # clean_names() for tidy column names

#===============================================================================
# Load and prepare the German PBC dataset (PBC276)
#===============================================================================

PBC276 <- read_csv("~/UTH Survival Analysis/PH1831_Project/PH1831 Project Data/PBC276.csv") %>% 
  clean_names() %>%                     # convert column names to snake_case
  mutate(
    stage = factor(stage),              # convert stage to factor
    drug  = factor(drug),               # convert drug to factor
    sex   = factor(sex)                 # convert sex to factor
  )

glimpse(PBC276)                         # inspect dataset structure

#===============================================================================
# Fit a Random Forest Survival model to PBC276
#===============================================================================

o <- rfsrc(
  Surv(futime, status) ~ .,             # survival formula using all predictors
  data = PBC276                         # German PBC dataset
)

print(o)                                # print RSF summary including OOB error

#===============================================================================
# Key RSF quantities (same as ASA slides)
#===============================================================================

o$time.interest                         # event times used internally by RSF
o$predicted                             # inbag mortality estimates
o$predicted.oob                         # OOB mortality estimates
o$survival                              # inbag survival curves
o$survival.oob                          # OOB survival curves
o$chf                                   # inbag cumulative hazard
o$chf.oob                               # OOB cumulative hazard

#===============================================================================
# Plot inbag vs OOB survival curves for selected patients
#===============================================================================

idx <- c(11, 34, 60)                    # choose three patients to visualize

matplot(
  o$time.interest,                      # x-axis: event times
  t(o$survival[idx, ]),                 # inbag survival curves
  type = "l", col = 4, lwd = 3,
  xlab = "Days", ylab = "Survival"
)

matlines(
  o$time.interest,                      # overlay OOB curves
  t(o$survival.oob[idx, ]),
  type = "l", col = 2, lwd = 3
)

legend(
  "bottomleft",
  legend = c("inbag", "oob"),           # legend labels
  fill = c(4, 2)
)

#===============================================================================
# Inspect RSF object dimensions (same as ASA slide)
#===============================================================================

o$time.interest[1:5]                    # first few event times
o$predicted[1:5]                        # first few inbag mortality estimates
o$predicted.oob[1:5]                    # first few OOB mortality estimates

dim(o$yvar)                             # number of observations
length(o$time.interest)                 # number of event times
dim(o$survival)                         # survival matrix dimensions
dim(o$survival.oob)                     # OOB survival matrix dimensions

#===============================================================================
# Prediction error for survival
#===============================================================================

get.cindex(o)                           # Harrell C-index (discrimination)

get.brier.survival(o)                   # time-varying Brier score

# plot.brier.auc(o)                     # workshop helper (if installed)

#===============================================================================
# Prediction on new data (same as ASA slides)
#===============================================================================

pred10 <- predict(o, PBC276[1:10, ])    # predict survival for first 10 patients
print(pred10)                           # print prediction object

#===============================================================================
# Factor-level robustness example (unseen factor level)
#===============================================================================

test_new <- PBC276[1:3, ]               # take three patients
test_new$drug <- factor(c("newdrug", "D-penicillamine", "Placebo"))
# introduce unseen level "newdrug"

pred_new <- predict(o, test_new)        # predict with unseen level
print(pred_new)                         # RSF handles unseen levels gracefully

print(pred_new$xvar)                    # show how RSF encoded the new level

#===============================================================================
# Restore mode example (same as ASA slides)
#===============================================================================

p_restore <- predict(o, perf.type = "brier")   # restore-mode Brier score
p_restore

#===============================================================================
# Custom estimator example (forest weights)
#===============================================================================

fwt <- predict(o, forest.wt = "oob")$forest.wt   # extract OOB forest weights

mort_hat <- c(fwt %*% o$yvar)                    # custom mortality estimate

summary(mort_hat - o$predicted.oob)              # compare to OOB ensemble

#===============================================================================
# Partial dependence plots for survival (adapted for PBC276)
#===============================================================================

plot.variable(
  o,
  surv.type = "surv",                 # plot survival curves
  xvar.names = "age",                 # variable of interest
  time = 1000,                        # survival time point
  partial = TRUE                      # partial dependence
)

plot.variable(
  o,
  surv.type = "surv",
  xvar.names = "age",
  smooth.lines = TRUE,                # smooth partial dependence
  time = 1000,
  partial = TRUE
)

#===============================================================================
# Partial dependence for log_bili (biochemical marker)
#===============================================================================

partial_obj <- partial(
  o,
  partial.type = "surv",              # survival partial dependence
  partial.xvar = "log_bili",          # variable of interest
  partial.values = quantile(PBC276$log_bili),  # evaluate at quantiles
  partial.time = o$time.interest
)

pdta <- get.partial.plot.data(partial_obj)      # extract plot data

matplot(
  pdta$partial.time,
  t(pdta$yhat),
  type = "l", lty = 1,
  xlab = "Days", ylab = "Adjusted Survival"
)

legend(
  "bottomleft",
  legend = paste0("log_bili = ", quantile(PBC276$log_bili)),
  bty = "n", cex = .75, fill = 1:5
)
