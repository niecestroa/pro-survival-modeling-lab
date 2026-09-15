###############################################################################
# Author:       Aaron Niecestro
# DATE CREATED: February 15, 2026
# DATE EDITED:  February 15, 2026
#
# FILE DESCRIPTION:
#
# This script adapts the ASA Traveling Course:
# "Tree-Based Machine Learning Methods — Part III: Variable Selection"
# for use with the German PBC dataset (PBC276). The original workshop examples
# used mtcars, iris, peakVO2, and high-dimensional microarray data. This version
# rewrites all variable-importance, minimal-depth, VarPro, and ivarpro examples
# so they operate directly on the PBC276 clinical survival dataset.
#
# Demonstrated methods:
#
#   • Permutation variable importance (VIMP)
#   • Block-permutation VIMP
#   • Joint VIMP for variable pairs
#   • Confidence intervals for VIMP via subsampling
#   • Minimal depth variable selection
#   • Guided trees using var.used split counts
#   • VarPro (Variable Priority) for survival
#   • Cross-validated VarPro
#   • ivarpro() — Individual Variable Priority
#   • shap.ivarpro() — SHAP-style variable priority
#   • uvarpro() — Unsupervised variable priority
#
# Origin: ASA Traveling Course — Part III: Variable Selection
# Adaptation: Aaron Niecestro (German PBC dataset)
###############################################################################

#===============================================================================
# Load required packages
#===============================================================================

library(randomForestSRC)   # Random Forest Survival + VIMP + minimal depth
library(survival)          # Surv() object
library(varPro)            # VarPro, ivarpro, uvarpro, subsample
library(tidyverse)         # data manipulation
library(janitor)           # clean_names()

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

#===============================================================================
# Fit a baseline RSF model with permutation VIMP enabled
#===============================================================================

o <- rfsrc(
  Surv(futime, status) ~ .,             # survival formula
  data = PBC276,
  importance = "permute"                # Breiman-Cutler permutation VIMP
)

print(o)                                # print RSF summary

#===============================================================================
# Extract permutation VIMP (same as ASA Slide 7)
#===============================================================================

vimp.grow <- o$importance               # VIMP computed during forest growth

vimp.grow                               # print VIMP values

#===============================================================================
# Block-permutation VIMP (larger blocks → more stable inference)
#===============================================================================

o_block <- rfsrc(
  Surv(futime, status) ~ ., 
  data = PBC276,
  importance = "permute",
  block.size = 10                       # block permutation size
)

vimp.grow.block <- o_block$importance   # block-permutation VIMP

vimp.grow.block                         # print block VIMP

#===============================================================================
# Restore-mode VIMP (same as ASA Slide 7)
#===============================================================================

vimp.restore <- predict(o, importance = "permute")$importance

vimp.restore                            # print restore-mode VIMP

vimp.restore.block <- predict(o, importance = "permute", block.size = 10)$importance

vimp.restore.block                      # print block restore-mode VIMP

#===============================================================================
# Dedicated vimp() interface (same as ASA Slide 7)
#===============================================================================

vimp(o, importance = "permute")$importance

vimp(o, importance = "permute", block.size = 10)$importance

#===============================================================================
# Joint VIMP for pairs of variables (adapted from iris example)
#===============================================================================

# Choose two continuous variables from PBC276
pair1 <- c("age", "albumin")
pair2 <- c("log_bili", "log_copper")

vimp(o, pair1, importance = "permute", joint = TRUE)$importance
vimp(o, pair2, importance = "permute", joint = TRUE)$importance

#===============================================================================
# Confidence intervals for VIMP via subsampling (ASA Slide 11)
#===============================================================================

oo <- subsample(o)                      # subsampling for VIMP CI

o$importance                            # print VIMP

plot.vimp.ci(oo, alpha = .05)           # 95% CI for VIMP

#===============================================================================
# Minimal depth variable selection (ASA Slide 18)
#===============================================================================

md <- max.subtree(o)$order[, 1]         # minimal depth values

barplot(
  sort(md),
  las = 2, horiz = TRUE, col = "cadetblue3",
  main = "Minimal Depth Variable Ranking (PBC276)"
)

#===============================================================================
# Guided trees using var.used split counts (ASA Slide 19)
#===============================================================================

xvar.used <- predict(o, var.used = "all.trees")$var.used

o_guided <- rfsrc(
  Surv(futime, status) ~ ., 
  data = PBC276,
  xvar.wt = xvar.used                   # weight variables by split frequency
)

md_guided <- max.subtree(o_guided)$order[, 1]

barplot(
  sort(md_guided),
  las = 2, horiz = TRUE, col = "cadetblue3",
  main = "Guided Minimal Depth (PBC276)"
)

#===============================================================================
# VarPro canonical illustration (ASA Slide 26)
#===============================================================================

vp <- varpro(
  Surv(futime, status) ~ ., 
  PBC276,
  ntree = 300                           # more stable VarPro
)

importance(vp)                          # VarPro importance

vp_cv <- cv.varpro(
  Surv(futime, status) ~ ., 
  PBC276
)

print(vp_cv)                            # cross-validated VarPro

#===============================================================================
# Examine cross-validated VarPro results (ASA Slide 27)
#===============================================================================

vp_cv$imp                               # raw importance
vp_cv$imp.conserve                      # conservative selection
vp_cv$imp.liberal                       # liberal selection
vp_cv$err                               # CV error
vp_cv$zcut                              # cutoff threshold

#===============================================================================
# VarPro + VIMP CI (ASA Slide 28)
#===============================================================================

vp_rfsrc <- rfsrc(
  Surv(futime, status) ~ ., 
  PBC276,
  importance = "permute"
)

vp_sub <- subsample(vp_rfsrc)

plot.vimp.ci(vp_sub, alpha = .05)

barplot(
  vp_cv$imp.liberal$z,
  names.arg = vp_cv$imp.liberal$variable,
  las = 2, horiz = TRUE, col = "coral2",
  main = "Liberal VarPro Importance (PBC276)"
)

#===============================================================================
# ivarpro() — Individual Variable Priority (ASA Slide 34)
#===============================================================================

ivp <- ivarpro(vp)                      # individual variable priority

print(ivp[1:5, 1:8])                    # print first few rows

imp_cut <- ivarpro(vp, cut.max = 2, adaptive = FALSE)

shap.ivarpro(imp_cut)                   # SHAP-style variable priority

#===============================================================================
# plot.ivarpro() (ASA Slide 35)
#===============================================================================

plot(
  imp_cut,
  var = "log_bili",                     # variable to visualize
  col.var = "age",                      # color by age
  size.var = "y"                        # size by response
)

#===============================================================================
# uvarpro() — Unsupervised Variable Priority (ASA Slide 38)
#===============================================================================

uvp <- uvarpro(PBC276)                  # unsupervised variable priority

head(importance(uvp))                   # print top unsupervised variables

beta <- get.beta.entropy(uvp)           # entropy-based dependency matrix

beta[1:4, 1:4]                          # print first block

sdependent(beta)                        # structural dependency summary

###############################################################################
# End of ASA Part III — PBC276 Adaptation
###############################################################################
