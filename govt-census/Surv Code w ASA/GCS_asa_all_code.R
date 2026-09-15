###############################################################################
# Author:       Aaron Niecestro
# DATE CREATED: February 15, 2026
# DATE EDITED:  February 15, 2026
#
# FILE DESCRIPTION:
#
# Master script that orchestrates:
#   • Part I  – Training (RSF on PBC276)
#   • Part II – Inference & Prediction
#   • Part III – Variable Selection
#   • Part IV – Advanced Topics (imputation, OOD)
#
# All parts are adapted from the ASA Traveling Course:
# "Tree-Based Machine Learning Methods" for the German PBC dataset (PBC276).
###############################################################################

#===============================================================================
# 0. Global setup
#===============================================================================

set.seed(1831)                      # global reproducibility

library(tidyverse)
library(janitor)
library(survival)
library(randomForestSRC)
library(varPro)

#===============================================================================
# 1. Load and prepare PBC276 once
#===============================================================================

PBC276 <- read_csv("~/UTH Survival Analysis/PH1831_Project/PH1831 Project Data/PBC276.csv") %>%
  clean_names() %>%
  mutate(
    stage = factor(stage),
    drug  = factor(drug),
    sex   = factor(sex)
  )

glimpse(PBC276)

#===============================================================================
# 2. Part I – Training (RSF baseline model)
#===============================================================================

run_part1_training <- function(data) {
  fit_rsf <- rfsrc(
    Surv(futime, status) ~ .,
    data = data,
    ntree = 1000,
    splitrule = "logrank",
    importance = "permute"
  )
  print(fit_rsf)
  fit_rsf
}

#===============================================================================
# 3. Part II – Inference & Prediction
#===============================================================================

run_part2_inference <- function(data, fit_rsf) {
  # OOB survival curves for selected cases
  idx <- c(11, 34, 60)
  matplot(
    fit_rsf$time.interest,
    t(fit_rsf$survival[idx, ]),
    type = "l", col = 4, lwd = 3,
    xlab = "Days", ylab = "Survival"
  )
  matlines(
    fit_rsf$time.interest,
    t(fit_rsf$survival.oob[idx, ]),
    type = "l", col = 2, lwd = 3
  )
  legend("bottomleft", legend = c("inbag", "oob"), fill = c(4, 2))
  
  # C-index and Brier
  print(get.cindex(fit_rsf))
  print(get.brier.survival(fit_rsf))
  
  # Prediction on first 10 patients
  pred10 <- predict(fit_rsf, data[1:10, ])
  print(pred10)
  
  # Unseen factor level example
  test_new <- data[1:3, ]
  test_new$drug <- factor(c("newdrug", as.character(test_new$drug[2:3])))
  pred_new <- predict(fit_rsf, test_new)
  print(pred_new$xvar)
}

#===============================================================================
# 4. Part III – Variable Selection (VIMP, VarPro, minimal depth)
#===============================================================================

run_part3_varselect <- function(data, fit_rsf) {
  # VIMP
  print(fit_rsf$importance)
  
  # Subsample CI for VIMP
  oo <- subsample(fit_rsf)
  plot.vimp.ci(oo, alpha = .05)
  
  # Minimal depth
  md <- max.subtree(fit_rsf)$order[, 1]
  barplot(sort(md), las = 2, horiz = TRUE, col = "cadetblue3",
          main = "Minimal Depth (PBC276)")
  
  # VarPro
  vp <- varpro(Surv(futime, status) ~ ., data, ntree = 300)
  print(importance(vp))
  
  vp_cv <- cv.varpro(Surv(futime, status) ~ ., data)
  print(vp_cv)
  
  # ivarpro
  ivp <- ivarpro(vp)
  print(ivp[1:5, 1:8])
}

#===============================================================================
# 5. Part IV – Advanced Topics (Imputation, OOD)
#===============================================================================

run_part4_advanced <- function(data) {
  # Train/test split
  id <- sample(seq_len(nrow(data)), 100)
  train <- data[id, ]
  test  <- data[-id, ]
  
  # Supervised imputation
  fit_imp <- impute.learn(
    data = train,
    mf.q = 1,
    max.iter = 5,
    full.sweep.options = list(ntree = 25, nsplit = 5),
    target.mode = "all",
    supervised.formula = Surv(futime, status) ~ .,
    supervised.args = list(ntree = 50, nsplit = 5),
    save.ood = TRUE
  )
  
  # Test-time imputation + OOD
  test_imp <- predict(fit_imp, test, max.predict.iter = 2)
  ood <- impute.ood(fit_imp, test)
  
  print(head(ood$score))
  print(head(ood$score.percentile))
}

#===============================================================================
# 6. Driver
#===============================================================================

main <- function() {
  fit_rsf <- run_part1_training(PBC276)
  run_part2_inference(PBC276, fit_rsf)
  run_part3_varselect(PBC276, fit_rsf)
  run_part4_advanced(PBC276)
}

if (sys.nframe() == 0) main()
