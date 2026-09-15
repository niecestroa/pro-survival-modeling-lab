# Author: Aaron Niecestro
# Created on January 8, 2024
# Last Editted on February 15, 2026

library(tidyverse)
library(survival)
library(asaur)
library(purrr)
library(splines)
library(broom)

#----------------------------------------------------------
# STEP 1: Load & Clean Data
#----------------------------------------------------------

gbcs <- read_csv("~/BIST0627 Applied Survival Data Analysis/BIST0627 Project/project_data_files/gbcs.csv")

gbcs2 <- gbcs %>%
  mutate(
    age          = as.numeric(age),
    menopause_f  = factor(if_else(menopause == "1", "Yes", "No")),
    hormone_f    = factor(if_else(hormone == "1", "Yes", "No")),
    size         = as.numeric(size),
    grade_f      = factor(grade),
    nodes        = as.numeric(nodes),
    prog_recp    = as.numeric(prog_recp),
    estrg_recp   = as.numeric(estrg_recp),
    rectime      = as.numeric(rectime),
    censrec_f    = factor(if_else(censrec == "1", "Recurrence", "Censored")),
    survtime     = as.numeric(survtime),
    censdead     = as.numeric(censdead),
    age.med_f    = factor(if_else(age > 53, "Above(>53)", "Below(<=53)")),
    rectime.med_f = factor(if_else(rectime > 1084, "Above", "Below"))
  ) %>%
  select(
    -id, -diagdateb, -recdate, -deathdate,
    -menopause, -hormone, -grade, -censrec
  )

#----------------------------------------------------------
# STEP 2: Fit Null Cox Model
#----------------------------------------------------------

null_model_gbcs <- gbcs2 %>%
  coxph(Surv(survtime, censdead) ~ 1, data = .)

# Compute martingale residuals
gbcs_resid <- gbcs2 %>%
  mutate(martingale = residuals(null_model_gbcs, type = "martingale"))

#----------------------------------------------------------
# STEP 3: Residual Diagnostics — Numeric Variables
#----------------------------------------------------------

numeric_vars <- c(
  "age", "size", "nodes", "prog_recp",
  "estrg_recp", "rectime", "survtime"
)

numeric_plots <- map(numeric_vars, function(var) {
  gbcs_resid %>%
    ggplot(aes_string(x = var, y = "martingale")) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue", se = FALSE) +
    labs(
      title = paste("Martingale Residuals vs", var),
      x = var,
      y = "Martingale Residuals"
    ) +
    theme_minimal()
})

# Special log-scale diagnostic for nodes
nodes_log_plot <- gbcs_resid %>%
  ggplot(aes(x = nodes, y = martingale)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  scale_x_log10() +
  labs(
    title = "Martingale Residuals vs log(nodes)",
    x = "Nodes (log scale)",
    y = "Martingale Residuals"
  ) +
  theme_minimal()

#----------------------------------------------------------
# STEP 4: Residual Diagnostics — Categorical Variables
#----------------------------------------------------------

categorical_vars <- c(
  "censrec_f", "menopause_f", "hormone_f",
  "grade_f", "censdead", "age.med_f", "rectime.med_f"
)

categorical_plots <- map(categorical_vars, function(var) {
  gbcs_resid %>%
    ggplot(aes_string(x = var, y = "martingale")) +
    geom_boxplot(fill = "lightgray") +
    geom_jitter(width = 0.2, alpha = 0.4) +
    labs(
      title = paste("Martingale Residuals vs", var),
      x = var,
      y = "Martingale Residuals"
    ) +
    theme_minimal()
})

#----------------------------------------------------------
# STEP 5: Fit Multivariable Models (Systematic Comparison)
#----------------------------------------------------------

# Candidate exclusion sets
model_formulas <- list(
  all_vars_1 = as.formula("Surv(survtime, censdead) ~ . - age.med_f - rectime.med_f"),
  all_vars_2 = as.formula("Surv(survtime, censdead) ~ . - age - rectime"),
  all_vars_3 = as.formula("Surv(survtime, censdead) ~ . - age - rectime.med_f"),
  all_vars_4 = as.formula("Surv(survtime, censdead) ~ . - age.med_f - rectime"),
  all_vars_5 = as.formula("Surv(survtime, censdead) ~ . - age.med_f - rectime.med_f")
)

# Fit all models in one pipeline
multi_models <- map(model_formulas, ~ coxph(.x, data = gbcs2))

# Extract significance summaries
multi_summaries <- map(multi_models, broom::tidy)

#----------------------------------------------------------
# STEP 6: Reduced Model Based on Significant Covariates
#----------------------------------------------------------

reduced_formulas <- list(
  reduced_1 = Surv(survtime, censdead) ~ size + prog_recp + rectime,
  reduced_2 = Surv(survtime, censdead) ~ size + prog_recp + rectime.med_f
)

reduced_models <- map(reduced_formulas, ~ coxph(.x, data = gbcs2))
reduced_summaries <- map(reduced_models, broom::tidy)

#----------------------------------------------------------
# STEP 7: 20% Change-in-Estimate Rule
#----------------------------------------------------------

# Function to compute percent change in coefficients
percent_change <- function(full, reduced) {
  full_coef <- coef(full)
  red_coef  <- coef(reduced)
  common    <- intersect(names(full_coef), names(red_coef))
  tibble(
    variable = common,
    full     = full_coef[common],
    reduced  = red_coef[common],
    pct_change = abs((full_coef[common] - red_coef[common]) / full_coef[common]) * 100
  )
}

change_results <- percent_change(
  full  = multi_models$all_vars_1,
  reduced = reduced_models$reduced_1
)

#----------------------------------------------------------
# STEP 8: Re-examine Excluded Variables (Automated)
#----------------------------------------------------------

excluded_vars <- c(
  "age", "nodes", "estrg_recp", "menopause_f",
  "hormone_f", "grade_f", "censrec_f", "age.med_f"
)

# Build formulas programmatically
reexam_formulas <- map(
  excluded_vars,
  ~ as.formula(paste("Surv(survtime, censdead) ~ size + prog_recp + rectime +", .x))
)

reexam_models <- map(reexam_formulas, ~ coxph(.x, data = gbcs2))
reexam_summaries <- map(reexam_models, broom::tidy)

#----------------------------------------------------------
# STEP 9: Transformation Checks (Splines, Log, Quadratic)
#----------------------------------------------------------

# Spline models
spline_models <- list(
  spline_size   = coxph(Surv(survtime, censdead) ~ pspline(size) + prog_recp + rectime, data = gbcs2),
  spline_rectime = coxph(Surv(survtime, censdead) ~ size + prog_recp + pspline(rectime), data = gbcs2)
)

# Individual spline checks
spline_individual <- map(
  c("age", "size", "nodes", "prog_recp", "estrg_recp", "rectime"),
  ~ coxph(as.formula(paste("Surv(survtime, censdead) ~ pspline(", .x, ")")), data = gbcs2)
)

# Log-transform diagnostics (ggplot version)
log_plots <- map(
  c("age", "size", "nodes", "prog_recp", "estrg_recp", "rectime"),
  ~ gbcs_resid %>%
    ggplot(aes_string(x = .x, y = "martingale")) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess", color = "blue", se = FALSE) +
    scale_x_log10() +
    labs(
      title = paste("Martingale Residuals vs log(", .x, ")"),
      x = paste(.x, "(log scale)"),
      y = "Martingale Residuals"
    ) +
    theme_minimal()
)

# Quadratic transformation example
quad_age <- coxph(Surv(survtime, censdead) ~ pspline(age, df = 2), data = gbcs2)

#----------------------------------------------------------
# STEP 10.1: Interaction Screening Pipeline
#----------------------------------------------------------

# Base main-effects model from earlier steps
base_model <- coxph(Surv(survtime, censdead) ~ size + prog_recp + rectime, data = gbcs2)

# Candidate interaction terms to test
interaction_terms <- list(
  size_grade      = "size:grade_f",
  size_menopause  = "size:menopause_f",
  size_hormone    = "size:hormone_f",
  size_censrec    = "size:censrec_f",
  size_lognodes   = "size:log(nodes)",
  menopause_grade = "menopause_f:grade_f",
  menopause_horm  = "menopause_f:hormone_f",
  menopause_cens  = "menopause_f:censrec_f"
)

# Build formulas programmatically
interaction_formulas <- map(
  interaction_terms,
  ~ as.formula(paste(
    "Surv(survtime, censdead) ~ size + prog_recp + rectime +", .x
  ))
)

# Fit all interaction models
interaction_models <- map(interaction_formulas, ~ coxph(.x, data = gbcs2))

# Compare each interaction model to the base model using ANOVA
interaction_anova <- map(
  interaction_models,
  ~ anova(base_model, .x)
)

#----------------------------------------------------------
# STEP 10.2: Full Interaction Model (All 2‑way interactions)
#----------------------------------------------------------

full_interaction_model <- coxph(
  Surv(survtime, censdead) ~ (. - rectime.med_f - age.med_f)^2,
  data = gbcs2
)

full_all_interactions <- coxph(
  Surv(survtime, censdead) ~ (.)^2,
  data = gbcs2
)

#----------------------------------------------------------
# STEP 10.3: Stepwise AIC/BIC Selection
#----------------------------------------------------------

# AIC stepwise
step_aic_inter <- step(full_interaction_model, direction = "both")

# BIC stepwise
step_bic_inter <- step(full_interaction_model, direction = "backward", k = log(nrow(gbcs2)))

#----------------------------------------------------------
# STEP 10.4: Final AIC/BIC Models (Cleaned)
#----------------------------------------------------------

final_aic_model <- step_aic_inter
final_bic_model <- step_bic_inter

# Tidy summaries
aic_summary <- tidy(final_aic_model)
bic_summary <- tidy(final_bic_model)

#----------------------------------------------------------
# STEP 10.5: PH Assumption Checks (Vectorized)
#----------------------------------------------------------

models_to_check <- list(
  final_interaction_1 = final_aic_model,
  final_interaction_2 = final_bic_model,
  base_no_interaction = base_model
)

ph_tests <- map(models_to_check, cox.zph)

#----------------------------------------------------------
# STEP 10.6: Final Candidate Models (Clean Professional Versions)
#----------------------------------------------------------

final_model_1 <- coxph(
  Surv(survtime, censdead) ~ size + prog_recp + rectime.med_f +
    menopause_f + hormone_f + size:hormone_f + prog_recp:rectime.med_f,
  data = gbcs2
)

final_model_2 <- coxph(
  Surv(survtime, censdead) ~ size + prog_recp + rectime +
    prog_recp:rectime,
  data = gbcs2
)

final_model_3 <- coxph(
  Surv(survtime, censdead) ~ size + prog_recp + rectime +
    menopause_f + hormone_f + size:hormone_f + prog_recp:rectime,
  data = gbcs2
)

# PH tests for final models
final_ph <- map(
  list(final_model_1, final_model_2, final_model_3),
  cox.zph
)

#----------------------------------------------------------
# STEP 11: Final Cox Model + PH Assumption + DFbeta Diagnostics
#----------------------------------------------------------

# Final Cox model (from prior selection steps)
fit_final_cox <- coxph(
  Surv(survtime, censdead) ~ size + prog_recp + rectime +
    hormone_f + age.med_f + prog_recp:rectime + size:hormone_f,
  data = gbcs2
)

# PH assumption
ph_final <- cox.zph(fit_final_cox)

# DFbeta residuals
dfbeta_final <- residuals(fit_final_cox, type = "dfbeta") %>%
  as_tibble() %>%
  mutate(obs = row_number())

# Example: DFbeta for a specific coefficient (e.g., 3rd)
dfbeta_plot <- dfbeta_final %>%
  ggplot(aes(x = obs, y = .data[[3]])) +
  geom_segment(aes(xend = obs, yend = 0)) +
  labs(
    title = "DFbeta for Final Cox Model (Coefficient 3)",
    x = "Observation Index",
    y = "DFbeta"
  ) +
  theme_minimal()

#----------------------------------------------------------
# STEP 12: Weibull Model Corresponding to Final Cox Model
#----------------------------------------------------------

fit_final_weib <- survreg(
  Surv(survtime, censdead) ~ size + prog_recp + rectime +
    hormone_f + age.med_f + prog_recp:rectime + size:hormone_f,
  dist = "weibull",
  data = gbcs2
)

#----------------------------------------------------------
# STEP 13: Compare Weibull and Cox Coefficients
#----------------------------------------------------------

weib_coef_all <- fit_final_weib$coefficients[-1]  # drop intercept
weib_scale    <- fit_final_weib$scale
weib_coef_ph  <- -weib_coef_all / weib_scale

cox_coef      <- coef(fit_final_cox)

coef_compare <- tibble(
  term          = names(weib_coef_ph),
  weibull_ph    = as.numeric(weib_coef_ph),
  coxph_coef    = as.numeric(cox_coef[term])
)

#----------------------------------------------------------
# STEP 14: Deviance Residuals for Weibull Model
#----------------------------------------------------------

dev_resid <- residuals(fit_final_weib, type = "deviance")

gbcs_dev <- gbcs2 %>%
  mutate(dev_resid = dev_resid)

dev_plot_hormone <- gbcs_dev %>%
  ggplot(aes(x = hormone_f, y = dev_resid)) +
  geom_boxplot(fill = "lightgray") +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    title = "Deviance Residuals vs Hormone Therapy",
    x = "Hormone Therapy",
    y = "Deviance Residuals"
  ) +
  theme_minimal()

dev_plot_age_med <- gbcs_dev %>%
  ggplot(aes(x = age.med_f, y = dev_resid)) +
  geom_boxplot(fill = "lightgray") +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    title = "Deviance Residuals vs Median Age at Diagnosis",
    x = "Median Age Group",
    y = "Deviance Residuals"
  ) +
  theme_minimal()

#----------------------------------------------------------
# STEP 15: Case-Deletion (DFbeta) for Weibull Model
#----------------------------------------------------------

dfbeta_weib <- residuals(fit_final_weib, type = "dfbeta") %>%
  as_tibble() %>%
  mutate(obs = row_number())

# DFbeta for size (assume column 2)
dfbeta_weib_size_plot <- dfbeta_weib %>%
  ggplot(aes(x = obs, y = .data[[2]])) +
  geom_segment(aes(xend = obs, yend = 0)) +
  labs(
    title = "Case Deletion Plot for Tumor Size (Weibull)",
    x = "Observation Index",
    y = "Change in Size Coefficient"
  ) +
  theme_minimal()

# Flag potential influential observations
influential_size <- dfbeta_weib %>%
  transmute(
    obs,
    dfb_size = .data[[2]]
  ) %>%
  filter(dfb_size > 0.0005 | dfb_size < -0.0005)

#----------------------------------------------------------
# STEP 16: Residual Plots for Slides (Martingale Residuals)
#----------------------------------------------------------

# Assuming `gbcs_resid` already has `martingale` from earlier steps
slide_resid_rectime <- gbcs_resid %>%
  ggplot(aes(x = rectime, y = martingale)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 1.2) +
  labs(
    title = "Martingale Residuals vs Time Until Recurrence",
    x = "Time Until Recurrence",
    y = "Martingale Residuals"
  ) +
  theme_minimal()

slide_resid_prog <- gbcs_resid %>%
  ggplot(aes(x = prog_recp, y = martingale)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", color = "red", se = FALSE, size = 1.2) +
  labs(
    title = "Martingale Residuals vs Progesterone Receptors",
    x = "Number of Progesterone Receptors",
    y = "Martingale Residuals"
  ) +
  theme_minimal()

#----------------------------------------------------------
# STEP 17: Treatment Comparison: Cox vs Weibull (Hormone Only)
#----------------------------------------------------------

# Cox model with hormone only
cox_hormone <- coxph(Surv(survtime, censdead) ~ hormone_f, data = gbcs2)

# Weibull model with hormone only
weib_hormone <- survreg(
  Surv(survtime, censdead) ~ hormone_f,
  dist = "weibull",
  data = gbcs2
)

# Extract Weibull parameters
mu0_hat    <- weib_hormone$coefficients[1]
sigma_hat  <- weib_hormone$scale
alpha_hat  <- 1 / sigma_hat
lambda0_hat <- exp(-mu0_hat)

# Baseline survival function
tt_vec   <- 0:2601
surv0_vec <- 1 - pweibull(tt_vec, shape = alpha_hat, scale = 1 / lambda0_hat)

gamma_hat <- weib_hormone$coefficients[2]
surv1_vec <- surv0_vec^(exp(-gamma_hat / sigma_hat))

# Cox survival estimates
cox_surv_est <- survfit(
  cox_hormone,
  newdata = data.frame(hormone_f = c("No", "Yes"))
)

# Plot option 1: Cox (red) vs Weibull (blue)
plot(cox_surv_est, col = c("red", "red"), lwd = 2,
     ylim = c(0.4, 1),
     xlab = "Survival Time (days)",
     ylab = "Survival Probability",
     main = "Cox and Weibull Hormone Therapy Comparison")

lines(tt_vec, surv0_vec, col = "blue", lwd = 1)
lines(tt_vec, surv1_vec, col = "blue", lwd = 1)

legend(
  "bottomleft",
  legend = c(
    "Cox: Hormone Therapy", "Cox: No Therapy",
    "Weibull: Hormone Therapy", "Weibull: No Therapy"
  ),
  col = c("red", "red", "blue", "blue"),
  lwd = c(2, 2, 1, 1),
  cex = 0.8
)

# Plot option 2: direct comparison by treatment
plot(cox_surv_est, col = c("red", "blue"), lwd = 2,
     xlab = "Survival Time (days)",
     ylab = "Survival Probability",
     main = "Comparisons of Hormone Therapy")

lines(tt_vec, surv0_vec, col = "red", lwd = 1)
lines(tt_vec, surv1_vec, col = "blue", lwd = 1)

legend(
  x = 250, y = 0.4,
  legend = c(
    "Cox: No Therapy", "Weibull: No Therapy",
    "Cox: Hormone Therapy", "Weibull: Hormone Therapy"
  ),
  col = c("red", "red", "blue", "blue"),
  lwd = c(2, 1, 2, 1),
  cex = 0.8
)
