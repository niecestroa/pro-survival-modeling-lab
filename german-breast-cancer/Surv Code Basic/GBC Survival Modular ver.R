############################################################
# GBCS Survival Analysis Pipeline — Full Modular Workflow
# Author:       Aaron Niecestro
# Date Created: September 15, 2026
# Last Edit:    September 15, 2026
# Description:
#
# This script implements a fully modular, production‑grade
# survival‑analysis pipeline for the GBCS breast‑cancer dataset.
# Each module is designed to be reusable, interpretable, and
# statistically rigorous. The workflow covers every major step
# required in modern clinical survival modeling.
#
# MODULE 0 — Utility Infrastructure
# Defines helper functions for formula construction, Cox/Weibull
# model fitting, residual diagnostics, DFbeta influence analysis,
# spline/log transformations, and interaction formula generation.
# These utilities eliminate duplication and ensure consistent,
# maintainable modeling throughout the pipeline.
#
# MODULE 1 — Data Loading & Cleaning
# Loads raw GBCS data, converts numeric and factor variables,
# creates median‑split versions of age and rectime, and removes
# unused ID/date fields. This ensures the dataset is clean,
# consistent, and ready for survival modeling.
#
# MODULE 2 — Null Cox Model & Martingale Residuals
# Fits a Cox model with no covariates and computes martingale
# residuals. These residuals are essential for diagnosing
# functional forms and identifying nonlinear relationships in
# continuous predictors.
#
# MODULE 3 — Residual Diagnostics
# Produces scatter plots (numeric predictors), boxplots
# (categorical predictors), log‑scale diagnostics, and special
# log(nodes) plots. These diagnostics reveal whether predictors
# require log transforms, splines, or categorical treatment.
#
# MODULE 4 — Multivariable Model Comparison
# Fits multiple Cox models with different exclusion sets to test
# sensitivity to age vs age.med_f and rectime vs rectime.med_f.
# This identifies robust predictors and highlights unstable or
# redundant variables.
#
# MODULE 5 — Reduced Models
# Fits reduced Cox models containing only size, progesterone
# receptors, and rectime (or rectime.med_f). These serve as base
# models for confounding assessment, interaction screening, and
# transformation testing.
#
# MODULE 6 — 20% Change‑in‑Estimate Rule
# Computes percent change in coefficients between full and reduced
# models. This identifies confounders—variables that meaningfully
# alter the effect size of core predictors.
#
# MODULE 7 — Re‑examining Excluded Variables
# Adds each excluded variable back into the reduced model and
# refits. This ensures no important predictor was mistakenly
# removed and evaluates each variable’s independent contribution.
#
# MODULE 8 — Transformation Checks
# Fits spline models, log‑scale diagnostics, and quadratic
# approximations. This determines whether predictors behave
# linearly, log‑linearly, or require spline modeling.
#
# MODULE 9 — Interaction Screening
# Tests biologically plausible interactions using ANOVA
# comparisons against a base model. Identifies effect modification
# such as hormone therapy modifying the effect of tumor size.
#
# MODULE 10 — Full Interaction Models + Stepwise AIC/BIC
# Fits full two‑way interaction models and performs AIC and BIC
# stepwise selection. This identifies the best interaction
# structure under different penalization regimes.
#
# MODULE 11 — Final Candidate Cox Models
# Fits curated final models based on biological plausibility,
# statistical significance, PH stability, and interpretability.
# These models represent clinically meaningful candidates.
#
# MODULE 12 — Final Cox Model + PH + DFbeta
# Fits the final chosen Cox model, performs PH assumption testing,
# and computes DFbeta influence diagnostics. Ensures the final
# model is valid, stable, and not driven by influential outliers.
#
# MODULE 13 — Weibull Model + Comparison
# Fits a Weibull parametric model with the same covariates,
# converts coefficients to PH‑equivalent form, and compares them
# to Cox coefficients. Tests whether the Weibull assumption is
# reasonable and whether parametric and semiparametric results
# align.
#
# MODULE 14 — Deviance Residuals + Case Deletion
# Computes Weibull deviance residuals, plots them against hormone
# therapy and age groups, computes DFbeta for the Weibull model,
# and identifies influential cases. Complements Cox diagnostics.
#
# MODULE 15 — Slide‑Ready Residual Plots
# Generates polished martingale residual plots optimized for
# presentation slides, highlighting rectime and progesterone
# receptor effects.
#
# MODULE 16 — Hormone Treatment Comparison
# Fits Cox and Weibull models with hormone therapy only, computes
# survival curves for both models, and overlays them. Provides a
# treatment‑effect visualization under both modeling frameworks.
#
# MODULE 17 — Full Pipeline Orchestrator
# Runs every module in sequence and returns a structured list
# containing all models, diagnostics, plots, and results. Enables
# reproducibility, automation, and integration into reports or
# dashboards.
#
############################################################


#----------------------------------------------------------
# 0. Package loading and global options
#----------------------------------------------------------

suppressPackageStartupMessages({                     # Suppress startup messages to keep console output clean and professional.
  library(tidyverse)                                 # Load tidyverse for data manipulation, piping, and plotting; core for modern R workflows.
  library(survival)                                  # Load survival for Cox models, parametric survival models, and survival objects.
  library(asaur)                                     # Load asaur for applied survival analysis utilities; supports teaching-style datasets and tools.
  library(purrr)                                     # Load purrr for functional programming; used for mapping over lists and formulas.
  library(splines)                                   # Load splines for pspline functions; enables flexible modeling of nonlinear effects.
  library(broom)                                     # Load broom for tidying model outputs; converts model objects to clean data frames.
})

set.seed(2026)                                       # Set a global seed for reproducibility; ensures consistent results across runs.

#----------------------------------------------------------
# 1. Utility helper functions
#----------------------------------------------------------

safe_coxph   <- safely(coxph)                        # Wrap coxph in safely; returns list with result or error, preventing pipeline crashes.
safe_survreg <- safely(survreg)                      # Wrap survreg in safely; similarly protects parametric survival fits from stopping execution.

make_surv_formula <- function(rhs) {                 # Define helper to build a Surv(...) ~ RHS formula from a right-hand side string.
  as.formula(paste("Surv(survtime, censdead) ~", rhs)) # Construct full formula string and convert to formula; standardizes survival model specification.
}

make_rhs <- function(terms) {                        # Define helper to collapse a vector of terms into a single RHS string.
  paste(terms, collapse = " + ")                     # Join terms with '+'; ensures consistent formula construction for additive models.
}

fit_cox <- function(formula, data) {                 # Define a general Cox fitting wrapper that returns model, tidy output, PH test, and DFbeta.
  m <- coxph(formula, data = data)                   # Fit Cox proportional hazards model using provided formula and dataset.
  list(
    model   = m,                                     # Store the fitted Cox model object for later use in predictions and diagnostics.
    tidy    = broom::tidy(m),                        # Tidy the model coefficients into a data frame; simplifies reporting and inspection.
    ph      = cox.zph(m),                            # Compute proportional hazards test using Schoenfeld residuals; checks PH assumption.
    dfbeta  = residuals(m, type = "dfbeta") |>       # Compute DFbeta residuals to assess influence of individual observations on coefficients.
      as_tibble() |>                                 # Convert DFbeta matrix to tibble for easier manipulation and plotting.
      mutate(obs = row_number())                     # Add observation index column; useful for plotting and identifying influential cases.
  )
}

fit_weibull <- function(formula, data, dist = "weibull") { # Define a general Weibull (or other dist) survival regression wrapper.
  m <- survreg(formula, dist = dist, data = data)    # Fit parametric survival model with specified distribution and formula.
  list(
    model   = m,                                     # Store the fitted parametric survival model object.
    tidy    = broom::tidy(m),                        # Tidy the parametric model coefficients into a data frame for interpretation.
    dfbeta  = residuals(m, type = "dfbeta") |>       # Compute DFbeta residuals for parametric model; assesses influence on coefficients.
      as_tibble() |>                                 # Convert DFbeta matrix to tibble for easier downstream processing.
      mutate(obs = row_number())                     # Add observation index for plotting and case identification.
  )
}

plot_martingale_scatter <- function(data, var) {     # Define helper to create scatter plot of martingale residuals vs a numeric predictor.
  ggplot(data, aes_string(x = var, y = "martingale")) + # Map predictor to x and martingale residuals to y; uses aes_string for dynamic variable.
    geom_point(alpha = 0.5) +                        # Add semi-transparent points to visualize residual distribution without overplotting.
    geom_smooth(method = "loess", color = "blue", se = FALSE) + # Add LOESS smooth to detect nonlinear patterns; no confidence band for clarity.
    labs(
      title = paste("Martingale Residuals vs", var),  # Set informative title indicating predictor being assessed.
      x = var,                                       # Label x-axis with predictor name for interpretability.
      y = "Martingale Residuals"                     # Label y-axis to indicate residual type; important for diagnostic context.
    ) +
    theme_minimal()                                  # Use minimal theme for clean, presentation-ready plots.
}

plot_martingale_box <- function(data, var) {         # Define helper to create boxplot of martingale residuals vs a categorical predictor.
  ggplot(data, aes_string(x = var, y = "martingale")) + # Map categorical predictor to x and martingale residuals to y.
    geom_boxplot(fill = "lightgray") +               # Draw boxplots to summarize residual distribution across categories.
    geom_jitter(width = 0.2, alpha = 0.4) +          # Add jittered points to show individual residuals; avoids overlap and reveals outliers.
    labs(
      title = paste("Martingale Residuals vs", var),  # Title indicates categorical variable being examined.
      x = var,                                       # X-axis label is the categorical predictor name.
      y = "Martingale Residuals"                     # Y-axis label clarifies residual type for diagnostic interpretation.
    ) +
    theme_minimal()                                  # Minimal theme for clarity and consistency across plots.
}

plot_dfbeta_term <- function(dfbeta_tbl, term_index, title_suffix = "") { # Define helper to plot DFbeta for a specific coefficient index.
  ggplot(dfbeta_tbl, aes(x = obs, y = .data[[term_index]])) + # Map observation index to x and DFbeta for chosen term to y.
    geom_segment(aes(xend = obs, yend = 0)) +         # Draw vertical segments from each DFbeta value to zero; highlights magnitude and direction.
    labs(
      title = paste("DFbeta for Coefficient", term_index, title_suffix), # Title indicates which coefficient’s influence is being visualized.
      x = "Observation Index",                        # X-axis label clarifies that each point corresponds to an observation.
      y = "DFbeta"                                    # Y-axis label indicates DFbeta values, representing change in coefficient if case removed.
    ) +
    theme_minimal()                                  # Minimal theme for clean diagnostic visualization.
}

percent_change <- function(full, reduced) {          # Define helper to compute percent change in coefficients between full and reduced models.
  full_coef <- coef(full)                            # Extract coefficient vector from full model; serves as baseline.
  red_coef  <- coef(reduced)                         # Extract coefficient vector from reduced model; used for comparison.
  common    <- intersect(names(full_coef), names(red_coef)) # Identify coefficients present in both models; ensures valid comparison.
  tibble(
    variable   = common,                             # Store variable names for which percent change is computed.
    full       = full_coef[common],                  # Store full model coefficients for these variables.
    reduced    = red_coef[common],                   # Store reduced model coefficients for these variables.
    pct_change = abs((full_coef[common] - red_coef[common]) / full_coef[common]) * 100 # Compute absolute percent change; used for confounding assessment.
  )
}

build_interaction_formula <- function(main_terms, interaction_term) { # Define helper to build formula with main effects plus one interaction term.
  rhs <- make_rhs(c(main_terms, interaction_term))   # Combine main terms and interaction term into a single RHS string.
  make_surv_formula(rhs)                             # Convert RHS string into a Surv(...) ~ RHS formula for Cox modeling.
}

build_additive_formula <- function(main_terms) {     # Define helper to build formula with only additive main effects.
  rhs <- make_rhs(main_terms)                        # Collapse main terms into RHS string with '+' separators.
  make_surv_formula(rhs)                             # Convert RHS string into a Surv(...) ~ RHS formula.
}

log_scale_plot <- function(data, var) {              # Define helper to plot martingale residuals vs predictor on log scale.
  ggplot(data, aes_string(x = var, y = "martingale")) + # Map predictor to x and martingale residuals to y.
    geom_point(alpha = 0.5) +                        # Add semi-transparent points to show residual distribution.
    geom_smooth(method = "loess", color = "blue", se = FALSE) + # Add LOESS smooth to detect nonlinear patterns on log scale.
    scale_x_log10() +                                 # Transform x-axis to log10 scale; useful for skewed predictors like counts.
    labs(
      title = paste("Martingale Residuals vs log(", var, ")"), # Title indicates log-transformed predictor being assessed.
      x = paste(var, "(log scale)"),                 # X-axis label clarifies that predictor is shown on log scale.
      y = "Martingale Residuals"                     # Y-axis label indicates residual type for diagnostic context.
    ) +
    theme_minimal()                                  # Minimal theme for clean, interpretable plots.
}

#----------------------------------------------------------
# 2. Load & clean data
#----------------------------------------------------------

load_and_clean_gbcs <- function(path) {              # Define function to load and preprocess GBCS dataset from a given file path.
  read_csv(path) %>%                                 # Read CSV file into a tibble; assumes standard comma-separated format.
    mutate(
      age          = as.numeric(age),                # Convert age to numeric; ensures proper handling in models and avoids factor issues.
      menopause_f  = factor(if_else(menopause == "1", "Yes", "No")), # Create factor for menopause status; improves interpretability in models.
      hormone_f    = factor(if_else(hormone == "1", "Yes", "No")),   # Create factor for hormone therapy; used as categorical treatment variable.
      size         = as.numeric(size),               # Convert tumor size to numeric; required for continuous modeling and transformations.
      grade_f      = factor(grade),                  # Convert tumor grade to factor; categorical representation for ordinal pathology.
      nodes        = as.numeric(nodes),              # Convert number of nodes to numeric; important for continuous or log-transformed modeling.
      prog_recp    = as.numeric(prog_recp),          # Convert progesterone receptor count to numeric; used as continuous biomarker.
      estrg_recp   = as.numeric(estrg_recp),         # Convert estrogen receptor count to numeric; similar continuous biomarker.
      rectime      = as.numeric(rectime),            # Convert recurrence time to numeric; key time-to-event predictor.
      censrec_f    = factor(if_else(censrec == "1", "Recurrence", "Censored")), # Create factor indicating recurrence vs censored; useful for descriptive plots.
      survtime     = as.numeric(survtime),           # Convert survival time to numeric; required for Surv object construction.
      censdead     = as.numeric(censdead),           # Convert censoring indicator for death to numeric; used as event indicator in Surv.
      age.med_f    = factor(if_else(age > 53, "Above(>53)", "Below(<=53)")), # Create median-based age group; supports categorical age analysis.
      rectime.med_f = factor(if_else(rectime > 1084, "Above", "Below")) # Create median-based recurrence time group; used for categorical modeling.
    ) %>%
    select(
      -id, -diagdateb, -recdate, -deathdate,         # Drop ID and date variables; not needed for modeling and may complicate formulas.
      -menopause, -hormone, -grade, -censrec         # Drop original coding variables now replaced by factor versions; avoids redundancy.
    )
}

#----------------------------------------------------------
# 3. Null Cox model + martingale residuals
#----------------------------------------------------------

compute_null_martingale <- function(data) {          # Define function to fit null Cox model and compute martingale residuals.
  null_model <- coxph(Surv(survtime, censdead) ~ 1, data = data) # Fit Cox model with no covariates; baseline hazard only.
  data %>%
    mutate(martingale = residuals(null_model, type = "martingale")) # Add martingale residuals to dataset; used for functional form diagnostics.
}

#----------------------------------------------------------
# 4. Residual diagnostics (numeric & categorical)
#----------------------------------------------------------

diagnostic_plots <- function(resid_data) {           # Define function to generate residual diagnostic plots for numeric and categorical variables.
  numeric_vars <- c("age", "size", "nodes", "prog_recp", # List numeric predictors to examine via scatter plots.
                    "estrg_recp", "rectime", "survtime")
  
  categorical_vars <- c(
    "censrec_f", "menopause_f", "hormone_f",         # List categorical predictors to examine via boxplots.
    "grade_f", "censdead", "age.med_f", "rectime.med_f"
  )
  
  numeric_plots <- map(numeric_vars, ~ plot_martingale_scatter(resid_data, .x)) # Generate scatter plots for each numeric variable.
  names(numeric_plots) <- numeric_vars           # Name list elements by variable for easier access and interpretation.
  
  categorical_plots <- map(categorical_vars, ~ plot_martingale_box(resid_data, .x)) # Generate boxplots for each categorical variable.
  names(categorical_plots) <- categorical_vars   # Name list elements by variable for clarity.
  
  nodes_log_plot <- resid_data %>%               # Create special log-scale plot for nodes; often skewed and benefits from log transform.
    ggplot(aes(x = nodes, y = martingale)) +
    geom_point(alpha = 0.5) +                    # Plot residuals vs nodes with semi-transparent points.
    geom_smooth(method = "loess", color = "red", se = FALSE) + # Add LOESS smooth to detect nonlinear patterns in nodes.
    scale_x_log10() +                            # Use log10 scale for nodes; helps reveal structure in skewed distributions.
    labs(
      title = "Martingale Residuals vs log(nodes)", # Title indicates log-transformed nodes being assessed.
      x = "Nodes (log scale)",                   # X-axis label clarifies log scaling.
      y = "Martingale Residuals"                 # Y-axis label indicates residual type.
    ) +
    theme_minimal()                              # Minimal theme for clean visualization.
  
  list(
    numeric      = numeric_plots,                # Return list of numeric residual plots for further inspection or saving.
    categorical  = categorical_plots,            # Return list of categorical residual plots.
    nodes_log    = nodes_log_plot                # Return special nodes log-scale plot.
  )
}

#----------------------------------------------------------
# 5. Multivariable models (systematic comparison)
#----------------------------------------------------------

fit_multivariable_models <- function(data) {         # Define function to fit multiple Cox models with different exclusion sets.
  model_formulas <- list(
    all_vars_1 = as.formula("Surv(survtime, censdead) ~ . - age.med_f - rectime.med_f"), # Model excluding median-based age and rectime groups.
    all_vars_2 = as.formula("Surv(survtime, censdead) ~ . - age - rectime"),             # Model excluding continuous age and rectime.
    all_vars_3 = as.formula("Surv(survtime, censdead) ~ . - age - rectime.med_f"),       # Model excluding continuous age and median rectime group.
    all_vars_4 = as.formula("Surv(survtime, censdead) ~ . - age.med_f - rectime"),       # Model excluding median age group and continuous rectime.
    all_vars_5 = as.formula("Surv(survtime, censdead) ~ . - age.med_f - rectime.med_f")  # Model excluding both median-based age and rectime groups.
  )
  
  multi_models <- map(model_formulas, ~ coxph(.x, data = data)) # Fit Cox models for each formula; systematic comparison of variable sets.
  multi_summaries <- map(multi_models, broom::tidy)             # Tidy each model’s coefficients; facilitates comparison of significance patterns.
  
  list(
    formulas  = model_formulas,                 # Return list of formulas used; useful for documentation and reproducibility.
    models    = multi_models,                   # Return list of fitted Cox models for further diagnostics.
    summaries = multi_summaries                 # Return list of tidy summaries for reporting and variable selection.
  )
}

#----------------------------------------------------------
# 6. Reduced models based on significant covariates
#----------------------------------------------------------

fit_reduced_models <- function(data) {               # Define function to fit reduced Cox models based on key significant covariates.
  reduced_formulas <- list(
    reduced_1 = Surv(survtime, censdead) ~ size + prog_recp + rectime,      # Reduced model with size, progesterone receptors, and rectime.
    reduced_2 = Surv(survtime, censdead) ~ size + prog_recp + rectime.med_f # Alternative reduced model using median-based rectime group.
  )
  
  reduced_models    <- map(reduced_formulas, ~ coxph(.x, data = data))      # Fit Cox models for each reduced formula.
  reduced_summaries <- map(reduced_models, broom::tidy)                     # Tidy reduced model coefficients for interpretation and comparison.
  
  list(
    formulas  = reduced_formulas,               # Return reduced formulas for documentation.
    models    = reduced_models,                 # Return fitted reduced models for diagnostics and PH checks.
    summaries = reduced_summaries               # Return tidy summaries for reporting and confounding assessment.
  )
}

#----------------------------------------------------------
# 7. 20% change-in-estimate rule
#----------------------------------------------------------

compute_change_in_estimate <- function(full_model, reduced_model) { # Define function to compute percent change in coefficients between full and reduced models.
  percent_change(full_model, reduced_model)      # Use previously defined helper to calculate percent change; supports confounding evaluation.
}

#----------------------------------------------------------
# 8. Re-examine excluded variables
#----------------------------------------------------------

reexamine_excluded <- function(data) {               # Define function to systematically reintroduce excluded variables into reduced model.
  excluded_vars <- c(
    "age", "nodes", "estrg_recp", "menopause_f",    # List variables previously excluded; candidates for re-examination.
    "hormone_f", "grade_f", "censrec_f", "age.med_f"
  )
  
  base_terms <- c("size", "prog_recp", "rectime")   # Define base model terms that remain in all re-examination models.
  
  reexam_formulas <- map(
    excluded_vars,
    ~ build_additive_formula(c(base_terms, .x))     # Build formulas adding one excluded variable at a time to base terms.
  )
  
  reexam_models    <- map(reexam_formulas, ~ coxph(.x, data = data)) # Fit Cox models for each re-examination formula.
  reexam_summaries <- map(reexam_models, broom::tidy)                # Tidy each model’s coefficients; assess significance of reintroduced variables.
  
  list(
    vars      = excluded_vars,                    # Return list of excluded variables considered.
    formulas  = reexam_formulas,                  # Return corresponding formulas for transparency.
    models    = reexam_models,                    # Return fitted models for diagnostics and PH checks.
    summaries = reexam_summaries                  # Return tidy summaries for interpretation and reporting.
  )
}

#----------------------------------------------------------
# 9. Transformation checks (splines, log, quadratic)
#----------------------------------------------------------

fit_spline_models <- function(data) {                # Define function to fit spline-based Cox models for key predictors.
  spline_models <- list(
    spline_size    = coxph(Surv(survtime, censdead) ~ pspline(size) + prog_recp + rectime, data = data), # Model with spline for size plus main effects.
    spline_rectime = coxph(Surv(survtime, censdead) ~ size + prog_recp + pspline(rectime), data = data)  # Model with spline for rectime plus main effects.
  )
  
  spline_individual <- map(
    c("age", "size", "nodes", "prog_recp", "estrg_recp", "rectime"),
    ~ coxph(as.formula(paste("Surv(survtime, censdead) ~ pspline(", .x, ")")), data = data) # Fit single-predictor spline models for each numeric variable.
  )
  
  quad_age <- coxph(Surv(survtime, censdead) ~ pspline(age, df = 2), data = data) # Fit spline with df=2 for age; approximates quadratic effect.
  
  list(
    spline_joint      = spline_models,             # Return joint spline models for size and rectime.
    spline_individual = spline_individual,         # Return individual spline models for each numeric predictor.
    quad_age          = quad_age                   # Return age spline model approximating quadratic relationship.
  )
}

log_diagnostic_plots <- function(resid_data) {       # Define function to generate log-scale residual plots for numeric predictors.
  vars <- c("age", "size", "nodes", "prog_recp", "estrg_recp", "rectime") # List numeric variables to examine on log scale.
  log_plots <- map(vars, ~ log_scale_plot(resid_data, .x))                # Generate log-scale plots for each variable.
  names(log_plots) <- vars                     # Name list elements by variable for easy access.
  log_plots                                   # Return list of log-scale diagnostic plots.
}

#----------------------------------------------------------
# 10. Interaction screening pipeline
#----------------------------------------------------------

interaction_screening <- function(data) {            # Define function to screen candidate interaction terms using ANOVA comparisons.
  base_model <- coxph(Surv(survtime, censdead) ~ size + prog_recp + rectime, data = data) # Define base main-effects Cox model.
  
  interaction_terms <- list(
    size_grade      = "size:grade_f",               # Candidate interaction between size and tumor grade.
    size_menopause  = "size:menopause_f",           # Candidate interaction between size and menopause status.
    size_hormone    = "size:hormone_f",             # Candidate interaction between size and hormone therapy.
    size_censrec    = "size:censrec_f",             # Candidate interaction between size and recurrence status.
    size_lognodes   = "size:log(nodes)",            # Candidate interaction between size and log-transformed nodes.
    menopause_grade = "menopause_f:grade_f",        # Candidate interaction between menopause status and tumor grade.
    menopause_horm  = "menopause_f:hormone_f",      # Candidate interaction between menopause status and hormone therapy.
    menopause_cens  = "menopause_f:censrec_f"       # Candidate interaction between menopause status and recurrence status.
  )
  
  main_terms <- c("size", "prog_recp", "rectime")   # Define main effects to include in all interaction models.
  
  interaction_formulas <- map(
    interaction_terms,
    ~ build_interaction_formula(main_terms, .x)     # Build formulas adding each interaction term to main effects.
  )
  
  interaction_models <- map(interaction_formulas, ~ coxph(.x, data = data)) # Fit Cox models including each interaction term.
  
  interaction_anova <- map(
    interaction_models,
    ~ anova(base_model, .x)                         # Compare each interaction model to base model via likelihood ratio test.
  )
  
  list(
    base_model          = base_model,               # Return base main-effects model.
    interaction_terms   = interaction_terms,        # Return list of interaction term specifications.
    interaction_formulas = interaction_formulas,    # Return formulas used for interaction models.
    interaction_models  = interaction_models,       # Return fitted interaction models.
    interaction_anova   = interaction_anova         # Return ANOVA comparisons for significance assessment.
  )
}

#----------------------------------------------------------
# 11. Full interaction models + stepwise selection
#----------------------------------------------------------

fit_full_interaction_models <- function(data) {      # Define function to fit full interaction models and perform stepwise selection.
  full_interaction_model <- coxph(
    Surv(survtime, censdead) ~ (. - rectime.med_f - age.med_f)^2,
    data = data
  )                                                 # Fit model with all two-way interactions among covariates except median-based age and rectime groups.
  
  full_all_interactions <- coxph(
    Surv(survtime, censdead) ~ (.)^2,
    data = data
  )                                                 # Fit model with all two-way interactions among all covariates; most saturated interaction structure.
  
  step_aic_inter <- step(full_interaction_model, direction = "both") # Perform stepwise selection using AIC; allows both forward and backward moves.
  step_bic_inter <- step(full_interaction_model, direction = "backward", k = log(nrow(data))) # Perform backward selection using BIC; penalizes complexity more strongly.
  
  list(
    full_interaction      = full_interaction_model, # Return full interaction model excluding median-based groups.
    full_all_interactions = full_all_interactions,  # Return fully saturated interaction model.
    step_aic              = step_aic_inter,         # Return AIC-selected interaction model.
    step_bic              = step_bic_inter          # Return BIC-selected interaction model.
  )
}

ph_assumption_checks <- function(models_list) {      # Define function to compute PH assumption tests for a list of models.
  map(models_list, cox.zph)                         # Apply cox.zph to each model; returns Schoenfeld-based PH diagnostics.
}

#----------------------------------------------------------
# 12. Final candidate Cox models
#----------------------------------------------------------

fit_final_candidate_models <- function(data) {       # Define function to fit manually curated final candidate Cox models.
  final_model_1 <- coxph(
    Surv(survtime, censdead) ~ size + prog_recp + rectime.med_f +
      menopause_f + hormone_f + size:hormone_f + prog_recp:rectime.med_f,
    data = data
  )                                                 # Candidate model including median-based rectime, menopause, hormone, and key interactions.
  
  final_model_2 <- coxph(
    Surv(survtime, censdead) ~ size + prog_recp + rectime +
      prog_recp:rectime,
    data = data
  )                                                 # Simpler candidate model with continuous rectime and interaction with progesterone receptors.
  
  final_model_3 <- coxph(
    Surv(survtime, censdead) ~ size + prog_recp + rectime +
      menopause_f + hormone_f + size:hormone_f + prog_recp:rectime,
    data = data
  )                                                 # Candidate model combining continuous rectime, menopause, hormone, and key interactions.
  
  final_ph <- map(
    list(final_model_1, final_model_2, final_model_3),
    cox.zph
  )                                                 # Compute PH assumption tests for each candidate model; ensures validity of Cox framework.
  
  list(
    models = list(
      final_model_1 = final_model_1,                # Return first candidate model.
      final_model_2 = final_model_2,                # Return second candidate model.
      final_model_3 = final_model_3                 # Return third candidate model.
    ),
    ph     = final_ph                               # Return PH diagnostics for all candidate models.
  )
}

#----------------------------------------------------------
# 13. Final Cox model + PH + DFbeta
#----------------------------------------------------------

fit_final_cox_model <- function(data) {              # Define function to fit chosen final Cox model and compute PH and DFbeta diagnostics.
  fit_final_cox <- coxph(
    Surv(survtime, censdead) ~ size + prog_recp + rectime +
      hormone_f + age.med_f + prog_recp:rectime + size:hormone_f,
    data = data
  )                                                 # Final Cox model including size, progesterone receptors, rectime, hormone, age group, and key interactions.
  
  ph_final <- cox.zph(fit_final_cox)                # Compute PH assumption test for final model; checks time-constancy of covariate effects.
  
  dfbeta_final <- residuals(fit_final_cox, type = "dfbeta") %>%
    as_tibble() %>%
    mutate(obs = row_number())                      # Compute DFbeta residuals and add observation index; used to detect influential cases.
  
  list(
    model  = fit_final_cox,                         # Return final Cox model object.
    ph     = ph_final,                              # Return PH diagnostics for final model.
    dfbeta = dfbeta_final                           # Return DFbeta residuals for influence analysis.
  )
}

#----------------------------------------------------------
# 14. Weibull model corresponding to final Cox model
#----------------------------------------------------------

fit_final_weibull_model <- function(data) {          # Define function to fit Weibull parametric model with same covariates as final Cox.
  fit_final_weib <- survreg(
    Surv(survtime, censdead) ~ size + prog_recp + rectime +
      hormone_f + age.med_f + prog_recp:rectime + size:hormone_f,
    dist = "weibull",
    data = data
  )                                                 # Fit Weibull survival regression; provides parametric counterpart to Cox model.
  
  fit_final_weib                                   # Return fitted Weibull model object.
}

compare_weibull_cox <- function(fit_final_weib, fit_final_cox) { # Define function to compare Weibull-derived PH coefficients to Cox coefficients.
  weib_coef_all <- fit_final_weib$coefficients[-1]  # Extract Weibull coefficients excluding intercept; focus on covariate effects.
  weib_scale    <- fit_final_weib$scale             # Extract Weibull scale parameter; used to convert to PH form.
  weib_coef_ph  <- -weib_coef_all / weib_scale      # Convert Weibull coefficients to PH-equivalent form; aligns with Cox interpretation.
  
  cox_coef      <- coef(fit_final_cox)              # Extract Cox model coefficients; baseline for comparison.
  
  tibble(
    term       = names(weib_coef_ph),               # Store term names for which comparison is made.
    weibull_ph = as.numeric(weib_coef_ph),          # Store PH-equivalent Weibull coefficients.
    coxph_coef = as.numeric(cox_coef[term])         # Store corresponding Cox coefficients for same terms.
  )
}

#----------------------------------------------------------
# 15. Deviance residuals for Weibull model
#----------------------------------------------------------

compute_weibull_deviance_resid <- function(data, fit_final_weib) { # Define function to compute and visualize deviance residuals for Weibull model.
  dev_resid <- residuals(fit_final_weib, type = "deviance")        # Compute deviance residuals; measure goodness-of-fit at observation level.
  
  gbcs_dev <- data %>%
    mutate(dev_resid = dev_resid)                  # Add deviance residuals to dataset; enables plotting against covariates.
  
  dev_plot_hormone <- gbcs_dev %>%
    ggplot(aes(x = hormone_f, y = dev_resid)) +
    geom_boxplot(fill = "lightgray") +             # Create boxplot of deviance residuals by hormone therapy status.
    geom_jitter(width = 0.2, alpha = 0.4) +        # Add jittered points to show individual residuals and potential outliers.
    labs(
      title = "Deviance Residuals vs Hormone Therapy", # Title indicates relationship between residuals and hormone therapy.
      x = "Hormone Therapy",                       # X-axis label clarifies treatment groups.
      y = "Deviance Residuals"                     # Y-axis label indicates residual type.
    ) +
    theme_minimal()                                # Minimal theme for clean visualization.
  
  dev_plot_age_med <- gbcs_dev %>%
    ggplot(aes(x = age.med_f, y = dev_resid)) +
    geom_boxplot(fill = "lightgray") +             # Create boxplot of deviance residuals by median age group.
    geom_jitter(width = 0.2, alpha = 0.4) +        # Add jittered points to show individual residuals and potential outliers.
    labs(
      title = "Deviance Residuals vs Median Age at Diagnosis", # Title indicates relationship between residuals and age group.
      x = "Median Age Group",                      # X-axis label clarifies age categories.
      y = "Deviance Residuals"                     # Y-axis label indicates residual type.
    ) +
    theme_minimal()                                # Minimal theme for clean visualization.
  
  list(
    data              = gbcs_dev,                  # Return dataset with deviance residuals attached.
    dev_plot_hormone  = dev_plot_hormone,          # Return deviance residual plot vs hormone therapy.
    dev_plot_age_med  = dev_plot_age_med           # Return deviance residual plot vs median age group.
  )
}

#----------------------------------------------------------
# 16. Case-deletion (DFbeta) for Weibull model
#----------------------------------------------------------

weibull_dfbeta_diagnostics <- function(fit_final_weib) { # Define function to compute DFbeta diagnostics for Weibull model.
  dfbeta_weib <- residuals(fit_final_weib, type = "dfbeta") %>%
    as_tibble() %>%
    mutate(obs = row_number())                      # Compute DFbeta residuals and add observation index; used for influence analysis.
  
  dfbeta_weib_size_plot <- dfbeta_weib %>%
    ggplot(aes(x = obs, y = .data[[2]])) +
    geom_segment(aes(xend = obs, yend = 0)) +
    labs(
      title = "Case Deletion Plot for Tumor Size (Weibull)", # Title indicates influence of each case on size coefficient.
      x = "Observation Index",                      # X-axis label clarifies that each point corresponds to an observation.
      y = "Change in Size Coefficient"              # Y-axis label indicates DFbeta magnitude for size coefficient.
    ) +
    theme_minimal()                                # Minimal theme for clean diagnostic visualization.
  
  influential_size <- dfbeta_weib %>%
    transmute(
      obs,
      dfb_size = .data[[2]]
    ) %>%
    filter(dfb_size > 0.0005 | dfb_size < -0.0005) # Flag observations with DFbeta magnitude above threshold; potential influential cases.
  
  list(
    dfbeta_tbl   = dfbeta_weib,                    # Return full DFbeta table for Weibull model.
    size_plot    = dfbeta_weib_size_plot,          # Return case-deletion plot for size coefficient.
    influential  = influential_size                # Return subset of potentially influential observations.
  )
}

#----------------------------------------------------------
# 17. Residual plots for slides (martingale)
#----------------------------------------------------------

slide_residual_plots <- function(resid_data) {       # Define function to generate slide-ready martingale residual plots.
  slide_resid_rectime <- resid_data %>%
    ggplot(aes(x = rectime, y = martingale)) +
    geom_point(alpha = 0.6) +                       # Plot residuals vs rectime with semi-transparent points.
    geom_smooth(method = "loess", color = "red", se = FALSE, size = 1.2) + # Add LOESS smooth with thicker line for presentation.
    labs(
      title = "Martingale Residuals vs Time Until Recurrence", # Title suitable for slides; highlights key predictor.
      x = "Time Until Recurrence",                 # X-axis label clarifies predictor.
      y = "Martingale Residuals"                   # Y-axis label indicates residual type.
    ) +
    theme_minimal()                                # Minimal theme for clean slide visuals.
  
  slide_resid_prog <- resid_data %>%
    ggplot(aes(x = prog_recp, y = martingale)) +
    geom_point(alpha = 0.6) +                       # Plot residuals vs progesterone receptors with semi-transparent points.
    geom_smooth(method = "loess", color = "red", se = FALSE, size = 1.2) + # Add LOESS smooth for visualizing functional form.
    labs(
      title = "Martingale Residuals vs Progesterone Receptors", # Title suitable for slides; highlights biomarker predictor.
      x = "Number of Progesterone Receptors",      # X-axis label clarifies biomarker.
      y = "Martingale Residuals"                   # Y-axis label indicates residual type.
    ) +
    theme_minimal()                                # Minimal theme for clean slide visuals.
  
  list(
    rectime = slide_resid_rectime,                 # Return slide-ready plot for rectime.
    prog    = slide_resid_prog                     # Return slide-ready plot for progesterone receptors.
  )
}

#----------------------------------------------------------
# 18. Treatment comparison: Cox vs Weibull (hormone only)
#----------------------------------------------------------

treatment_comparison_hormone <- function(data) {     # Define function to compare Cox and Weibull survival curves for hormone therapy.
  cox_hormone <- coxph(Surv(survtime, censdead) ~ hormone_f, data = data) # Fit Cox model with hormone therapy as sole predictor.
  
  weib_hormone <- survreg(
    Surv(survtime, censdead) ~ hormone_f,
    dist = "weibull",
    data = data
  )                                                 # Fit Weibull parametric model with hormone therapy as sole predictor.
  
  mu0_hat     <- weib_hormone$coefficients[1]       # Extract Weibull intercept; corresponds to baseline log-scale parameter.
  sigma_hat   <- weib_hormone$scale                 # Extract Weibull scale parameter; controls spread of survival times.
  alpha_hat   <- 1 / sigma_hat                      # Compute Weibull shape parameter; inverse of scale in survreg parameterization.
  lambda0_hat <- exp(-mu0_hat)                      # Compute baseline hazard scale parameter; derived from intercept.
  
  tt_vec    <- 0:2601                               # Define time grid for survival curve evaluation; covers observed range.
  surv0_vec <- 1 - pweibull(tt_vec, shape = alpha_hat, scale = 1 / lambda0_hat) # Compute baseline survival curve under Weibull model.
  
  gamma_hat <- weib_hormone$coefficients[2]         # Extract hormone effect coefficient from Weibull model.
  surv1_vec <- surv0_vec^(exp(-gamma_hat / sigma_hat)) # Compute survival curve for hormone group using proportional hazards relationship.
  
  cox_surv_est <- survfit(
    cox_hormone,
    newdata = data.frame(hormone_f = c("No", "Yes"))
  )                                                 # Compute Cox survival estimates for hormone = No and Yes; used for comparison.
  
  list(
    cox_model     = cox_hormone,                    # Return Cox model object for hormone-only analysis.
    weib_model    = weib_hormone,                   # Return Weibull model object for hormone-only analysis.
    tt_vec        = tt_vec,                         # Return time grid used for survival curve evaluation.
    surv0_vec     = surv0_vec,                      # Return baseline Weibull survival curve.
    surv1_vec     = surv1_vec,                      # Return Weibull survival curve for hormone therapy group.
    cox_surv_est  = cox_surv_est                    # Return Cox survival estimates for both hormone groups.
  )
}

#----------------------------------------------------------
# 19. Orchestrator: full pipeline
#----------------------------------------------------------

run_gbcs_survival_pipeline <- function(path) {       # Define main orchestrator function to run full survival analysis pipeline on GBCS data.
  gbcs2 <- load_and_clean_gbcs(path)                # Load and preprocess dataset from given path; returns cleaned analysis-ready data.
  
  gbcs_resid <- compute_null_martingale(gbcs2)      # Fit null Cox model and compute martingale residuals; used for diagnostic plots.
  
  diag_plots      <- diagnostic_plots(gbcs_resid)   # Generate residual diagnostic plots for numeric and categorical predictors.
  multi_models    <- fit_multivariable_models(gbcs2) # Fit multivariable Cox models with different exclusion sets.
  reduced_models  <- fit_reduced_models(gbcs2)      # Fit reduced Cox models based on key covariates.
  
  change_results <- compute_change_in_estimate(
    full    = multi_models$models$all_vars_1,
    reduced = reduced_models$models$reduced_1
  )                                                 # Compute percent change in coefficients between full and reduced models; assess confounding.
  
  reexam_results  <- reexamine_excluded(gbcs2)      # Reintroduce excluded variables one at a time and refit models; evaluate their impact.
  spline_results  <- fit_spline_models(gbcs2)       # Fit spline-based models for key predictors; assess nonlinear effects.
  log_plots       <- log_diagnostic_plots(gbcs_resid) # Generate log-scale residual plots for numeric predictors.
  
  interaction_res <- interaction_screening(gbcs2)   # Screen candidate interaction terms using ANOVA comparisons.
  full_inter      <- fit_full_interaction_models(gbcs2) # Fit full interaction models and perform stepwise AIC/BIC selection.
  
  models_to_check <- list(
    final_interaction_1 = full_inter$step_aic,
    final_interaction_2 = full_inter$step_bic,
    base_no_interaction = interaction_res$base_model
  )                                                 # Define list of models for PH assumption checks; includes stepwise-selected and base models.
  
  ph_tests <- ph_assumption_checks(models_to_check) # Compute PH diagnostics for selected models; ensures validity of Cox assumptions.
  
  final_candidates <- fit_final_candidate_models(gbcs2) # Fit manually curated final candidate Cox models.
  final_cox        <- fit_final_cox_model(gbcs2)   # Fit chosen final Cox model and compute PH and DFbeta diagnostics.
  
  final_weib       <- fit_final_weibull_model(gbcs2) # Fit Weibull parametric model corresponding to final Cox model.
  coef_compare     <- compare_weibull_cox(final_weib, final_cox$model) # Compare Weibull-derived PH coefficients to Cox coefficients.
  
  dev_resid_res    <- compute_weibull_deviance_resid(gbcs2, final_weib) # Compute and plot deviance residuals for Weibull model.
  weib_dfbeta_res  <- weibull_dfbeta_diagnostics(final_weib)            # Compute DFbeta diagnostics for Weibull model.
  
  slide_plots      <- slide_residual_plots(gbcs_resid) # Generate slide-ready martingale residual plots for key predictors.
  treat_comp       <- treatment_comparison_hormone(gbcs2) # Compute Cox and Weibull survival curves for hormone therapy comparison.
  
  list(
    data              = gbcs2,                     # Return cleaned dataset used in analysis.
    resid_data        = gbcs_resid,                # Return dataset with martingale residuals attached.
    diagnostics       = diag_plots,                # Return residual diagnostic plots.
    multivariable     = multi_models,              # Return multivariable model fits and summaries.
    reduced           = reduced_models,            # Return reduced model fits and summaries.
    change_estimate   = change_results,            # Return percent change-in-estimate results.
    reexam_excluded   = reexam_results,            # Return re-examination of excluded variables.
    spline            = spline_results,            # Return spline-based model fits.
    log_plots         = log_plots,                 # Return log-scale residual plots.
    interactions      = interaction_res,           # Return interaction screening results.
    full_interactions = full_inter,                # Return full interaction model fits and stepwise results.
    ph_tests          = ph_tests,                  # Return PH diagnostics for selected models.
    final_candidates  = final_candidates,          # Return final candidate Cox models and PH checks.
    final_cox         = final_cox,                 # Return final Cox model, PH diagnostics, and DFbeta.
    final_weib        = final_weib,                # Return final Weibull model.
    coef_compare      = coef_compare,              # Return comparison of Weibull and Cox coefficients.
    deviance_resid    = dev_resid_res,             # Return deviance residuals and plots for Weibull model.
    weib_dfbeta       = weib_dfbeta_res,           # Return DFbeta diagnostics for Weibull model.
    slide_plots       = slide_plots,               # Return slide-ready residual plots.
    treatment_comp    = treat_comp                 # Return Cox vs Weibull hormone treatment comparison results.
  )
}

# Example call (uncomment and adjust path as needed):
# results <- run_gbcs_survival_pipeline(
#   "~/BIST0627 Applied Survival Data Analysis/BIST0627 Project/project_data_files/gbcs.csv"
# )
