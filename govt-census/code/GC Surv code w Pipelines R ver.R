# Author: Aaron Niecestro
# Created on June 15, 2023
# Last Editted on February 14, 2026

#======================================================================
# Load Required Packages
#======================================================================

library(tidyverse)     # data manipulation, pipes, plotting
library(survival)      # Cox proportional hazards models
library(survminer)     # ggplot-based survival visualizations
library(asaur)         # additional survival analysis utilities
library(janitor)       # clean column names


#======================================================================
# Import and Prepare Data
#======================================================================

PBC276 <- read_csv("~/UTH Survival Analysis/PH1831_Project/PH1831 Project Data/PBC276.csv") %>%
  clean_names() %>%                     # convert column names to snake_case
  mutate(
    stage = factor(stage),              # convert stage to factor
    drug  = factor(drug),               # convert drug to factor
    sex   = factor(sex)                 # convert sex to factor
  )


#======================================================================
# Null Cox Model + Martingale Residuals
#======================================================================

null_model <- PBC276 %>%
  coxph(Surv(futime, status) ~ 1, data = .)   # baseline hazard only

rr <- residuals(null_model, type = "martingale")  # martingale residuals


#======================================================================
# Martingale Residual Plots (Linear + Log Scale)
#======================================================================

# Helper function: linear-scale residual plot
plot_resid <- function(varname) {
  x <- PBC276[[varname]]                       # extract variable
  plot(x, rr,                                  # scatterplot
       xlab = varname,
       ylab = "Martingale Residuals")
  lines(lowess(x, rr), col = "blue", lwd = 2)  # smooth trend line
}

# Helper function: log-scale residual plot
plot_resid_log <- function(varname) {
  x <- PBC276[[varname]]
  plot(x, rr,
       log = "x",                              # log-transform x-axis
       xlab = paste0("log(", varname, ")"),
       ylab = "Martingale Residuals")
  lines(lowess(x, rr), col = "red", lwd = 2)
}

# Variables to check
vars <- c("age", "bili", "copper", "albumin", "protime")

# Generate both sets of plots
walk(vars, plot_resid)
walk(vars, plot_resid_log)


#======================================================================
# Spline Fits to Assess Nonlinearity
#======================================================================

fit_spline <- function(varname) {
  formula <- as.formula(paste0("Surv(futime, status) ~ pspline(", varname, ")"))
  
  PBC276 %>%
    coxph(formula, data = .) %>%               # fit spline model
    { termplot(., term = 1, se = TRUE,         # plot spline effect
               col.term = 1, col.se = "blue") }
}

walk(vars, fit_spline)


#======================================================================
# Build Full Cox Model + Stepwise Selection
#======================================================================

model_full <- PBC276 %>%
  coxph(Surv(futime, status) ~ bili + copper + albumin + protime +
          age + stage + drug, data = .)

step(model_full, direction = "both")           # AIC-based selection


#======================================================================
# Add Log Transformations + Refit
#======================================================================

model_transformed <- PBC276 %>%
  mutate(
    log_bili   = log(bili),
    log_copper = log(copper)
  ) %>%
  coxph(Surv(futime, status) ~ albumin + protime + age +
          log_bili + log_copper + stage + drug, data = .)

step(model_transformed, direction = "both")


#======================================================================
# Final Models
#======================================================================

finalmodel_pbc1 <- PBC276 %>%
  mutate(
    log_bili   = log(bili),
    log_copper = log(copper)
  ) %>%
  coxph(Surv(futime, status) ~ albumin + protime + age +
          log_bili + log_copper + drug, data = .)

summary(finalmodel_pbc1)

finalmodel_pbc2 <- PBC276 %>%
  mutate(
    log_bili   = log(bili),
    log_copper = log(copper)
  ) %>%
  coxph(Surv(futime, status) ~ albumin + protime + age +
          log_bili + log_copper, data = .)

summary(finalmodel_pbc2)


#======================================================================
# Proportional Hazards Assumption
#======================================================================

ph1 <- cox.zph(finalmodel_pbc1)
ph2 <- cox.zph(finalmodel_pbc2)

ph1
ph2


#======================================================================
# DFbeta Residuals
#======================================================================

dfb <- residuals(finalmodel_pbc1, type = "dfbeta") %>%
  as_tibble() %>%
  mutate(obs = row_number())                    # add observation index

coef_names <- names(coef(finalmodel_pbc1))      # coefficient names


# Helper function: DFbeta plot
plot_dfbeta <- function(df, var_index, var_label) {
  ggplot(df, aes(x = obs, y = .data[[var_index]])) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = paste("DFbeta for", var_label),
      x = "Observation Index",
      y = paste("DFbeta:", var_label)
    ) +
    theme_bw()
}

# Generate DFbeta plots
dfbeta_plots <- map2(seq_along(coef_names), coef_names,
                     ~ plot_dfbeta(dfb, .x, .y))


#======================================================================
# Deviance Residuals
#======================================================================

dev_resid <- tibble(
  dev = residuals(finalmodel_pbc1, type = "deviance"),
  age = PBC276$age,
  protime = PBC276$protime,
  albumin = PBC276$albumin,
  log_bili = log(PBC276$bili),
  log_copper = log(PBC276$copper)
)

# Helper function: deviance residual plot
plot_dev <- function(xvar) {
  ggplot(dev_resid, aes(x = .data[[xvar]], y = dev)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(
      title = paste("Deviance Residuals vs", xvar),
      x = xvar,
      y = "Deviance Residuals"
    ) +
    theme_bw()
}

dev_plots <- map(c("age", "protime", "albumin", "log_bili", "log_copper"),
                 plot_dev)


#======================================================================
# Case-Deletion Plots
#======================================================================

plot_case_deletion <- function(df, var_index, var_label) {
  ggplot(df, aes(x = obs, y = .data[[var_index]])) +
    geom_segment(aes(xend = obs, yend = 0), alpha = 0.7) +
    geom_hline(yintercept = 0, color = "red") +
    labs(
      title = paste("Case Deletion Plot for", var_label),
      x = "Observation Index",
      y = paste("Change in", var_label, "Coefficient")
    ) +
    theme_bw()
}

case_deletion_plots <- map2(seq_along(coef_names), coef_names,
                            ~ plot_case_deletion(dfb, .x, .y))

#======================================================================
# Survival Curves for Drug × Sex
#======================================================================

fit2 <- survfit(Surv(futime, status) ~ drug + sex, data = PBC276)

ggsurvplot(
  fit2,
  pval = TRUE,
  break.time.by = 400,
  risk.table = FALSE,
  palette = "Dark2",
  title = "Survival Curves by Treatment and Sex"
)

#======================================================================
# PH Diagnostics (Plots)
#======================================================================

ggcoxzph(ph1)     # ggplot version of PH assumption plots
