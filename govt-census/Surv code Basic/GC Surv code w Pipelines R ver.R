#======================================================================
# Author: Aaron Niecestro
# Created:      June 15, 2023
# Last Edit:    September 14, 2026
#
# Project: US Government Survival Analysis – PBC276
#======================================================================

#----------------------------------------------------------------------
# 0. Setup
#----------------------------------------------------------------------

set.seed(1831)                     # ensures reproducibility for smoothers and resampling

library(tidyverse)                 # data manipulation, pipes, ggplot2
library(survival)                  # Cox proportional hazards models
library(survminer)                 # ggplot-based survival visualizations
library(asaur)                    # additional survival analysis utilities
library(janitor)                   # clean column names
library(broom)                     # tidy model outputs (tidy/glance/augment)

#----------------------------------------------------------------------
# 1. Import and Prepare Data
#----------------------------------------------------------------------

PBC <- read_csv("~/UTH Survival Analysis/PH1831_Project/PH1831 Project Data/PBC276.csv") %>%
  clean_names() %>%                # convert column names to snake_case
  mutate(
    stage      = factor(stage),    # convert stage to factor
    drug       = factor(drug),     # convert drug to factor
    sex        = factor(sex),      # convert sex to factor
    log_bili   = log(bili),        # log-transform bilirubin
    log_copper = log(copper)       # log-transform copper
  )

glimpse(PBC)                       # quick structural check of dataset

#----------------------------------------------------------------------
# 2. Null Cox Model + Martingale Residuals
#   Purpose: assess functional form of continuous covariates
#----------------------------------------------------------------------

null_model <- coxph(Surv(futime, status) ~ 1, data = PBC)   # baseline hazard only
martingale_resid <- residuals(null_model, type = "martingale")  # extract martingale residuals

# Build long-format tibble for faceted plotting
mart_df <- PBC %>%
  select(age, bili, copper, albumin, protime) %>%            # continuous covariates
  mutate(rr = martingale_resid) %>%                          # add residuals
  pivot_longer(-rr, names_to = "variable", values_to = "x")  # reshape for ggplot

# Plot martingale residuals (linear scale)
ggplot(mart_df, aes(x = x, y = rr)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Martingale Residuals vs Continuous Covariates",
    x = "Covariate Value",
    y = "Martingale Residuals"
  ) +
  theme_bw()

# Build log-transformed version for skewed variables
mart_df_log <- PBC %>%
  transmute(
    age,
    log_bili,
    log_copper,
    albumin,
    protime,
    rr = martingale_resid
  ) %>%
  pivot_longer(-rr, names_to = "variable", values_to = "x")

# Plot martingale residuals (log scale)
ggplot(mart_df_log, aes(x = x, y = rr)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Martingale Residuals vs Log/Transformed Covariates",
    x = "Transformed Covariate",
    y = "Martingale Residuals"
  ) +
  theme_bw()

#----------------------------------------------------------------------
# 3. Spline Fits to Assess Nonlinearity
#----------------------------------------------------------------------

vars_cont <- c("age", "bili", "copper", "albumin", "protime")  # continuous variables

# Function to fit spline and plot effect
fit_spline <- function(varname) {
  formula <- as.formula(paste0("Surv(futime, status) ~ pspline(", varname, ")"))  # build formula
  m <- coxph(formula, data = PBC)                                                 # fit spline model
  termplot(m, term = 1, se = TRUE, col.term = 1, col.se = "blue",                # plot spline effect
           main = paste("Spline Effect of", varname))
}

walk(vars_cont, fit_spline)                  # apply spline diagnostic to all variables

#----------------------------------------------------------------------
# 4. Full Cox Model + Stepwise Selection
#   Note: stepwise is included for teaching; penalized Cox is preferred
#----------------------------------------------------------------------

model_full <- coxph(
  Surv(futime, status) ~ bili + copper + albumin + protime +
    age + stage + drug,
  data = PBC
)

model_full_step <- step(model_full, direction = "both", trace = FALSE)  # AIC-based selection
summary(model_full_step)                                                # model summary

#----------------------------------------------------------------------
# 5. Log-Transformed Model + Stepwise Selection
#----------------------------------------------------------------------

model_transformed <- coxph(
  Surv(futime, status) ~ albumin + protime + age +
    log_bili + log_copper + stage + drug,
  data = PBC
)

model_transformed_step <- step(model_transformed, direction = "both", trace = FALSE)
summary(model_transformed_step)

#----------------------------------------------------------------------
# 6. Final Models (Nested)
#   finalmodel_pbc1: includes drug
#   finalmodel_pbc2: excludes drug
#----------------------------------------------------------------------

finalmodel_pbc1 <- coxph(
  Surv(futime, status) ~ albumin + protime + age +
    log_bili + log_copper + drug,
  data = PBC,
  robust = TRUE                     # robust standard errors
)

finalmodel_pbc2 <- coxph(
  Surv(futime, status) ~ albumin + protime + age +
    log_bili + log_copper,
  data = PBC,
  robust = TRUE
)

summary(finalmodel_pbc1)            # detailed model output
summary(finalmodel_pbc2)

tidy(finalmodel_pbc1)               # tidy coefficient table
glance(finalmodel_pbc1)             # model-level statistics

#----------------------------------------------------------------------
# 7. Proportional Hazards Assumption
#----------------------------------------------------------------------

ph1 <- cox.zph(finalmodel_pbc1)     # Schoenfeld residual test
ph2 <- cox.zph(finalmodel_pbc2)

print(ph1)                          # print PH test results
print(ph2)

ggcoxzph(ph1)                       # ggplot version of PH diagnostics

#----------------------------------------------------------------------
# 8. DFbeta Residuals (Influence Diagnostics)
#----------------------------------------------------------------------

dfb <- residuals(finalmodel_pbc1, type = "dfbeta") %>%  # extract dfbeta residuals
  as_tibble() %>%                                       # convert to tibble
  mutate(obs = row_number())                            # add observation index

dfb_long <- dfb %>%
  pivot_longer(-obs, names_to = "coef", values_to = "dfb")  # reshape for plotting

# Plot DFbeta residuals
ggplot(dfb_long, aes(x = obs, y = dfb)) +
  geom_point(size = 0.7) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~ coef, scales = "free_y") +
  labs(
    title = "DFbeta Residuals by Coefficient",
    x = "Observation Index",
    y = "DFbeta"
  ) +
  theme_bw()

#----------------------------------------------------------------------
# 9. Deviance Residuals vs Covariates
#----------------------------------------------------------------------

dev_resid <- tibble(
  dev        = residuals(finalmodel_pbc1, type = "deviance"),  # deviance residuals
  age        = PBC$age,
  protime    = PBC$protime,
  albumin    = PBC$albumin,
  log_bili   = PBC$log_bili,
  log_copper = PBC$log_copper
)

dev_resid_long <- dev_resid %>%
  pivot_longer(-dev, names_to = "variable", values_to = "x")  # reshape for plotting

# Plot deviance residuals
ggplot(dev_resid_long, aes(x = x, y = dev)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Deviance Residuals vs Covariates",
    x = "Covariate",
    y = "Deviance Residuals"
  ) +
  theme_bw()

#----------------------------------------------------------------------
# 10. Case-Deletion Visualization (via DFbeta)
#----------------------------------------------------------------------

ggplot(dfb_long, aes(x = obs, y = dfb)) +
  geom_segment(aes(xend = obs, yend = 0), alpha = 0.7) +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(~ coef, scales = "free_y") +
  labs(
    title = "Case Deletion Visualization via DFbeta",
    x = "Observation Index",
    y = "Change in Coefficient"
  ) +
  theme_bw()

#----------------------------------------------------------------------
# 11. Survival Curves for Drug × Sex
#----------------------------------------------------------------------

fit_drug_sex <- survfit(Surv(futime, status) ~ drug + sex, data = PBC)  # fit survival curves

ggsurvplot(
  fit_drug_sex,
  pval          = TRUE,             # show log-rank p-value
  break.time.by = 400,              # x-axis breaks
  risk.table    = FALSE,            # hide risk table
  palette       = "Dark2",          # color palette
  title         = "Survival Curves by Treatment and Sex"
)

#----------------------------------------------------------------------
# 12. Concordance (Model Discrimination)
#----------------------------------------------------------------------

concordance_pbc1 <- survConcordance(
  Surv(futime, status) ~ predict(finalmodel_pbc1),  # predicted risk scores
  data = PBC
)

print(concordance_pbc1)             # print concordance index

#----------------------------------------------------------------------
# 13. Save Final Models (Optional)
#----------------------------------------------------------------------

saveRDS(finalmodel_pbc1, file = "finalmodel_pbc1.rds")   # save model with drug
saveRDS(finalmodel_pbc2, file = "finalmodel_pbc2.rds")   # save model without drug
saveRDS(ph1, file = "ph1_diagnostics.rds")               # save PH diagnostics
