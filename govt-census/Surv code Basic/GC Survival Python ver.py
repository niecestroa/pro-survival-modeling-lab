# Author: Aaron Niecestro
# Created on June 15, 2023
# Last Editted on February 14, 2026

# ======================================================
# Imports
# ======================================================
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")

# ======================================================
# Import and Prepare Data
# ======================================================

# Update path as needed
pbc = pd.read_csv(
    r"C:\Users\aniec\UTH Survival Analysis\PH1831_Project\PH1831 Project Data\PBC276.csv"
)

# Clean column names to snake_case (simple version)
pbc.columns = (
    pbc.columns.str.strip()
    .str.lower()
    .str.replace(" ", "_")
    .str.replace(".", "_")
)

# Convert to categorical
for col in ["stage", "drug", "sex"]:
    pbc[col] = pbc[col].astype("category")

# ======================================================
# Null Cox Model + Martingale Residuals
# ======================================================

cph_null = CoxPHFitter()
cph_null.fit(pbc, duration_col="futime", event_col="status", formula="1")

pbc["martingale"] = cph_null.compute_residuals(pbc, "martingale")

# ======================================================
# Martingale Residual Plots (Linear + Log Scale)
# ======================================================

def plot_resid(varname):
    x = pbc[varname]
    y = pbc["martingale"]
    plt.figure()
    sns.scatterplot(x=x, y=y, alpha=0.6)
    sns.regplot(x=x, y=y, scatter=False, lowess=True, color="blue")
    plt.xlabel(varname)
    plt.ylabel("Martingale Residuals")
    plt.title(f"Martingale Residuals vs {varname}")
    plt.show()

def plot_resid_log(varname):
    x = pbc[varname]
    y = pbc["martingale"]
    plt.figure()
    plt.xscale("log")
    sns.scatterplot(x=x, y=y, alpha=0.6)
    sns.regplot(x=x, y=y, scatter=False, lowess=True, color="red")
    plt.xlabel(f"log({varname})")
    plt.ylabel("Martingale Residuals")
    plt.title(f"Martingale Residuals vs log({varname})")
    plt.show()

vars_to_check = ["age", "bili", "copper", "albumin", "protime"]

for v in vars_to_check:
    plot_resid(v)
for v in vars_to_check:
    plot_resid_log(v)

# ======================================================
# Spline Fits to Assess Nonlinearity (approximate)
# ======================================================
# lifelines does not have pspline directly; you’d typically
# use patsy or manually create spline bases. Placeholder:
# this section is conceptually acknowledged but not fully
# replicated here.

# ======================================================
# Build Full Cox Model + Stepwise Selection (AIC-style)
# ======================================================

# Full model
cph_full = CoxPHFitter()
cph_full.fit(
    pbc,
    duration_col="futime",
    event_col="status",
    formula="bili + copper + albumin + protime + age + stage + drug",
)

# Simple manual stepwise-style helper (backward by AIC)
def backward_stepwise_aic(df, duration, event, start_formula):
    current_formula = start_formula
    best_aic = np.inf
    improved = True

    while improved:
        improved = False
        terms = [t.strip() for t in current_formula.split("+")]
        candidates = []

        for term in terms:
            reduced_terms = [t for t in terms if t != term]
            if not reduced_terms:
                continue
            formula = " + ".join(reduced_terms)
            model = CoxPHFitter()
            model.fit(df, duration_col=duration, event_col=event, formula=formula)
            candidates.append((formula, model.AIC_))

        if not candidates:
            break

        best_candidate = min(candidates, key=lambda x: x[1])
        if best_candidate[1] + 1e-6 < best_aic:
            best_aic = best_candidate[1]
            current_formula = best_candidate[0]
            improved = True

    final_model = CoxPHFitter()
    final_model.fit(df, duration_col=duration, event_col=event, formula=current_formula)
    return final_model, current_formula

full_formula = "bili + copper + albumin + protime + age + stage + drug"
step_model_full, step_formula_full = backward_stepwise_aic(
    pbc, "futime", "status", full_formula
)

# ======================================================
# Add Log Transformations + Refit
# ======================================================

pbc["log_bili"] = np.log(pbc["bili"])
pbc["log_copper"] = np.log(pbc["copper"])

trans_formula = "albumin + protime + age + log_bili + log_copper + stage + drug"
step_model_trans, step_formula_trans = backward_stepwise_aic(
    pbc, "futime", "status", trans_formula
)

# ======================================================
# Final Models
# ======================================================

final_formula_1 = "albumin + protime + age + log_bili + log_copper + drug"
finalmodel_pbc1 = CoxPHFitter()
finalmodel_pbc1.fit(pbc, duration_col="futime", event_col="status", formula=final_formula_1)
print(finalmodel_pbc1.summary)

final_formula_2 = "albumin + protime + age + log_bili + log_copper"
finalmodel_pbc2 = CoxPHFitter()
finalmodel_pbc2.fit(pbc, duration_col="futime", event_col="status", formula=final_formula_2)
print(finalmodel_pbc2.summary)

# ======================================================
# Proportional Hazards Assumption
# ======================================================

ph1 = proportional_hazard_test(finalmodel_pbc1, pbc, time_transform="rank")
ph2 = proportional_hazard_test(finalmodel_pbc2, pbc, time_transform="rank")

print(ph1.summary)
print(ph2.summary)

# ======================================================
# DFbeta Residuals
# ======================================================

dfb = finalmodel_pbc1.compute_residuals(pbc, "dfbeta")
dfb = dfb.reset_index(drop=True)
dfb["obs"] = dfb.index + 1

coef_names = list(finalmodel_pbc1.params_.index)

def plot_dfbeta(df, var_label):
    plt.figure()
    sns.scatterplot(x=df["obs"], y=df[var_label], s=10)
    plt.axhline(0, color="red", linestyle="--")
    plt.title(f"DFbeta for {var_label}")
    plt.xlabel("Observation Index")
    plt.ylabel(f"DFbeta: {var_label}")
    plt.show()

for name in coef_names:
    plot_dfbeta(dfb, name)

# ======================================================
# Deviance Residuals
# ======================================================

dev = finalmodel_pbc1.compute_residuals(pbc, "deviance")
dev_resid = pd.DataFrame({
    "dev": dev,
    "age": pbc["age"],
    "protime": pbc["protime"],
    "albumin": pbc["albumin"],
    "log_bili": pbc["log_bili"],
    "log_copper": pbc["log_copper"],
})

def plot_dev(xvar):
    plt.figure()
    sns.scatterplot(x=dev_resid[xvar], y=dev_resid["dev"], alpha=0.6)
    sns.regplot(x=dev_resid[xvar], y=dev_resid["dev"], scatter=False, lowess=True, color="blue")
    plt.title(f"Deviance Residuals vs {xvar}")
    plt.xlabel(xvar)
    plt.ylabel("Deviance Residuals")
    plt.show()

for v in ["age", "protime", "albumin", "log_bili", "log_copper"]:
    plot_dev(v)

# ======================================================
# Case-Deletion Plots (using DFbeta as proxy)
# ======================================================

def plot_case_deletion(df, var_label):
    plt.figure()
    for i in range(len(df)):
        plt.plot([df["obs"].iloc[i], df["obs"].iloc[i]],
                 [0, df[var_label].iloc[i]],
                 color="black", alpha=0.7)
    plt.axhline(0, color="red")
    plt.title(f"Case Deletion Plot for {var_label}")
    plt.xlabel("Observation Index")
    plt.ylabel(f"Change in {var_label} Coefficient")
    plt.show()

for name in coef_names:
    plot_case_deletion(dfb, name)

# ======================================================
# Survival Curves for Drug × Sex
# ======================================================

km = KaplanMeierFitter()

# Create combined group variable like drug + sex
pbc["drug_sex"] = pbc["drug"].astype(str) + "_" + pbc["sex"].astype(str)

plt.figure()
for group, dfg in pbc.groupby("drug_sex"):
    km.fit(durations=dfg["futime"], event_observed=dfg["status"], label=group)
    km.plot(ci_show=False)
plt.title("Survival Curves by Treatment and Sex")
plt.xlabel("Time")
plt.ylabel("Survival Probability")
plt.show()

# ======================================================
# PH Diagnostics Plots (approximate ggcoxzph)
# ======================================================
# lifelines doesn't have ggcoxzph, but you can inspect scaled Schoenfeld
# residuals via proportional_hazard_test and custom plotting if desired.
