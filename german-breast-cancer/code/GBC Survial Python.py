# Author: Aaron Niecestro
# Created on January 2, 2024
# Last Editted on February 15, 2026

# ---------------------------------------------------------
# 0. Imports
# ---------------------------------------------------------
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, WeibullAFTFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------
# 1. Load & Clean Data  (R: mutate() / SAS: DATA step)
# ---------------------------------------------------------
df = pd.read_csv(
    r"C:\Users\aniec\BIST0627 Applied Survival Data Analysis\BIST0627 Project\project_data_files\gbcs.csv"
)

# Convert numeric columns
num_cols = ["age","size","nodes","prog_recp","estrg_recp","rectime","survtime","censdead"]
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")

# Factor variables
df["menopause_f"] = np.where(df["menopause"]=="1","Yes","No")
df["hormone_f"]   = np.where(df["hormone"]=="1","Yes","No")
df["censrec_f"]   = np.where(df["censrec"]=="1","Recurrence","Censored")
df["grade_f"]     = df["grade"].astype("category")

# Median splits
df["age_med_f"]     = np.where(df["age"]>53,"Above(>53)","Below(<=53)")
df["rectime_med_f"] = np.where(df["rectime"]>1084,"Above","Below")

# Drop unused columns
df = df.drop(columns=["id","diagdateb","recdate","deathdate","menopause","hormone","grade","censrec"])

# ---------------------------------------------------------
# 2. Null Cox Model  (R: coxph(~1) / SAS: PROC PHREG)
# ---------------------------------------------------------
cph_null = CoxPHFitter()
cph_null.fit(df, duration_col="survtime", event_col="censdead", formula="1")

# Martingale residuals
df["martingale"] = cph_null.compute_residuals(df, "martingale")

# ---------------------------------------------------------
# 3. Residual Diagnostics (R: ggplot + lowess / SAS: SGPLOT)
# ---------------------------------------------------------
def martingale_plot(var):
    sns.scatterplot(x=df[var], y=df["martingale"], alpha=0.5)
    sns.regplot(x=df[var], y=df["martingale"], scatter=False, lowess=True, color="red")
    plt.title(f"Martingale Residuals vs {var}")
    plt.show()

for v in ["age","size","nodes","prog_recp","estrg_recp","rectime","survtime"]:
    martingale_plot(v)

# ---------------------------------------------------------
# 4. Full Cox Model  (R: coxph(~.) / SAS: PHREG)
# ---------------------------------------------------------
cph_full = CoxPHFitter()
cph_full.fit(df, duration_col="survtime", event_col="censdead")

# ---------------------------------------------------------
# 5. Reduced Model (R: size + prog_recp + rectime)
# ---------------------------------------------------------
cph_reduced = CoxPHFitter()
cph_reduced.fit(df, duration_col="survtime", event_col="censdead",
                formula="size + prog_recp + rectime")

# ---------------------------------------------------------
# 6. Interaction Screening (R: size*hormone_f / SAS: size*hormone_f)
# ---------------------------------------------------------
interaction_terms = [
    "size:hormone_f",
    "size:menopause_f",
    "size:grade_f",
    "prog_recp:rectime",
    "size:prog_recp",
]

interaction_models = {}

for term in interaction_terms:
    formula = f"size + prog_recp + rectime + {term}"
    model = CoxPHFitter()
    model.fit(df, duration_col="survtime", event_col="censdead", formula=formula)
    interaction_models[term] = model

# ---------------------------------------------------------
# 7. PH Assumption Tests (R: cox.zph / SAS: ASSESS PH)
# ---------------------------------------------------------
results_ph = proportional_hazard_test(cph_reduced, df, time_transform="rank")
print(results_ph.summary)

# ---------------------------------------------------------
# 8. Final Model (matches your R/SAS final)
# ---------------------------------------------------------
final_formula = """
size + prog_recp + rectime + hormone_f + age_med_f +
prog_recp:rectime + size:hormone_f
"""

cph_final = CoxPHFitter()
cph_final.fit(df, duration_col="survtime", event_col="censdead",
              formula=final_formula)

# DFbeta residuals
dfbeta = cph_final.compute_residuals(df, "dfbeta")

# ---------------------------------------------------------
# 9. Weibull Model (R: survreg / SAS: LIFEREG)
# ---------------------------------------------------------
weib = WeibullAFTFitter()
weib.fit(df, duration_col="survtime", event_col="censdead",
         formula=final_formula)

# Compare Cox vs Weibull coefficients
cox_coef = cph_final.params_
weib_coef = weib.params_.iloc[1:]  # drop intercept

comparison = pd.DataFrame({
    "cox_coef": cox_coef,
    "weibull_ph_coef": -weib_coef / weib.lambda_
})

print(comparison)

# ---------------------------------------------------------
# 10. Cox vs Weibull Survival Curves (R: survfit / SAS: LIFETEST)
# ---------------------------------------------------------
km = KaplanMeierFitter()

for group in ["Yes","No"]:
    mask = df["hormone_f"] == group
    km.fit(df.loc[mask,"survtime"], df.loc[mask,"censdead"], label=f"Cox KM: {group}")
    km.plot()

# Weibull curves
timeline = np.arange(0, df["survtime"].max(), 10)
surv_weib = weib.predict_survival_function(df, times=timeline)
surv_weib.mean(axis=1).plot(label="Weibull Mean Curve")

plt.legend()
plt.title("Cox vs Weibull Survival Comparison")
plt.show()
