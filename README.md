# Survival Analysis Modeling Lab  
Advanced, multi‑language survival‑analysis pipelines using ASA‑aligned methods, including Cox PH modeling, Kaplan–Meier estimation, Weibull regression, diagnostics, and reproducible time‑to‑event workflows.

---

# **Overview**

This repository contains **professional, production‑ready survival‑analysis pipelines** built in **R, SAS, and Python**, demonstrating modern statistical modeling, reproducible engineering, and cross‑language mastery. Each project began as academic coursework and was systematically refactored into a modular, automated, and clinically aligned workflow.

The repo currently includes:

- **German Breast Cancer Survival Analysis**  
- **Government Census Survival Analysis**

Both projects follow the same engineering pattern:

1. **Academic R version** — exploratory, concept‑focused  
2. **Professional R pipeline** — modular, automated, reproducible  
3. **SAS translation** — clinical‑grade, regulatory‑aligned  
4. **Python implementation** — modern, scriptable analytics  

This structure highlights the progression from learning statistical methods to building **industry‑grade survival‑analysis systems**.

---

# **Project Evolution — Academic → Professional → Cross‑Language Pipelines**

## **1. Academic Foundations — R (Exploratory Coursework)**  
All survival projects began as graduate‑level assignments written in base R:

- Sequential, line‑by‑line scripts  
- Manual Cox PH fitting and diagnostics  
- Spline checks, interactions, Weibull modeling  
- Presentation‑oriented figures  
- No automation or modularity  

These versions emphasize **learning the methods**, not engineering reproducible workflows.

---

## **2. Professional Refactor — R Pipeline (Modular, Automated, Reproducible)**  
Each project was rebuilt into a **clinical‑grade R pipeline**:

- Modular scripts (`01_load_data.R`, `02_clean_data.R`, etc.)  
- Automated execution via a pipeline runner  
- Tidyverse‑based data preparation  
- purrr‑driven model loops  
- ggplot‑based diagnostics  
- Structured model selection (AIC/BIC, PH tests, interactions)  
- Standardized folder structure and reproducible outputs  

This version represents the transition from academic exploration to **production‑ready workflow engineering**.

---

## **3. SAS Translation — Clinical‑Grade, Regulatory‑Aligned Workflow**  
The R pipeline was translated into **SAS**, mirroring workflows used in pharmaceutical and regulatory environments:

- Controlled data preparation (`PROC IMPORT`, `DATA` steps)  
- `PROC PHREG` for Cox modeling and influence diagnostics  
- `PROC LIFEREG` for Weibull regression  
- `PROC LIFETEST` for Kaplan–Meier estimation  
- SGPLOT‑based diagnostics  
- FitStatistics‑based AIC/BIC model comparison  
- Case‑deletion and DFbeta influence analysis  

This version demonstrates the ability to implement survival analysis in an **audit‑ready clinical setting**.

---

## **4. Python Implementation — Modern, Scriptable Analytics**  
The pipeline was then ported into **Python**, completing the cross‑language workflow:

- `pandas` for data ingestion and cleaning  
- `lifelines` for Cox PH, Weibull, and diagnostics  
- `matplotlib` / `seaborn` for visualization  
- Automated model loops and reproducible reporting  
- Cross‑validated survival‑curve comparisons  

This version highlights **modern analytics engineering** and cross‑platform reproducibility.

---

# **Cross‑Language Comparison**

| Component | Academic R | Professional R | Professional SAS | Professional Python |
|----------|-------------|----------------|------------------|---------------------|
| Data Import | `read.csv()` | `readr::read_csv()` | `PROC IMPORT` | `pd.read_csv()` |
| Data Cleaning | Inline mutate | Modular script | `DATA` step | `df.assign()` |
| Factor Handling | `factor()` | `forcats` | `FORMAT` / `CLASS` | `astype('category')` |
| KM Estimation | `survfit()` | Modular KM | `PROC LIFETEST` | `KaplanMeierFitter()` |
| Cox PH | `coxph()` | Pipeline module | `PROC PHREG` | `CoxPHFitter()` |
| Residuals | `residuals()` | ggplot diagnostics | `OUTPUT resmart=` | `compute_residuals()` |
| PH Assumption | `cox.zph()` | Automated checks | `ASSESS PH` | `proportional_hazard_test()` |
| DFbeta | `dfbeta` residuals | Automated loops | `OUTPUT dfbeta=` | `compute_residuals("dfbeta")` |
| Splines | `pspline()` | Modular script | `EFFECT spl=Spline()` | `patsy.bs()` |
| Interactions | Manual | purrr loops | `size*hormone_f` | Formula interactions |
| Weibull | `survreg()` | Pipeline module | `PROC LIFEREG` | `WeibullAFTFitter()` |
| Model Selection | Manual | AIC/BIC pipeline | FitStatistics | `.AIC_`, `.BIC_` |
| Visualization | Base R | ggplot2 | SGPLOT | matplotlib / seaborn |
| Automation | None | `purrr::map()` | SAS macros | Python functions |

---

# **ASA Tree‑Based Machine Learning Methods (Completed)**

This repository also includes a full implementation of the  
**ASA Traveling Course: Tree‑Based Machine Learning Methods**  
by **Hemant Ishwaran** and **Min Lu** (University of Miami).

The ASA guide states:

> “RF‑SRC supplies unified forests for regression, classification, survival, and competing risks; VarPro supplies observed‑data variable priority; SGT supplies multivariate geometric tree splits; and RHF extends forest modeling to time‑varying hazard estimation.”

All components listed in the ASA ecosystem have been **successfully integrated** into the master workflow.

> **Note:**  
> All ASA implementations in this repository are currently written **exclusively in R**.  
> Equivalent Python implementations have **not yet been developed**, as the ASA ecosystem (RF‑SRC, VarPro, SGT, RHF) is presently available only in R.

---

## **Random Survival Forest Extensions (RF‑SRC) — Completed**
- Log‑rank, Brier‑score gradient, and log‑rank‑score split rules  
- OOB CRPS and time‑dependent performance metrics  
- Survival, regression, classification, multivariate, and competing‑risk forests  
- Full compatibility with GBCS dataset via ASA‑ready preprocessing  

---

## **VarPro Variable Selection — Completed**
- Observed‑data rule‑release variable priority  
- Case‑specific iVarPro importance  
- SHAP‑style iVarPro summaries  
- Unsupervised UVarPro  
- Cross‑validated VarPro (cv.varpro)  
- Liberal, conservative, and entropy‑based importance measures  

---

## **Super Greedy Trees (SGT) — Completed**
- Geometric split dictionaries  
- Hyperplane and higher‑order split families  
- Local polynomial and lasso‑based contributions  
- hcut tuning and basis‑function filtering  
- SGT model explainers (beta, partial effects)  

---

## **Random Hazard Forests (RHF) — Completed**
- Time‑varying hazard estimation  
- Counting‑process survival modeling  
- Time‑localized variable importance  
- Smoothed hazard curves  
- Time‑dependent AUC‑t evaluation  

---

## **Status**
All ASA Traveling Course components — RF‑SRC, VarPro, SGT, and RHF — are now fully implemented, documented, and integrated into the master workflow. The repository demonstrates **tree‑based survival modeling across regression, classification, survival, competing risks, geometric trees, and hazard forests**, providing a complete end‑to‑end ecosystem aligned with the ASA curriculum.

> **R‑only note:**  
> The ASA ecosystem is currently implemented only in R.  
> Python support will be added if and when official ASA packages become available.

---

# **Credits**

### **ASA Traveling Course Contributors**
- **Hemant Ishwaran, PhD** — University of Miami  
- **Min Lu, PhD** — University of Miami  

### **Source Materials**
- ASA Traveling Course Slide Guide  
- randomForestSRC (CRAN)  
- randomForestSRC.run (GitHub)  
- varPro, randomForestSGT, randomForestRHF  
