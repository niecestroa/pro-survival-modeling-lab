# German Breast Cancer Survival Analysis — Professional Pipeline  
# Cross‑Language Survival Analysis Workflow

## Overview
This project contains the **professional, pipeline‑driven implementation** of the German Breast Cancer survival analysis. The work began as a traditional academic project written entirely in **R**, then evolved into a **modular R pipeline**, and was later **translated into SAS and Python** to demonstrate cross‑platform reproducibility and industry‑aligned statistical programming.

The result is a fully engineered, multi‑language survival‑analysis workflow that mirrors real‑world clinical and regulatory environments.

---

## Project Evolution  
### **1. Academic Beginning — R (Exploratory Coursework)**
The project originally started as a graduate‑level survival analysis assignment written in **base R**:

- Sequential, line‑by‑line scripts  
- Manual model fitting  
- Stepwise exploration of Cox PH models  
- Spline checks, interactions, and Weibull modeling  
- Presentation‑oriented figures and diagnostics  

This version focused on learning statistical concepts, not engineering a reproducible workflow.

---

### **2. Professional Refactor — R Pipeline (Modular, Automated)**
After completing the academic work, the entire analysis was rebuilt into a **professional R pipeline**:

- Modular scripts (`01_load_data.R`, `02_clean_data.R`, etc.)  
- Automated execution via a pipeline runner  
- Tidyverse‑based data cleaning  
- purrr‑driven model loops  
- ggplot diagnostics and residual plots  
- Structured model selection (AIC/BIC, PH tests, interactions)  
- Reproducible outputs and standardized folder structure  

This version represents the first major step toward a production‑ready workflow.

---

### **3. SAS Translation — Clinical‑Grade Workflow**
The R pipeline was then **translated into SAS**, mirroring workflows used in clinical trials and regulatory submissions:

- `PROC IMPORT` and `DATA` steps for controlled data preparation  
- `PROC PHREG` for Cox modeling and influence diagnostics  
- `PROC LIFEREG` for Weibull regression  
- `PROC LIFETEST` for Kaplan–Meier estimation  
- SGPLOT‑based diagnostic visualizations  
- FitStatistics‑based AIC/BIC model comparison  
- Case‑deletion and DFbeta influence analysis  

The SAS version demonstrates the ability to implement survival analysis in a validated, audit‑ready environment.

---

### **4. Python Implementation — Reproducible, Scripted Analytics**
Finally, the pipeline was ported into **Python** using:

- `pandas` for data ingestion and cleaning  
- `lifelines` for Cox PH and Weibull models  
- `matplotlib` / `seaborn` for diagnostics  
- Automated model loops and reproducible reporting  

The Python version completes the cross‑language workflow and shows the ability to translate statistical logic across modern analytics ecosystems.

---

## Pipeline Components (Implemented in R → SAS → Python)
- **Data ingestion and formatting**  
- **Data cleaning and preprocessing**  
- **Kaplan–Meier estimation**  
- **Cox proportional hazards modeling**  
- **Spline, log, and quadratic transformations**  
- **Interaction screening and model selection**  
- **PH assumption diagnostics**  
- **DFbeta and case‑deletion influence analysis**  
- **Weibull regression and Cox–Weibull comparison**  
- **Automated figure and table generation**  

Each step is implemented consistently across all three languages to ensure analytical integrity.

---

## Intended Use
This professional version is designed for:

- Reproducible research workflows  
- Cross‑language validation (R → SAS → Python)  
- Demonstrating industry‑aligned statistical programming practices  
- Clinical‑style documentation and traceability  
- Portfolio presentation of production‑ready survival analysis  

---

## Notes
This pipeline is intentionally paired with the academic version to highlight the evolution from:

- **Exploratory, presentation‑oriented R scripts** →  
- **Modular R pipeline** →  
- **Clinical‑grade SAS workflow** →  
- **Reproducible Python implementation**

This progression demonstrates:

- mastery of survival analysis across platforms  
- ability to translate statistical logic into production environments  
- readiness for biostatistics, clinical research, and data science roles  

---

# **Cross‑Language Function Comparison Table**  
### *R → SAS → Python equivalents for your survival‑analysis pipeline*

This table shows how each major analytical step in your project maps across the three languages you used.

| Task / Concept | **R Function** | **SAS Procedure** | **Python (lifelines / pandas)** |
|----------------|----------------|-------------------|----------------------------------|
| Import CSV | `read_csv()` | `PROC IMPORT` | `pd.read_csv()` |
| Data cleaning | `mutate()`, `case_when()` | `DATA` step | `df.assign()`, `np.where()` |
| Factor variables | `factor()` | `FORMAT` + `CLASS` | `astype('category')` |
| Kaplan–Meier | `survfit()` | `PROC LIFETEST` | `KaplanMeierFitter()` |
| Cox PH model | `coxph()` | `PROC PHREG` | `CoxPHFitter()` |
| Martingale residuals | `residuals(type="martingale")` | `OUTPUT resmart=` | `CoxPHFitter().compute_residuals(..., "martingale")` |
| Schoenfeld PH test | `cox.zph()` | `ASSESS PH` | `check_assumptions()` |
| DFbeta influence | `residuals(type="dfbeta")` | `OUTPUT dfbeta=` | `compute_residuals(..., "dfbeta")` |
| Spline terms | `pspline()` | `EFFECT spl=Spline()` | `patsy.bs()` or custom spline basis |
| Interaction terms | `size*hormone_f` | `size*hormone_f` | `'size:hormone_f'` in formula |
| Weibull regression | `survreg(dist="weibull")` | `PROC LIFEREG` | `WeibullAFTFitter()` |
| KM plots | `ggplot2` | `SGPLOT` | `matplotlib` / `seaborn` |
| Model comparison | AIC/BIC via `extractAIC()` | FitStatistics table | `.AIC_`, `.BIC_` attributes |
| Pipeline automation | `purrr::map()` | SAS macros | Python loops / functions |
