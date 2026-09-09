# Government Census Survival Analysis — Professional Pipeline  

# Abstract
This project demonstrates the full lifecycle transformation of an academic survival‑analysis assignment into a production‑ready, cross‑language analytical pipeline suitable for clinical research, regulatory workflows, and modern data‑science environments. Beginning with exploratory R scripts focused on statistical learning, the work was systematically refactored into a modular R pipeline emphasizing automation, reproducibility, and engineering discipline. The workflow was then translated into SAS to align with validated, audit‑ready practices used in clinical trials, and subsequently implemented in Python to highlight contemporary, scriptable analytics. Across all three languages, the pipeline performs standardized data preparation, Kaplan–Meier estimation, Cox proportional‑hazards modeling, Weibull regression, diagnostic evaluation, and automated model selection. The result is a transparent, maintainable, and cross‑platform survival‑analysis framework that illustrates both statistical expertise and professional‑grade programming maturity.

---

### R -> SAS -> Python Cross‑Language Workflow

This folder contains the **production‑ready version** of the Government Census survival‑analysis project. The work began as a traditional academic assignment written entirely in **R**, then evolved into a **modular R pipeline**, and was later **translated into SAS** and **Python** to demonstrate cross‑platform reproducibility and industry‑aligned statistical programming.

This professional version showcases how exploratory coursework code can be transformed into a **clean, maintainable, and fully automated analysis pipeline** suitable for real‑world clinical research, data science, and regulatory environments.

---

# Project Evolution  
### **1. Academic Beginning — R (Exploratory Coursework)**  
The project originally started as a graduate‑level survival‑analysis assignment written in base R:

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
- ggplot‑based diagnostics  
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

This version demonstrates the ability to implement survival analysis in a validated, audit‑ready environment.

---

### **4. Python Implementation — Modern, Scriptable Analytics**  
Finally, the pipeline was ported into **Python** using:

- `pandas` for data ingestion and cleaning  
- `lifelines` for Cox PH and Weibull models  
- `matplotlib` / `seaborn` for diagnostics  
- Automated model loops and reproducible reporting  

The Python version completes the cross‑language workflow and shows the ability to translate statistical logic across modern analytics ecosystems.

---

# Cross‑Language Comparison  
### Academic R → Professional R → SAS → Python

| Analytical Component | **Academic R** | **Professional R** | **Professional SAS** | **Professional Python** |
|----------------------|----------------|---------------------|------------------------|--------------------------|
| **Data Import** | `read.csv()` | `readr::read_csv()` | `PROC IMPORT` | `pd.read_csv()` |
| **Data Cleaning** | Inline `mutate()` | Modular cleaning script | `DATA` step | `df.assign()`, `np.where()` |
| **Factor Handling** | `factor()` | `forcats` utilities | `FORMAT` + `CLASS` | `astype('category')` |
| **KM Estimation** | `survfit()` | Modular KM module | `PROC LIFETEST` | `KaplanMeierFitter()` |
| **Cox PH Model** | `coxph()` | Pipeline‑based modeling | `PROC PHREG` | `CoxPHFitter()` |
| **Residuals** | `residuals()` | ggplot diagnostics | `OUTPUT resmart=` | `compute_residuals()` |
| **PH Assumption** | `cox.zph()` | Automated PH checks | `ASSESS PH` | `proportional_hazard_test()` |
| **DFbeta Influence** | `residuals(type="dfbeta")` | Automated loops | `OUTPUT dfbeta=` | `compute_residuals("dfbeta")` |
| **Splines** | `pspline()` | Modular spline script | `EFFECT spl=Spline()` | `patsy.bs()` or custom |
| **Interactions** | `size*hormone_f` | purrr‑driven loops | `size*hormone_f` | Formula interactions |
| **Weibull Model** | `survreg()` | Pipeline module | `PROC LIFEREG` | `WeibullAFTFitter()` |
| **Model Selection** | Manual | AIC/BIC pipelines | FitStatistics | `.AIC_`, `.BIC_` |
| **Visualization** | Base R | ggplot2 | SGPLOT | matplotlib / seaborn |
| **Pipeline Automation** | None | `purrr::map()` | SAS macros | Python functions / loops |

This table highlights the engineering progression from exploratory R scripts to a fully cross‑validated, multi‑language survival‑analysis pipeline.

---

# Purpose of This Professional Folder

This professional version demonstrates:

- How academic R code can be refactored into **industry‑aligned pipelines**  
- How the same workflow can be implemented in **R, SAS, and Python**  
- Cross‑language mastery of survival‑analysis modeling  
- Reproducible, auditable, and transparent statistical programming  
- Readiness for clinical, regulatory, and data‑science environments  

This folder serves as a **portfolio example** of both statistical expertise and engineering maturity.
