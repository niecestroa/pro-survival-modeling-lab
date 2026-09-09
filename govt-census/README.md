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
