# Survival Analysis Project  

# Abstract
This project traces the full evolution of a survival‑analysis workflow from its origins in graduate‑level coursework to a fully engineered, cross‑language analytical pipeline. Beginning with exploratory R scripts focused on learning statistical concepts, the work was systematically refactored into a modular, automated, and reproducible R pipeline aligned with professional data‑science and clinical‑research standards. The pipeline was then translated into SAS to mirror validated, audit‑ready workflows used in pharmaceutical and regulatory environments, and subsequently implemented in Python to demonstrate modern, scriptable analytics and cross‑platform reproducibility. Across all implementations, the project performs standardized data preparation, Kaplan–Meier estimation, Cox proportional‑hazards modeling, Weibull regression, diagnostic evaluation, and structured model selection. The result is a transparent, maintainable, and multi‑language survival‑analysis framework that highlights both statistical expertise and engineering maturity.

---

### Academic → Professional → Cross‑Language Workflow Evolution

This project demonstrates my progression from **academic survival‑analysis coursework** to a **professional, modular, reproducible analysis pipeline** implemented across **R, SAS, and Python**.  

The work began entirely in **R**, using exploratory, sequential scripts typical of graduate‑level coursework. After completing the academic version, I rebuilt the entire analysis into a **professional R pipeline**, then translated that pipeline into **SAS** to mirror clinical‑research workflows, and finally implemented the full workflow in **Python** to demonstrate cross‑platform reproducibility.

Each version analyzes the same dataset and addresses the same statistical questions, but the **workflow, engineering discipline, and reproducibility** improve dramatically across stages.

---

# Academic Version  
### Exploratory, Sequential, Concept‑Focused (R)

The **academic** folder contains the original survival‑analysis work completed during graduate coursework. It includes:

- A full written report  
- Presentation slides  
- Line‑by‑line R scripts  
- Manual model fitting and diagnostics  
- Stepwise exploration of Cox PH models, splines, interactions, and Weibull models  

This version emphasizes **learning the methods**, understanding model behavior, and interpreting survival outcomes. It intentionally preserves the exploratory workflow typical of academic training.

---

# Professional Version — R  
### Modular, Automated, Reproducible Pipeline

The **professional R** folder contains a fully engineered survival‑analysis pipeline built after the academic work was completed. It restructures the entire analysis into a modern, production‑ready workflow:

### **Pipeline Features**
- Modular R scripts for each analytical step  
- Automated execution via `run_pipeline.R`  
- Reproducible outputs saved to standardized folders  
- ggplot‑based diagnostics and residual plots  
- Functional programming (`purrr`) for model loops  
- Clean model‑selection workflow (AIC/BIC, PH tests, interactions)  
- Case‑wise deletion diagnostics and influence analysis  
- Parallel Cox and Weibull model comparison  

This version demonstrates the ability to transform exploratory academic code into a **clinical‑grade R pipeline**.

---

# Professional Version — SAS  
### Clinical‑Style, Regulatory‑Aligned Workflow

The **professional SAS** version translates the R pipeline into a workflow aligned with clinical‑trial and regulatory environments:

- `PROC IMPORT` and `DATA` steps for controlled data preparation  
- `PROC PHREG` for Cox modeling and influence diagnostics  
- `PROC LIFEREG` for Weibull regression  
- `PROC LIFETEST` for Kaplan–Meier estimation  
- SGPLOT‑based diagnostic visualizations  
- FitStatistics‑based AIC/BIC model comparison  
- Case‑deletion and DFbeta influence analysis  

This version mirrors the structure and rigor expected in pharmaceutical and CRO settings.

---

# Professional Version — Python  
### Modern, Scriptable, Cross‑Validated Analytics

The **professional Python** version implements the entire workflow using:

- `pandas` for data ingestion and cleaning  
- `lifelines` for Cox PH, Weibull, and diagnostics  
- `matplotlib` / `seaborn` for visualization  
- Automated model loops and reproducible reporting  
- Survival‑curve comparison between Cox and Weibull models  

This version demonstrates the ability to translate statistical logic into modern analytics ecosystems.

---

# Cross‑Language Comparison  
### Academic R → Professional R → SAS → Python

| Analytical Component | **Academic R** | **Professional R** | **Professional SAS** | **Professional Python** |
|----------------------|----------------|---------------------|------------------------|--------------------------|
| **Data Import** | `read.csv()` | `readr::read_csv()` | `PROC IMPORT` | `pd.read_csv()` |
| **Data Cleaning** | `mutate()` inline | Modular cleaning script | `DATA` step | `df.assign()`, `np.where()` |
| **Factor Handling** | `factor()` | `forcats` utilities | `FORMAT` + `CLASS` | `astype('category')` |
| **KM Estimation** | `survfit()` | Modular KM script | `PROC LIFETEST` | `KaplanMeierFitter()` |
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

# Purpose of This Repository

This project serves as:

- A transparent record of my **growth** as a biostatistician and data scientist  
- A demonstration of how academic work can be refactored into **industry‑aligned pipelines**  
- A portfolio example of survival‑analysis modeling, diagnostics, and workflow engineering  
- A showcase of **cross‑language mastery** (R → SAS → Python)  
- A template for future projects following the same academic → professional → multi‑language structure  
