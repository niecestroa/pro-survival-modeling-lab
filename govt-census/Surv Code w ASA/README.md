# **Government Survival Code with ASA**

## **Overview**
This directory contains a fully‑modular survival‑analysis pipeline adapted from the **ASA Traveling Course: “Tree‑Based Machine Learning Methods”**, originally demonstrated on the **German PBC276 dataset**.  
Here, the workflow has been extended and modified to support **government census survival modeling**, including:

- Random Survival Forest (RSF) training  
- Inference & prediction workflows  
- Variable selection using VIMP, minimal depth, VarPro, and iVarPro  
- Advanced topics such as supervised imputation and out‑of‑distribution (OOD) detection  

The folder includes a **master orchestration script** plus **four ASA‑derived modules** that were customized for your census‑based survival dataset.

---

## **Folder Contents**
> *Based on your description: 1 master script + 4 modified ASA parts.*

| File | Description |
|------|-------------|
| **master_surv_driver.R** | Full pipeline orchestrator. Runs Parts I–IV end‑to‑end. Loads PBC276 or your census dataset. |
| **part1_training.R** | RSF training code adapted from ASA materials. Handles model fitting, hyperparameters, and baseline diagnostics. |
| **part2_inference.R** | Inference utilities: OOB curves, C‑index, Brier score, prediction examples, unseen factor‑level handling. |
| **part3_varselect.R** | Variable selection methods: VIMP, subsampled CI, minimal depth, VarPro, cross‑validated VarPro, iVarPro. |
| **part4_advanced.R** | Advanced topics: supervised imputation, OOD scoring, test‑time imputation, iterative prediction. |

If you want, I can rewrite these filenames to match your actual file names exactly — just paste them.

---

## **Master Script Summary**
The master script (the one you pasted) performs:

### **Part I — RSF Training**
- Loads and cleans PBC276 (or your census dataset)
- Fits a Random Survival Forest using `randomForestSRC`
- Uses logrank splitting, permutation importance, and 1000 trees

### **Part II — Inference & Prediction**
- Plots in‑bag vs OOB survival curves for selected individuals  
- Computes C‑index and Brier score  
- Predicts survival for first 10 patients  
- Demonstrates unseen factor‑level handling (e.g., `drug = "newdrug"`)

### **Part III — Variable Selection**
- VIMP extraction  
- Subsampled confidence intervals  
- Minimal depth ranking  
- VarPro importance and cross‑validation  
- iVarPro table for interpretable variable grouping

### **Part IV — Advanced Topics**
- Train/test split  
- Supervised imputation using `impute.learn`  
- Test‑time imputation  
- OOD scoring and percentile ranking

---

## **How to Run**
From R:

```r
source("master_surv_driver.R")
main()
```

Or run individual parts:

```r
fit_rsf <- run_part1_training(PBC276)
run_part2_inference(PBC276, fit_rsf)
run_part3_varselect(PBC276, fit_rsf)
run_part4_advanced(PBC276)
```

---

## **Dependencies**
The pipeline uses:

- **tidyverse** — data manipulation  
- **janitor** — column cleaning  
- **survival** — Surv objects  
- **randomForestSRC** — RSF models, VIMP, minimal depth  
- **varPro** — VarPro, cv.VarPro, iVarPro  
- **matplot/matlines** — base R plotting for survival curves  

---

## **Dataset Notes**
The original ASA course uses **PBC276**.  
Your modified version supports:

- Government census survival datasets  
- Custom factor levels  
- OOD detection for new demographic categories  
- Supervised imputation for missing socioeconomic variables  

If you want, I can add a section describing your census dataset specifically (variables, censoring definition, time‑to‑event construction).

---

## **Citation**
If you use this code in publications:

```
ASA Traveling Course: Tree-Based Machine Learning Methods (2024–2026)
RandomForestSRC: Ishwaran & Kogalur
VarPro: Functional Variable Selection for Survival Models
```
