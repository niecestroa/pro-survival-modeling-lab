# **ASA Traveling Course Master Workflow (Parts I–IV)**  
Use this as `README_ASA.md`.

---

# **ASA Traveling Course — Master Workflow (Parts I–IV)**  
### **GBCS Machine‑Learning Forest Pipeline — RandomForestSRC, VarPro, SGT, RHF**

This repository contains a unified master script integrating **all four modules** of the ASA Traveling Course: *Tree‑Based Machine Learning Methods*. The workflow is fully adapted to the **GBCS breast‑cancer dataset**, with every line of code annotated using medium‑length (S2) end‑of‑line comments for clarity and teaching value.

---

## **📦 Contents**
The master script is organized into five modules:

### **Module 0 — Header + Dataset Preparation**
- Loads required packages  
- Imports the GBCS dataset  
- Converts raw variables into ASA‑compatible formats  
- Creates:
  - `time`, `status` for survival forests  
  - `y_class` for classification forests  
  - `y_reg` for regression forests  
- Removes unsupported columns (IDs, dates, raw encodings)  
- Ensures no missing values  
- Produces `gbcs_asa`, the unified dataset used across all ASA modules  

---

### **Module 1 — ASA Part I: Training**
Implements all training examples from ASA Part I, including:

- Regression forests  
- Classification forests  
- Survival forests  
- Quantile regression forests  
- CART‑style single‑tree interface  
- SID clustering  
- Integrated analysis via `run.rfsrc`  

All examples are preserved exactly as in the ASA companion code, with S2 explanations added.

---

### **Module 2 — ASA Part II: Inference & Prediction**
Implements ASA inference and prediction topics:

- OOB classification inference  
- OOB survival inference  
- Survival curve extraction  
- Prediction on new data  
- Factor‑level mismatch handling  
- Restore‑mode prediction  
- Custom OOB‑weighted estimators  
- Partial dependence plots (mortality + survival)  

All examples are annotated and adapted to the GBCS dataset where appropriate.

---

### **Module 3 — ASA Part III: Variable Selection**
Implements ASA variable‑selection methods:

- Permutation VIMP  
- Block‑size VIMP  
- Joint VIMP  
- Subsampling inference + VIMP confidence intervals  
- Minimal depth  
- Guided trees using split weights  
- VarPro (supervised)  
- cv.varpro (cross‑validated)  
- ivarpro + shap.ivarpro  
- uvarpro (unsupervised)  

All examples are preserved and explained line‑by‑line.

---

### **Module 4 — ASA Part IV: Advanced Topics**
Implements ASA advanced topics:

- Imbalanced classification (RFQ, BRF)  
- G‑mean VIMP + subsampling  
- Missing‑data imputation (supervised + unsupervised)  
- missForest / mForest grouped regressions  
- impute.learn + OOD scoring  
- Super Greedy Trees (SGT)  
- hcut tuning  
- SGT model explainers (beta, partial effects)  
- Random Hazard Forests (RHF)  
- Counting‑process conversion  
- Time‑dependent AUC  
- Smoothed hazards  
- Time‑localized RHF variable importance  

All examples are annotated and adapted to GBCS where applicable.

---

### **Module 5 — Unified Orchestrator**
Provides a single callable function:

```r
run_master_asa()
```

This function:

- Loads the ASA‑ready dataset  
- Runs selected workflows from Parts I–IV  
- Returns a structured list containing:
  - dataset  
  - part1 results  
  - part2 results  
  - part3 results  
  - part4 results  

---

## **Citation**
This workflow incorporates examples from:

**ASA Traveling Course: Tree‑Based Machine Learning Methods**  
American Statistical Association  
Student R‑Code Companion

---

## **Purpose**
This master script serves as:

- A complete teaching resource  
- A reproducible machine‑learning pipeline  
- A unified reference for all ASA forest‑based methods  
- A fully annotated codebase for training, inference, variable selection, and advanced modeling  
