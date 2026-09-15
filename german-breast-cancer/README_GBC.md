# **GBCS Master Workflow — ASA Forest Methods + Classical Survival Analysis**  
### **Unified Machine‑Learning + Statistical Survival Pipeline**

This repository contains a comprehensive master workflow integrating:

1. **ASA Traveling Course (Parts I–IV)**  
2. **Full classical survival‑analysis pipeline (Cox, Weibull, diagnostics)**  

Both workflows operate on the **GBCS breast‑cancer dataset**, with every line of code annotated using medium‑length (S2) end‑of‑line comments.

---

## **Contents**
The master workflow is divided into two major sections:

---

# **SECTION A — ASA Traveling Course (Parts I–IV)**  
*(See README_ASA.md for full details.)*

### Includes:
- Regression, classification, survival forests  
- Quantile forests  
- SID clustering  
- OOB inference  
- Prediction on new data  
- Partial dependence plots  
- Permutation VIMP  
- Minimal depth  
- VarPro (supervised + unsupervised)  
- Super Greedy Trees (SGT)  
- Random Hazard Forests (RHF)  
- Missing‑data imputation  
- OOD scoring  
- Unified orchestrator  

All ASA modules are fully annotated and adapted to the GBCS dataset.

---

# **SECTION B — Classical Survival Analysis Pipeline**  
Your full survival‑analysis workflow is integrated alongside ASA methods, including:

### **1. Data Preparation**
- Cleaning  
- Factor conversion  
- Numeric standardization  
- Median splits  
- Removal of unsupported fields  

### **2. Null Cox Model + Martingale Residuals**
- Functional‑form diagnostics  
- Nonlinearity detection  
- Transformation discovery  

### **3. Residual Diagnostics**
- Martingale residuals  
- Deviance residuals  
- Log‑scale diagnostics  
- Spline diagnostics  

### **4. Multivariable Model Comparison**
- Full vs reduced models  
- Sensitivity analysis  
- Confounder detection  

### **5. Reduced Models**
- Core predictors  
- Confounding checks  
- Effect‑size stability  

### **6. 20% Change‑in‑Estimate Rule**
- Confounder identification  
- Effect‑size comparison  

### **7. Transformation Checks**
- Log transforms  
- Quadratic terms  
- Splines  

### **8. Interaction Screening**
- Biologically plausible interactions  
- ANOVA comparison  
- Effect modification detection  

### **9. Stepwise AIC/BIC**
- Model selection  
- Penalized complexity control  

### **10. Final Cox Model**
- PH assumption testing  
- DFbeta influence diagnostics  
- Outlier detection  

### **11. Weibull Model**
- Parametric alternative  
- Hazard‑ratio comparison  
- Deviance residuals  

### **12. Survival Curves**
- Treatment‑effect visualization  
- Cox vs Weibull comparison  

### **13. Unified Survival Orchestrator**
A single function that:

- Runs the entire survival pipeline  
- Returns all models, diagnostics, and plots  
- Produces a structured output for reporting  

---

# **SECTION C — Combined Orchestrator**
The final orchestrator merges both ASA and classical survival workflows:

```r
run_master_all()
```

This function returns:

- ASA Part I–IV results  
- Classical survival‑analysis results  
- Unified dataset  
- All diagnostics  
- All models  
- All plots  

---

## **Citation**
ASA Traveling Course: *Tree‑Based Machine Learning Methods*  
American Statistical Association  
Student R‑Code Companion

---

## **Purpose**
This combined workflow is designed for:

- Teaching  
- Research  
- Clinical modeling  
- Machine‑learning benchmarking  
- Survival‑analysis reproducibility  
- Cross‑method comparison (Cox vs forests vs RHF vs SGT)  

It is the most complete and annotated version of the GBCS modeling pipeline.
