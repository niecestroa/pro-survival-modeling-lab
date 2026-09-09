/*
Author: Aaron Niecestro
Created on January 5, 2024
Last Editted on February 15, 2026
*/

/*----------------------------------------------------------
  0. Setup: Libraries and Options
----------------------------------------------------------*/
options nodate nonumber;
ods graphics on;

libname proj "C:\Users\aniec\BIST0627 Applied Survival Data Analysis\BIST0627 Project\project_data_files\gbcs.csv";

/*----------------------------------------------------------
  1. Import and Prepare Data
----------------------------------------------------------*/
proc import datafile="C:\Users\aniec\BIST0627 Applied Survival Data Analysis\BIST0627 Project\project_data_files\gbcs.csv"
    out=gbcs_raw
    dbms=csv
    replace;
    guessingrows=max;
run;

data gbcs;
    set gbcs_raw;

    /* Numeric conversions */
    age        = input(age, best12.);
    size       = input(size, best12.);
    nodes      = input(nodes, best12.);
    prog_recp  = input(prog_recp, best12.);
    estrg_recp = input(estrg_recp, best12.);
    rectime    = input(rectime, best12.);
    survtime   = input(survtime, best12.);
    censdead   = input(censdead, best12.);

    /* Factors / formats */
    if menopause = "1" then menopause_f = "Yes";
    else menopause_f = "No";

    if hormone = "1" then hormone_f = "Yes";
    else hormone_f = "No";

    if censrec = "1" then censrec_f = "Recurrence";
    else censrec_f = "Censored";

    grade_f = grade;

    /* Median splits */
    if age > 53 then age_med_f = "Above(>53)";
    else age_med_f = "Below(<=53)";

    if rectime > 1084 then rectime_med_f = "Above";
    else rectime_med_f = "Below";

    /* Drop unused raw vars */
    drop id diagdateb recdate deathdate menopause hormone grade censrec;
run;

/* Optional formats */
proc format;
    value $ynfmt "Yes"="Yes" "No"="No";
    value $recfmt "Recurrence"="Recurrence" "Censored"="Censored";
run;

/*----------------------------------------------------------
  2. Null Cox Model + Martingale Residuals
----------------------------------------------------------*/
proc phreg data=gbcs;
    model survtime*censdead(0) = / ties=efron;
    output out=gbcs_resid
        resmart = martingale;
run;

/*----------------------------------------------------------
  3. Residual Diagnostics (Numeric) with SGPLOT
----------------------------------------------------------*/
%macro mart_plot_num(var=);
proc sgplot data=gbcs_resid;
    scatter x=&var y=martingale / markerattrs=(symbol=circlefilled);
    loess x=&var y=martingale / lineattrs=(color=blue);
    xaxis label="&var";
    yaxis label="Martingale Residuals";
    title "Martingale Residuals vs &var";
run;
%mend;

%mart_plot_num(var=age);
%mart_plot_num(var=size);
%mart_plot_num(var=nodes);
%mart_plot_num(var=prog_recp);
%mart_plot_num(var=estrg_recp);
%mart_plot_num(var=rectime);
%mart_plot_num(var=survtime);

/* Log-scale example for nodes */
proc sgplot data=gbcs_resid;
    scatter x=nodes y=martingale;
    loess x=nodes y=martingale / lineattrs=(color=red);
    xaxis type=log label="Nodes (log scale)";
    yaxis label="Martingale Residuals";
    title "Martingale Residuals vs log(nodes)";
run;

/*----------------------------------------------------------
  4. Multivariable Cox Models (Full / Reduced)
----------------------------------------------------------*/
/* Full-ish model example */
proc phreg data=gbcs;
    class menopause_f hormone_f grade_f censrec_f age_med_f rectime_med_f / param=ref;
    model survtime*censdead(0) =
        age size nodes prog_recp estrg_recp rectime
        menopause_f hormone_f grade_f censrec_f
        age_med_f rectime_med_f;
run;

/* Reduced model: size + prog_recp + rectime */
proc phreg data=gbcs;
    model survtime*censdead(0) = size prog_recp rectime;
    ods output ParameterEstimates=PE_reduced;
run;

/* Alternative reduced model with rectime_med_f */
proc phreg data=gbcs;
    class rectime_med_f / param=ref;
    model survtime*censdead(0) = size prog_recp rectime_med_f;
    ods output ParameterEstimates=PE_reduced_med;
run;

/*----------------------------------------------------------
  5. Interaction Models (Example)
----------------------------------------------------------*/
/* Base main-effects model */
proc phreg data=gbcs;
    class menopause_f hormone_f grade_f / param=ref;
    model survtime*censdead(0) = size prog_recp rectime;
    ods output FitStatistics=Fit_base;
run;

/* Example interaction: size*hormone_f */
proc phreg data=gbcs;
    class hormone_f / param=ref;
    model survtime*censdead(0) = size prog_recp rectime size*hormone_f;
    ods output FitStatistics=Fit_int_size_horm;
run;

/* You can compare -2 Log L or AIC/BIC from FitStatistics tables */

/* Final chosen interaction model (from your R work) */
proc phreg data=gbcs;
    class hormone_f age_med_f / param=ref;
    model survtime*censdead(0) =
        size prog_recp rectime hormone_f age_med_f
        prog_recp*rectime size*hormone_f;
    ods output ParameterEstimates=PE_final;
run;

/* PH assumption via ASSESS (approximate) */
proc phreg data=gbcs plots(only)=cumhaz;
    class hormone_f age_med_f / param=ref;
    model survtime*censdead(0) =
        size prog_recp rectime hormone_f age_med_f
        prog_recp*rectime size*hormone_f;
    assess ph / resample;
run;

/*----------------------------------------------------------
  6. DFbeta / Influence Diagnostics for Final Cox Model
----------------------------------------------------------*/
proc phreg data=gbcs;
    class hormone_f age_med_f / param=ref;
    model survtime*censdead(0) =
        size prog_recp rectime hormone_f age_med_f
        prog_recp*rectime size*hormone_f;
    id _n_;
    output out=cox_dfbeta
        dfbeta = dfb_size dfb_prog dfb_rectime dfb_horm dfb_age dfb_int1 dfb_int2;
run;

proc sgplot data=cox_dfbeta;
    series x=_n_ y=dfb_size;
    refline 0 / axis=y;
    xaxis label="Observation Index";
    yaxis label="DFbeta (Size)";
    title "DFbeta for Size in Final Cox Model";
run;

/*----------------------------------------------------------
  7. Weibull Model (Survreg Equivalent)
----------------------------------------------------------*/
proc lifereg data=gbcs;
    class hormone_f age_med_f / param=ref;
    model survtime*censdead(0) =
        size prog_recp rectime hormone_f age_med_f
        prog_recp*rectime size*hormone_f
        / dist=weibull;
    ods output ParameterEstimates=PE_weib;
run;

/* Compare Weibull (PH-transformed) vs Cox coefficients */
data coef_compare;
    merge PE_weib(where=(Parameter ne 'Intercept') rename=(Estimate=Weib_Est))
          PE_final(rename=(Estimate=Cox_Est));
    by Parameter;
run;

/*----------------------------------------------------------
  8. Deviance Residuals (Weibull) and Plots
----------------------------------------------------------*/
proc lifereg data=gbcs;
    class hormone_f age_med_f / param=ref;
    model survtime*censdead(0) =
        size prog_recp rectime hormone_f age_med_f
        prog_recp*rectime size*hormone_f
        / dist=weibull;
    output out=weib_resid
        resdev = dev_resid;
run;

proc sgplot data=weib_resid;
    vbox dev_resid / category=hormone_f;
    yaxis label="Deviance Residuals";
    xaxis label="Hormone Therapy";
    title "Deviance Residuals vs Hormone Therapy";
run;

proc sgplot data=weib_resid;
    vbox dev_resid / category=age_med_f;
    yaxis label="Deviance Residuals";
    xaxis label="Median Age Group";
    title "Deviance Residuals vs Median Age at Diagnosis";
run;

/*----------------------------------------------------------
  9. Case-Deletion (DFbeta) for Weibull Model (Approx)
----------------------------------------------------------*/
/* SAS doesn’t give dfbeta directly for LIFEREG like R,
   but you can approximate influence via custom resampling
   or focus on Cox DFbetas as primary. */

/*----------------------------------------------------------
  10. Hormone-Only Models: Cox vs Weibull
----------------------------------------------------------*/
proc phreg data=gbcs;
    class hormone_f / param=ref;
    model survtime*censdead(0) = hormone_f;
    baseline out=cox_horm_surv
        survival=SurvProb
        / covariates=hormone_f;
run;

/* For Weibull hormone-only */
proc lifereg data=gbcs;
    class hormone_f / param=ref;
    model survtime*censdead(0) = hormone_f / dist=weibull;
    ods output ParameterEstimates=PE_horm_weib;
run;

/* You can reconstruct parametric survival curves from PE_horm_weib
   or use PROC LIFETEST with strata=hormone_f for nonparametric curves. */

/* Example: Nonparametric KM by hormone */
proc lifetest data=gbcs plots=survival;
    time survtime*censdead(0);
    strata hormone_f;
run;
