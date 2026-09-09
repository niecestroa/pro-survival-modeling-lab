/*
Author: Aaron Niecestro
Created on July 24, 2024
Last Editted on February 14, 2026
*/

/*===============================================================
  Load Data
===============================================================*/

proc import datafile="C:\Users\aniec\UTH Survival Analysis\PH1831_Project\PH1831 Project Data\PBC276.csv"
    out=pbc_raw
    dbms=csv
    replace;
    guessingrows=max;
run;

/* Clean column names and convert variables */
data pbc;
    set pbc_raw;

    /* Convert to numeric if needed */
    age      = input(age, best12.);
    bili     = input(bili, best12.);
    copper   = input(copper, best12.);
    albumin  = input(albumin, best12.);
    protime  = input(protime, best12.);
    futime   = input(futime, best12.);
    status   = input(status, best12.);

    /* Convert to categorical */
    stage = put(stage, best12.);
    drug  = put(drug, best12.);
    sex   = put(sex, best12.);

    stage = strip(stage);
    drug  = strip(drug);
    sex   = strip(sex);

    /* Log transforms */
    log_bili   = log(bili);
    log_copper = log(copper);
run;

/*===============================================================
  Null Cox Model + Martingale Residuals
===============================================================*/

proc phreg data=pbc;
    model futime*status(0) = ;
    output out=pbc_resid resmart=martingale;
run;

/*===============================================================
  Martingale Residual Plots (Linear + Log Scale)
===============================================================*/

%macro mart_plot(var);
proc sgplot data=pbc_resid;
    scatter x=&var y=martingale;
    loess x=&var y=martingale / lineattrs=(color=blue thickness=2);
    title "Martingale Residuals vs &var";
run;
%mend;

%mart_plot(age);
%mart_plot(bili);
%mart_plot(copper);
%mart_plot(albumin);
%mart_plot(protime);

/* Log-scale versions */
%macro mart_plot_log(var);
proc sgplot data=pbc_resid;
    scatter x=&var y=martingale;
    loess x=&var y=martingale / lineattrs=(color=red thickness=2);
    xaxis type=log;
    title "Martingale Residuals vs log(&var)";
run;
%mend;

%mart_plot_log(age);
%mart_plot_log(bili);
%mart_plot_log(copper);
%mart_plot_log(albumin);
%mart_plot_log(protime);

/*===============================================================
  Spline Fits (PROC TRANSREG)
===============================================================*/

%macro spline(var);
proc transreg data=pbc;
    model identity(futime) = spline(&var);
run;
%mend;

%spline(age);
%spline(bili);
%spline(copper);
%spline(albumin);
%spline(protime);

/*===============================================================
  Full Cox Model + Stepwise Selection
===============================================================*/

proc phreg data=pbc;
    class stage drug sex / param=ref;
    model futime*status(0) =
        bili copper albumin protime age stage drug;
    selection=stepwise slentry=0.05 slstay=0.05;
run;

/*===============================================================
  Log-Transformed Model + Stepwise
===============================================================*/

proc phreg data=pbc;
    class stage drug sex / param=ref;
    model futime*status(0) =
        albumin protime age log_bili log_copper stage drug;
    selection=stepwise slentry=0.05 slstay=0.05;
run;

/*===============================================================
  Final Models
===============================================================*/

proc phreg data=pbc;
    class drug / param=ref;
    model futime*status(0) =
        albumin protime age log_bili log_copper drug;
run;

proc phreg data=pbc;
    model futime*status(0) =
        albumin protime age log_bili log_copper;
run;

/*===============================================================
  Proportional Hazards Assumption
===============================================================*/

proc phreg data=pbc plots(only)=cumhaz;
    class drug / param=ref;
    model futime*status(0) =
        albumin protime age log_bili log_copper drug;
    assess ph / resample;
run;

/*===============================================================
  DFbeta Residuals
===============================================================*/

proc phreg data=pbc;
    class drug / param=ref;
    model futime*status(0) =
        albumin protime age log_bili log_copper drug;
    id _n_;
    output out=dfb_resid dfbeta=dfb_albumin dfb_protime dfb_age dfb_logbili dfb_logcopper dfb_drug;
run;

/* Plot DFbeta for each coefficient */
%macro dfb_plot(var,label);
proc sgplot data=dfb_resid;
    scatter x=_n_ y=&var;
    refline 0 / axis=y lineattrs=(color=red);
    title "DFbeta for &label";
run;
%mend;

%dfb_plot(dfb_albumin, Albumin);
%dfb_plot(dfb_protime, Protime);
%dfb_plot(dfb_age, Age);
%dfb_plot(dfb_logbili, Log(Bili));
%dfb_plot(dfb_logcopper, Log(Copper));
%dfb_plot(dfb_drug, Drug);

/*===============================================================
  Deviance Residuals
===============================================================*/

proc phreg data=pbc;
    class drug / param=ref;
    model futime*status(0) =
        albumin protime age log_bili log_copper drug;
    output out=dev_resid resdev=dev;
run;

%macro dev_plot(var);
proc sgplot data=dev_resid;
    scatter x=&var y=dev / alpha=0.6;
    loess x=&var y=dev / lineattrs=(color=blue);
    title "Deviance Residuals vs &var";
run;
%mend;

%dev_plot(age);
%dev_plot(protime);
%dev_plot(albumin);
%dev_plot(log_bili);
%dev_plot(log_copper);

/*===============================================================
  Case-Deletion Plots (DFbeta proxy)
===============================================================*/

%macro case_del(var,label);
proc sgplot data=dfb_resid;
    needle x=_n_ y=&var / markers;
    refline 0 / axis=y;
    title "Case Deletion Plot for &label";
run;
%mend;

%case_del(dfb_albumin, Albumin);
%case_del(dfb_protime, Protime);
%case_del(dfb_age, Age);
%case_del(dfb_logbili, Log(Bili));
%case_del(dfb_logcopper, Log(Copper));
%case_del(dfb_drug, Drug);

/*===============================================================
  Survival Curves for Drug × Sex
===============================================================*/

proc lifetest data=pbc plots=survival;
    time futime*status(0);
    strata drug*sex;
run;
