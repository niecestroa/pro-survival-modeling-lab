###############################################################################
# Title: GBCS Survival Analysis + ASA RandomForestSRC Master Script (Parts I–IV)
# Author:       Aaron Niecestro
# Created on:   September 15, 2024
# Last Edit:    September 15, 2026
#
# Description:
# This master script integrates all four ASA Traveling Course modules
# (Training, Inference & Prediction, Variable Selection, Advanced Topics)
# into a single unified workflow using the GBCS breast‑cancer dataset.
# Every line includes a medium‑length descriptive explanation to support
# teaching, reproducibility, and transparent machine‑learning analysis.
#
# Citation:
# Portions of this workflow reference the ASA Traveling Course:
# “Tree-Based Machine Learning Methods,” American Statistical Association.
# Code examples adapted from the ASA student R‑code companion.
###############################################################################

library(tidyverse)    # Loads tidyverse for data manipulation, piping, and general preprocessing.
library(survival)     # Loads survival for Surv() objects and survival‑model compatibility.
library(randomForestSRC)    # Loads randomForestSRC for all forest‑based methods used in ASA modules.
library(varPro)       # Loads VarPro for supervised and unsupervised variable‑priority methods.
library(randomForestSGT)    # Loads Super Greedy Trees for advanced splitting and model explainers.
library(randomForestRHF)    # Loads Random Hazard Forests for counting‑process survival modeling.

# Load original GBCS dataset from your project directory.
gbcs <- read_csv("~/BIST0627 Applied Survival Data Analysis/BIST0627 Project/project_data_files/gbcs.csv")    # Reads the raw GBCS dataset into memory for preprocessing.

# Prepare ASA‑compatible version of the dataset.
gbcs_asa <- gbcs %>%
  mutate(
    time = as.numeric(survtime),    # Converts survival time to numeric and renames it to 'time' for ASA survival forests.
    status = case_when(
      censdead == 1 ~ 1,            # Assigns event indicator 1 when death occurred for survival modeling.
      censdead == 0 ~ 0,            # Assigns censoring indicator 0 when no death occurred.
      TRUE ~ 0                      # Provides a fallback value to ensure status is always defined.
    ),
    y_class = factor(if_else(censdead == 1, "Dead", "Alive")),    # Creates a binary classification target for ASA classification examples.
    y_reg = as.numeric(survtime),    # Creates a regression target using survival time for ASA regression examples.
    age = as.numeric(age),           # Ensures age is numeric for modeling and forest algorithms.
    size = as.numeric(size),         # Ensures tumor size is numeric for regression and survival forests.
    nodes = as.numeric(nodes),       # Ensures lymph‑node count is numeric for modeling.
    prog_recp = as.numeric(prog_recp),    # Converts progesterone receptor count to numeric for modeling.
    estrg_recp = as.numeric(estrg_recp),  # Converts estrogen receptor count to numeric for modeling.
    rectime = as.numeric(rectime),        # Converts recurrence time to numeric for modeling.
    menopause_f = factor(if_else(menopause == "1", "Yes", "No")),    # Converts menopause indicator to a factor for classification and survival forests.
    hormone_f = factor(if_else(hormone == "1", "Yes", "No")),        # Converts hormone‑therapy indicator to a factor for modeling.
    grade_f = factor(grade),             # Converts tumor grade to a factor for classification and survival forests.
    censrec_f = factor(if_else(censrec == "1", "Recurrence", "Censored"))    # Converts recurrence indicator to a factor for modeling.
  ) %>%
  select(
    -id, -diagdateb, -recdate, -deathdate,    # Removes ID and date fields that ASA forests cannot process.
    -menopause, -hormone, -grade, -censrec     # Removes raw categorical encodings now replaced by factor versions.
  ) %>%
  drop_na()    # Removes rows with missing values because ASA forest methods require complete data.

str(gbcs_asa)    # Displays the structure of the ASA‑ready dataset to confirm numeric and factor types.
anyNA(gbcs_asa)  # Checks for remaining missing values to ensure full compatibility with ASA methods.

###############################################################################
# MODULE 1 — ASA Traveling Course Part I: Training
###############################################################################

library(randomForestSRC)    # Loads the randomForestSRC package which provides all forest-based training functions used in Part I.
library(survival)           # Loads the survival package to support Surv() objects required for survival forests.

# Slide 6: Brief Overview
# These lines show example syntax used in the ASA slides.
# They are kept as comments because they illustrate the interface rather than being executed.
# rfsrc(Surv(time, status)~., data = veteran)
# rfsrc(Surv(time, status)~., data = wihs)
# rfsrc(Ozone~., data = airquality)
# quantreg(mpg~., data = mtcars)
# rfsrc(Species~., data = iris)
# imbalanced(status~., data = breast)
# rfsrc(Multivar(mpg, cyl)~., data = mtcars)
# rfsrc(cbind(Species,Sepal.Length)~.,data=iris)
# quantreg(cbind(mpg, cyl)~., data = mtcars)
# quantreg(cbind(Species,Sepal.Length)~.,data=iris)
# rfsrc(data = mtcars)
# sidClustering(data = mtcars)
# sidClustering(data = mtcars, method = "sh")

# Slide 7: Quick Start — Iowa Housing
data(housing, package = "randomForestSRC")    # Loads the Iowa housing dataset used in ASA training examples.
dim(housing)    # Displays the dimensions of the housing dataset to show its size and structure.

# Slide 9: Quick Start
library(randomForestSRC)    # Reloads randomForestSRC to ensure the namespace is active for training examples.
o <- rfsrc(SalePrice ~ ., data = housing)    # Trains a random forest to predict SalePrice using all housing predictors.

# Slide 10: Quick Start
library(randomForestSRC)    # Ensures randomForestSRC is loaded before running the forest again.
o <- rfsrc(SalePrice ~ ., data = housing)    # Fits the same random forest model to demonstrate printed output.
print(o)    # Prints the forest summary including OOB error and performance metrics.

# Slide 11: Quick Start — OOB R-squared emphasis
print(o)    # Prints the forest again to highlight OOB R-squared and error lines shown in the ASA slide.

# Slide 12: Quick Start — Log-transform example
o <- rfsrc(SalePrice ~ ., data = housing)    # Fits the forest using raw SalePrice values.
o    # Displays the forest summary for the untransformed target.
housing$SalePrice <- log(housing$SalePrice)    # Applies a log transform to SalePrice to stabilize variance.
o <- rfsrc(SalePrice ~ ., data = housing)    # Fits the forest using log-transformed SalePrice.
o    # Displays the forest summary for the transformed target.

# Slide 13: Quick Start — ntree example
data(housing, package = "randomForestSRC")    # Reloads the housing dataset to reset SalePrice.
housing$SalePrice <- log(housing$SalePrice)    # Applies log transform again for consistency.
o <- rfsrc(SalePrice ~ ., data = housing)    # Fits the forest using default ntree settings.
print(o)    # Prints the forest summary to show default ntree behavior.
# rfsrc(..., ntree = 500)    # Illustrative syntax showing how to increase the number of trees.

# Slide 14: nodesize example
# rfsrc(..., nodesize = 5)    # Illustrative syntax showing how to adjust terminal node size.

# Slide 15: mtry example
# rfsrc(..., mtry = NULL)    # Illustrative syntax showing how to modify the number of variables tried at each split.

# Slide 16: samptype example
# rfsrc(..., samptype = "swor")    # Illustrative syntax showing sampling without replacement.

# Slide 17: Quick Start — resampling emphasis
data(housing, package = "randomForestSRC")    # Reloads housing dataset to demonstrate resampling details.
housing$SalePrice <- log(housing$SalePrice)    # Applies log transform again for consistency.
o <- rfsrc(SalePrice ~ ., data = housing)    # Fits the forest to show resampling method in printed output.
print(o)    # Prints the forest summary emphasizing resample size and method.

# Slide 18: splitrule example
# rfsrc(..., splitrule = "mse")    # Illustrative syntax showing how to use mean squared error split rule.

# Slide 19: nsplit example
# rfsrc(..., nsplit = 10)    # Illustrative syntax showing how to specify number of random splits per node.

# Slide 20: Prediction example
o.pred <- predict(o, newdata = housing[c(1:10),])    # Predicts SalePrice for the first 10 housing cases.
head(o.pred$predicted)    # Displays the first few predicted values for inspection.

# Slide 22: CART-style interface
# rfsrc.cart(formula, data, ntree = 1, mtry = ncol(data), bootstrap = "none", ...)    # Illustrative syntax for single-tree CART mode.

# Slide 23: Nonparametric regression examples
# rfsrc(Ozone~., data = airquality)    # Example of regression forest on airquality dataset.
# rfsrc(Species~., data = iris)        # Example of classification forest on iris dataset.
# imbalanced(status~., data = breast)  # Example of imbalanced classification.
# rfsrc(Surv(time, status)~., data = veteran)    # Example of survival forest.
# rfsrc(Surv(time, status)~., data = wihs)       # Example of survival forest on WIHS dataset.
# rfsrc(Multivar(mpg, cyl)~., data = mtcars)     # Example of multivariate forest.
# rfsrc(cbind(Species,Sepal.Length)~.,data=iris) # Example of mixed outcome forest.
# quantreg(cbind(mpg, cyl)~., data = mtcars)     # Example of quantile regression forest.
# quantreg(cbind(Species,Sepal.Length)~.,data=iris) # Example of quantile forest with mixed outcomes.
# rfsrc(data = mtcars)    # Example of unsupervised forest.
# sidClustering(data = mtcars)    # Example of SID clustering.
# sidClustering(data = mtcars, method = "sh")    # Example of SID clustering with alternative method.

# Slide 25: Regression example — quantile regression
o <- quantreg(SalePrice ~ ., housing, splitrule = "mse", ntree = 250)    # Fits quantile regression forest using MSE split rule.
o <- quantreg(SalePrice ~ ., housing, splitrule = "quantile.regr", ntree = 250)    # Fits quantile regression forest using quantile split rule.
o <- quantreg(SalePrice ~ ., housing, splitrule = "la.quantile.regr", ntree = 250)    # Fits quantile regression forest using default LA quantile rule.
o    # Prints the quantile regression forest summary.
plot.quantreg(o)    # Plots quantile regression forest results for visualization.

# Slide 28: Optional installation
# install.packages("devtools")    # Shows how to install devtools if needed.
# devtools::install_github("kogalur/randomForestSRC.run")    # Shows how to install randomForestSRC.run from GitHub.

# Slide 29: run.rfsrc example
library(randomForestSRC.run)    # Loads the run.rfsrc extension for integrated forest analysis.
run.rfsrc(SalePrice ~ ., housing, ntree = 500)    # Runs integrated forest analysis with 500 trees.

# Slide 31: Classification examples
# rfsrc(Ozone~., data = airquality)    # Regression example.
# quantreg(mpg~., data = mtcars)       # Quantile regression example.
# imbalanced(status~., data = breast)  # Imbalanced classification example.
# rfsrc(Surv(time, status)~., data = veteran)    # Survival forest example.
# rfsrc(Surv(time, status)~., data = wihs)       # Survival forest example.
# rfsrc(Multivar(mpg, cyl)~., data = mtcars)     # Multivariate forest example.
# rfsrc(cbind(Species,Sepal.Length)~.,data=iris) # Mixed outcome forest example.
# quantreg(cbind(mpg, cyl)~., data = mtcars)     # Quantile regression example.
# quantreg(cbind(Species,Sepal.Length)~.,data=iris) # Quantile regression example.
# rfsrc(data = mtcars)    # Unsupervised forest example.
# sidClustering(data = mtcars)    # SID clustering example.
# sidClustering(data = mtcars, method = "sh")    # SID clustering example with alternative method.

# Slide 32: Glioma classification example
library(varPro)    # Loads VarPro for variable priority methods used later.
data(glioma, package = "varPro")    # Loads glioma dataset for classification example.
dim(glioma)    # Displays dimensions of glioma dataset to show its size.

# Slide 33: Glioma classification forest
o <- rfsrc(y ~ ., data = glioma)    # Fits a random forest classifier to predict glioma subtype.
o    # Prints the forest summary including OOB error.

# Slide 35: Glioma classification with Gini split rule
o <- rfsrc(y ~ ., data = glioma, splitrule = "gini")    # Fits classifier using Gini impurity split rule.
o    # Prints the forest summary.

# Slide 36: Glioma classification with AUC split rule
o <- rfsrc(y ~ ., data = glioma, splitrule = "auc")    # Fits classifier using AUC-based split rule.
o    # Prints the forest summary.

# Slide 37: Glioma classification with entropy split rule
o <- rfsrc(y ~ ., data = glioma, splitrule = "entropy")    # Fits classifier using entropy-based split rule.
o    # Prints the forest summary.

# Slide 38: run.rfsrc integrated analysis
run.rfsrc(y ~ ., data = glioma)    # Runs integrated forest analysis on glioma dataset.

# Slide 40: Survival example — PBC Mayo Clinic
# rfsrc(Ozone~., data = airquality)    # Regression example.
# quantreg(mpg~., data = mtcars)       # Quantile regression example.
# rfsrc(Species~., data = iris)        # Classification example.
# imbalanced(status~., data = breast)  # Imbalanced classification example.
# rfsrc(Surv(time, status)~., data = wihs)    # Survival forest example.
# rfsrc(Multivar(mpg, cyl)~., data = mtcars)  # Multivariate forest example.
# rfsrc(cbind(Species,Sepal.Length)~.,data=iris) # Mixed outcome forest example.
# quantreg(cbind(mpg, cyl)~., data = mtcars)     # Quantile regression example.
# quantreg(cbind(Species,Sepal.Length)~.,data=iris) # Quantile regression example.
# rfsrc(data = mtcars)    # Unsupervised forest example.
# sidClustering(data = mtcars)    # SID clustering example.
# sidClustering(data = mtcars, method = "sh")    # SID clustering example with alternative method.

# Slide 41: PBC survival example
data(pbc, package = "survival")    # Loads PBC dataset for survival forest example.
dim(pbc)    # Displays dimensions of PBC dataset.

# Slide 42: PBC survival forest
pbc$id <- NULL    # Removes ID column because forests cannot use identifiers as predictors.
pbc.cr <- pbc    # Stores competing-risk version for later use.
pbc$status[pbc$status > 0] <- 1    # Converts competing-risk status to binary death indicator.
o <- rfsrc(Surv(time, status) ~ ., data = pbc)    # Fits survival forest using right-censored death indicator.
o    # Prints the forest summary.

# Slide 44: PBC survival forest with logrank split rule
o <- rfsrc(Surv(time, status) ~ ., data = pbc, splitrule = "logrank")    # Fits survival forest using logrank split rule.
o    # Prints the forest summary.

# Slide 45: PBC survival forest with bs.gradient split rule
o <- rfsrc(Surv(time, status) ~ ., data = pbc, splitrule = "bs.gradient")    # Fits survival forest using gradient-based split rule.
o    # Prints the forest summary.

# Slide 46: PBC survival forest with logrankscore split rule
o <- rfsrc(Surv(time, status) ~ ., data = pbc, splitrule = "logrankscore")    # Fits survival forest using logrank score split rule.
o    # Prints the forest summary.

# Slide 47: run.rfsrc integrated survival analysis
run.rfsrc(Surv(time, status) ~ ., data = pbc)    # Runs integrated survival forest analysis on PBC dataset.

###############################################################################
# MODULE 2 — ASA Traveling Course Part II: Inference & Prediction
###############################################################################

library(randomForestSRC)    # Loads randomForestSRC to access forest-based inference and prediction functions.
library(survival)           # Loads survival to support Surv() objects used in survival inference.

# Slide 5: Key quantities for classification
# These are comments showing what the ASA forest object returns.
# o$predicted     --->   inbag estimated probabilities
# o$predicted.oob --->   OOB estimated probabilities
# o$class         --->   inbag class predictions
# o$class.oob     --->   OOB class predictions

# Slide 6: OOB classification example — Glioma
data(glioma, package = "varPro")    # Loads the glioma dataset used for classification inference examples.
o <- rfsrc(y ~ ., data = glioma)    # Fits a random forest classifier to predict glioma subtype.
print(o)    # Prints the forest summary including OOB error and class probabilities.

mean(o$class != o$yvar)    # Computes inbag misclassification rate for the fitted classifier.
mean(o$class.oob != o$yvar)    # Computes OOB misclassification rate to assess generalization performance.

# Slide 7: Key quantities — classification
o$predicted[1:5, ]    # Shows the first five rows of inbag class-probability estimates.
o$predicted.oob[1:5, ]    # Shows the first five rows of OOB class-probability estimates.
o$class[1:5]    # Displays the first five inbag class predictions.
o$class.oob[1:5]    # Displays the first five OOB class predictions.

# Slide 8: Key quantities for survival
# These comments describe survival-specific outputs returned by rfsrc().
# o$time.interest ---> event times used for survival estimation
# o$predicted     ---> inbag estimated mortality
# o$predicted.oob ---> OOB estimated mortality
# o$survival      ---> inbag survival estimator for each case
# o$survival.oob  ---> OOB survival estimator for each case
# o$chf           ---> inbag cumulative hazard function
# o$chf.oob       ---> OOB cumulative hazard function

# Slide 9: OOB survival example — PBC Mayo Clinic
data(pbc, package="survival")    # Loads the PBC dataset for survival inference examples.
pbc$id <- NULL    # Removes ID column because forests cannot use identifiers as predictors.
pbc$status[pbc$status > 0] <- 1    # Converts competing-risk status to binary death indicator for survival modeling.

o <- rfsrc(Surv(time, status)~., pbc)    # Fits a survival forest using right-censored death indicator.
idx <- c(11,34,60)    # Selects three cases for plotting survival curves.

matplot(o$time.interest, t(o$survival[idx,]), type = "l", col=4, lwd=3,
        xlab = "Days", ylab = "Survival")    # Plots inbag survival curves for selected cases.
matlines(o$time.interest, t(o$survival.oob[idx,]), type = "l", col=2, lwd=3)    # Adds OOB survival curves for comparison.
legend("bottomleft", legend = c("inbag", "oob"), fill = c(4,2))    # Adds legend distinguishing inbag vs OOB curves.

# Slide 10: Key quantities — survival
o$time.interest[1:5]    # Displays first five event times used for survival estimation.
o$predicted[1:5]        # Shows first five inbag mortality estimates.
o$predicted.oob[1:5]    # Shows first five OOB mortality estimates.

dim(o$yvar)             # Displays dimensions of the response variable matrix.
length(o$time.interest) # Shows number of unique event times.
dim(o$survival)         # Displays dimensions of inbag survival matrix.
dim(o$survival.oob)     # Displays dimensions of OOB survival matrix.

# Slide 12: Prediction error for classification
# Comments showing helper functions discussed in ASA slides.
# get.misclass.error(object)
# get.brier.error(object)
# get.logloss(object)
# get.auc(object)

# Slide 13: Classification example — Glioma
o.glioma <- rfsrc(y ~ ., data = glioma)    # Re-fits glioma classifier because 'o' was reused for PBC.
tail(o.glioma$err.rate, 1)    # Displays final OOB error rate for the glioma classifier.

# Slide 14: Prediction error for survival
# Comments showing helper functions for survival performance.
# object$err.rate or get.cindex(object)
# get.brier.survival(object)
# plot.brier.auc(object)

# Slide 16: Prediction — syntax examples
# predict(object, testdata)
# predict(object, testdata, outcome = "test")
# predict(object, ...)

# Slide 18: Prediction — canonical example
data(veteran, package = "randomForestSRC")    # Loads veteran dataset for prediction examples.
dim(veteran)    # Displays dimensions of veteran dataset.

# Slide 19: Prediction — canonical example (train/test split)
data(veteran, package = "randomForestSRC")    # Reloads veteran dataset to ensure clean state.
veteran2 <- data.frame(lapply(veteran, factor))    # Converts all predictors to factors for demonstration.
veteran2$time <- veteran$time    # Restores numeric survival time column.
veteran2$status <- veteran$status    # Restores numeric event indicator column.

train <- sample(1:nrow(veteran2), round(nrow(veteran2) * .25))    # Creates 25% training sample for prediction demonstration.

summary(veteran2[train,])    # Summarizes training subset to show factor levels.
summary(veteran2[-train,])   # Summarizes test subset to show factor levels.

# Slide 20: Prediction — canonical example (train → test)
o <- rfsrc(Surv(time, status) ~ ., veteran2[train, ])    # Trains survival forest on training subset.
pred <- predict(o, veteran2[-train , ])    # Predicts survival outcomes on test subset.

print(o)    # Prints summary of trained forest.
print(pred) # Prints summary of predictions on test data.

# Slide 21: Prediction — factor level mismatch
veteran3 <- veteran2[1:3, ]    # Selects three cases for demonstration.
veteran3$celltype <- factor(c("newlevel", "1", "3"))    # Introduces unseen factor level to test prediction robustness.
pred2 <- predict(o, veteran3)    # Predicts survival outcomes despite unseen factor level.
print(pred2)    # Prints prediction results including handling of unseen level.
print(pred2$xvar)    # Shows how unseen level is treated internally.

# Slide 22: Restore mode — syntax
# predict(object, ...)

# Slide 23: Restore mode — brier example
o <- rfsrc(y~., data = glioma)    # Fits glioma classifier for restore-mode demonstration.
p <- predict(o, perf.type = "brier")    # Computes Brier score using restore-mode prediction.
p    # Prints restore-mode prediction object.

# Slide 24: Restore mode — custom regression estimator
o <- rfsrc(mpg~.,mtcars)    # Fits regression forest on mtcars dataset.
fwt <- predict(o, forest.wt="oob")$forest.wt    # Extracts OOB forest weights for custom estimator.
yhat <- c(fwt %*% o$yvar)    # Computes custom weighted prediction using OOB weights.

print(summary(yhat - o$predicted.oob))    # Compares custom predictions to OOB ensemble predictions.
print(o)    # Prints forest summary for reference.
print(mean((yhat - o$yvar)^2))    # Computes custom prediction error for comparison.

# Slide 28: Survival example — peakVO2
data(peakVO2, package = "randomForestSRC")    # Loads peakVO2 dataset for survival prediction examples.
dim(peakVO2)    # Displays dimensions of peakVO2 dataset.

# Slide 29: Partial effect for age
o <- rfsrc(Surv(ttodead,died)~., peakVO2)    # Fits survival forest on peakVO2 dataset.
plot.variable(o, surv.type = "surv", xvar.names = "age",
              time = 2, partial = TRUE)    # Plots partial survival effect of age at time=2.

# Slide 30: Partial effect with smoothing
plot.variable(o, surv.type = "surv", xvar.names = "age",
              smooth.lines = TRUE,
              time = 2, partial = TRUE)    # Plots smoothed partial survival effect of age.

# Slide 31: Relative-frequency partial effect
plot.variable(o, surv.type = "rel.freq", xvar.names = "age",
              smooth.lines = TRUE,
              partial = TRUE)    # Plots partial effect using relative-frequency survival measure.

# Slide 32: Partial effect — peak VO2
partial.o <- partial(o,
                     partial.type = "mort",
                     partial.xvar = "peak.vo2",
                     partial.values = o$xvar$peak.vo2,
                     partial.time = o$time.interest)    # Computes partial mortality effect for peak VO2.
pdta.m <- get.partial.plot.data(partial.o)    # Extracts partial-plot data for mortality.

pvo2 <- quantile(o$xvar$peak.vo2)    # Computes quantiles of peak VO2 for survival partial plots.
partial.o <- partial(o,
                     partial.type = "surv",
                     partial.xvar = "peak.vo2",
                     partial.values = pvo2,
                     partial.time = o$time.interest)    # Computes partial survival effect for peak VO2.
pdta.s <- get.partial.plot.data(partial.o)    # Extracts partial-plot data for survival.

par(mfrow=c(1,2))    # Sets plotting layout to two side-by-side panels.

plot(lowess(pdta.m$x, pdta.m$yhat, f = 2/3),
     type = "l", xlab = "peak VO2", ylab = "adjusted mortality")    # Plots smoothed partial mortality effect.
rug(o$xvar$peak.vo2)    # Adds rug plot showing distribution of peak VO2 values.

matplot(pdta.s$partial.time, t(pdta.s$yhat), type = "l", lty = 1,
        xlab = "years", ylab = "peak VO2 adjusted survival")    # Plots partial survival curves for peak VO2 quantiles.
legend("bottomleft", legend = paste0("peak VO2 = ", pvo2),
       bty = "n", cex = .75, fill = 1:5)    # Adds legend showing peak VO2 quantile values.

###############################################################################
# MODULE 3 — ASA Traveling Course Part III: Variable Selection
###############################################################################

library(randomForestSRC)    # Loads randomForestSRC to access permutation VIMP, minimal depth, and guided-tree methods.
library(survival)           # Loads survival to support Surv() objects used in survival VIMP and VarPro examples.
library(varPro)             # Loads VarPro for supervised and unsupervised variable-priority methods.

# Slide 6: Different VIMP types
# importance = c("anti", "permute", "random")    # Shows available VIMP types in randomForestSRC.
# importance = TRUE       --> anti-VIMP           # Indicates anti-VIMP is used when importance=TRUE.
# importance = "permute"  --> Breiman-Cutler      # Indicates permutation VIMP is used when importance="permute".
# importance = "random"   --> random-VIMP         # Indicates random VIMP is used when importance="random".

# Slide 7: Obtaining VIMP using the package
vimp.grow <- rfsrc(mpg ~ ., data = mtcars, importance = "permute")$importance    # Computes permutation VIMP during forest growth.
vimp.grow.block <- rfsrc(mpg ~ ., data = mtcars, importance = "permute", block.size = 10)$importance    # Computes block-size permutation VIMP.

obj <- rfsrc(mpg ~ ., data = mtcars)    # Fits a regression forest to demonstrate restore-mode VIMP.
vimp.restore <- predict(obj, importance = "permute")$importance    # Computes permutation VIMP using restore-mode prediction.
vimp.restore.block <- predict(obj, importance = "permute", block.size = 10)$importance    # Computes block-size VIMP using restore-mode prediction.

vimp(obj, importance = "permute")    # Computes permutation VIMP using the dedicated vimp() interface.
vimp(obj, importance = "permute", block.size = 10)$importance    # Computes block-size permutation VIMP using vimp().

iris.obj <- rfsrc(Species ~ ., data = iris)    # Fits a classification forest on iris dataset for joint VIMP demonstration.
vimp(iris.obj, iris.obj$xvar.names[1:2], importance = "permute", joint = TRUE)$importance    # Computes joint VIMP for first two iris predictors.
vimp(iris.obj, iris.obj$xvar.names[3:4], importance = "permute", joint = TRUE)$importance    # Computes joint VIMP for last two iris predictors.

# Slide 8: General call to vimp
iris.obj <- rfsrc(Species ~ ., data = iris)    # Fits iris classifier again for VIMP demonstration.
print(vimp(iris.obj, importance = "permute")$importance)    # Prints permutation VIMP for all predictors.
print(vimp(iris.obj, c("Petal.Length", "Petal.Width"), joint = TRUE, importance = "permute")$importance)    # Prints joint VIMP for two predictors.

# Slide 11: Confidence intervals for VIMP
library(varPro)    # Loads VarPro to access subsampling-based VIMP confidence intervals.

data(peakVO2, package = "randomForestSRC")    # Loads peakVO2 dataset for survival VIMP demonstration.
o <- rfsrc(Surv(ttodead, died)~., peakVO2, importance="permute")    # Fits survival forest with permutation VIMP.
oo <- subsample(o)    # Performs subsampling to compute VIMP confidence intervals.

o$importance    # Displays raw permutation VIMP values.
plot.vimp.ci(oo, alpha=.05)    # Plots 95% confidence intervals for VIMP.

# Slide 18: Minimal depth illustration
md <- max.subtree(o)$order[, 1]    # Extracts minimal depth values for each predictor.
barplot(sort(md), las=2, horiz = TRUE, col = "cadetblue3")    # Plots sorted minimal depth values horizontally.

# Slide 19: Guided trees using split weights
xvar.used <- predict(o, var.used="all.trees")$var.used    # Extracts variable usage counts across all trees.
os <- rfsrc(Surv(ttodead, died)~., peakVO2, xvar.wt = xvar.used)    # Fits guided survival forest using variable usage as split weights.
mds <- max.subtree(os)$order[, 1]    # Computes minimal depth for guided forest.
barplot(sort(mds), las=2, horiz = TRUE, col = "cadetblue3")    # Plots guided minimal depth values.

# Slide 21: VarPro motivation
summary(peakVO2[,c("bun","interval", "peak.vo2")])    # Summarizes selected peakVO2 predictors to motivate VarPro usage.

# Slide 26: VarPro canonical illustration
o <- varpro(Surv(ttodead, died) ~ ., peakVO2)    # Fits supervised VarPro model for variable priority.
importance(o)    # Displays variable importance from VarPro.

o.cv <- cv.varpro(Surv(ttodead, died) ~ ., peakVO2)    # Performs cross-validated VarPro to select cutoff.
print(o.cv)    # Prints cross-validation results including cutoff selection.

# Slide 27: VarPro cross-validation details
o.cv$imp    # Displays variable importance from cross-validation.
o.cv$imp.conserve    # Displays conservative variable importance.
o.cv$imp.liberal     # Displays liberal variable importance.
o.cv$err             # Displays cross-validation error.
o.cv$zcut            # Displays selected cutoff value.

# Slide 28: VarPro + VIMP CI
o <- rfsrc(Surv(ttodead, died)~., peakVO2, importance="permute")    # Fits survival forest with permutation VIMP.
oo <- subsample(o)    # Computes subsampling-based VIMP confidence intervals.
plot.vimp.ci(oo, alpha=.05)    # Plots VIMP confidence intervals.

o.cv <- cv.varpro(Surv(ttodead, died)~., peakVO2)    # Performs cross-validated VarPro again.
barplot(o.cv$imp.liberal$z, names.arg=o.cv$imp.liberal$variable,
        las=2, horiz = TRUE, col = "coral2")    # Plots liberal VarPro importance values.

# Slide 29: Continuation of Slide 28
# Re-run Slide 28 code to reproduce incremental display.

# Slide 30: VarPro high-dimensional example
data(vdv, package = "randomForestSRC")    # Loads high-dimensional van de Vijver dataset.
dim(vdv)    # Displays dimensions of high-dimensional dataset.

# Slide 31: High-dimensional VarPro with split-weight methods
f <- as.formula(Surv(Time, Censoring)~.)    # Defines survival formula for high-dimensional VarPro.

importance(varpro(f, vdv, split.weight.method = "lasso"))    # Computes importance using lasso split weights.
importance(varpro(f, vdv, split.weight.method = "lasso vimp"))    # Computes importance using lasso + VIMP split weights.
importance(varpro(f, vdv, split.weight.method = "lasso vimp tree"))    # Computes importance using lasso + VIMP + shallow-tree split weights.

rO <- lapply(1:25, function(b) {    # Repeats high-dimensional VarPro comparison 25 times.
  cat("replication:", b, "\n")    # Prints replication number for progress tracking.
  o1 <- varpro(f, vdv, split.weight.method = "lasso")    # Fits VarPro using lasso split weights.
  o2 <- varpro(f, vdv, split.weight.method = "lasso vimp")    # Fits VarPro using lasso + VIMP split weights.
  o3 <- varpro(f, vdv, split.weight.method = "lasso vimp tree")    # Fits VarPro using lasso + VIMP + shallow-tree split weights.
  o4 <- varpro(f, vdv, split.weight.method = "lasso vimp", sparse = FALSE)    # Fits VarPro using dense split weights.
  list("lasso"=intersect(nms,get.orgvimp(o1)$variable),
       "lasso.vimp"=intersect(nms,get.orgvimp(o2)$variable),
       "lasso.vimp.tree"=intersect(nms,get.orgvimp(o3)$variable),
       "lasso.vimp.sparseoff"=intersect(nms,get.orgvimp(o4)$variable))    # Stores overlap with reference signature for each method.
})

# Slide 32: Inspect overlap with reference signature
rO[[1]]    # Displays overlap results for first replication.

# Slide 34: ivarpro — Individual Variable Priority
data(peakVO2, package = "randomForestSRC")    # Reloads peakVO2 dataset for ivarpro demonstration.
o <- varpro(Surv(ttodead, died) ~ ., peakVO2, ntree = 50)    # Fits VarPro model with 50 trees.
ivp <- ivarpro(o)    # Computes individual variable priority matrix.
print(ivp[1:5, 1:8])    # Prints first five rows and eight columns of IVP matrix.

imp <- ivarpro(o, cut.max = 2, adaptive = FALSE)    # Computes IVP with cutoff and no adaptive adjustment.
shap.ivarpro(imp)    # Computes SHAP-style variable priority from IVP.

# Slide 35: plot.ivarpro
plot(imp, var = "peak.vo2", col.var = "interval", size.var = "y")    # Plots IVP results for peak VO2 with color and size encoding.

# Slide 37: Core VarPro functions
# varpro()            # Main supervised variable-priority function.
# partialpro()        # Partial effect computation for VarPro.
# plot()              # Plotting interface for VarPro objects.
# importance()        # Variable importance extraction.
# cv.varpro()         # Cross-validated VarPro.
# uvarpro()           # Unsupervised variable priority.
# sdependent()        # Dependency analysis for unsupervised VarPro.
# get.beta.entropy()  # Entropy-based beta extraction.
# ivarpro()           # Individual variable priority.
# shap.ivarpro()      # SHAP-style variable priority.
# partial.ivarpro()   # Partial effect for IVP.

# Slide 38: uvarpro — Unsupervised Variable Priority
data(BostonHousing, package = "mlbench")    # Loads Boston Housing dataset for unsupervised VarPro example.
uvp <- uvarpro(BostonHousing[-(1:5), ])    # Computes unsupervised variable priority on dataset excluding first five rows.
print(head(importance(uvp)))    # Prints top unsupervised variable importance values.

beta <- get.beta.entropy(uvp)    # Computes entropy-based beta coefficients for unsupervised VarPro.
print(beta[1:4, 1:4])    # Prints first four rows and columns of entropy-based beta matrix.

sdependent(beta)    # Computes dependency structure among variables using entropy-based beta.

# Slide 39: Summary
# varpro()    # Supervised variable priority.
# cv.varpro() # Cross-validated variable priority.
# uvarpro()   # Unsupervised variable priority.
# ivarpro()   # Individual variable priority.
# outpro()    # Out-of-sample variable priority.
# isopro()    # Iso-priority analysis.

###############################################################################
# MODULE 4 — ASA Traveling Course Part IV: Advanced Topics
###############################################################################

library(randomForestSRC)    # Loads randomForestSRC to support advanced forest methods including imbalanced classification and imputation.
library(survival)           # Loads survival to support Surv() objects used in survival and RHF examples.
library(varPro)             # Loads VarPro for supervised and unsupervised variable-priority methods.
library(randomForestSGT)    # Loads Super Greedy Trees for advanced splitting and model explanation.
library(randomForestRHF)    # Loads Random Hazard Forests for counting-process survival modeling.

# Slide 3: Imbalanced classification — syntax examples
# rfsrc(Ozone ~ ., data = airquality)    # Regression forest example.
# quantreg(mpg ~ ., data = mtcars)       # Quantile regression forest example.
# rfsrc(Species ~ ., data = iris)        # Classification forest example.
# imbalanced(status ~ ., data = breast)  # Imbalanced classification example.
# rfsrc(Surv(time, status) ~ ., data = veteran)    # Survival forest example.
# rfsrc(Surv(time, status) ~ ., data = wihs)       # Survival forest example.
# rfsrc(Multivar(mpg, cyl) ~ ., data = mtcars)     # Multivariate forest example.
# rfsrc(cbind(Species, Sepal.Length) ~ ., data = iris)    # Mixed outcome forest example.
# quantreg(cbind(mpg, cyl) ~ ., data = mtcars)     # Quantile regression forest example.
# quantreg(cbind(Species, Sepal.Length) ~ ., data = iris) # Quantile regression forest example.
# rfsrc(data = mtcars)    # Unsupervised forest example.
# sidClustering(data = mtcars)    # SID clustering example.
# sidClustering(data = mtcars, method = "sh")    # SID clustering example with alternative method.

# Slide 11: Imbalanced classification — Glioma
library(varPro)    # Loads VarPro to access glioma dataset and imbalanced classification tools.
data(glioma, package = "varPro")    # Loads glioma dataset for imbalanced classification example.
table(glioma$y)    # Displays class distribution to illustrate imbalance.

class.combine <- c("Classic-like", "Codel", "G-CIMP-high", "Mesenchymal-like")    # Defines majority-class labels to collapse.
ynew <- factor(1 * !is.element(glioma$y, class.combine))    # Creates binary outcome where minority class is labeled 1.

glioma2 <- glioma    # Copies original glioma dataset for modification.
glioma2$y <- ynew    # Replaces multiclass outcome with binary imbalanced outcome.
table(glioma2$y)    # Displays new binary class distribution.

# Slide 12: Standard random forest classifier
o1 <- rfsrc(y ~ ., data = glioma2)    # Fits standard random forest classifier on imbalanced glioma data.
print(o1)    # Prints classifier summary including OOB error.

# Slide 13: RFQ classifier
o2 <- imbalanced(y ~ ., data = glioma2)    # Fits RFQ classifier designed for imbalanced classification.
print(o2)    # Prints RFQ classifier summary.

# Slide 15: Balanced random forest (BRF)
o3 <- imbalanced(y ~ ., data = glioma2, method = "brf")    # Fits balanced random forest using BRF method.
print(o3)    # Prints BRF classifier summary.

# Slide 16: G-mean VIMP with subsampling
o2 <- imbalanced(y ~ ., data = glioma2, importance = "permute", block.size = 20)    # Computes permutation VIMP with block size for imbalanced data.
oo2 <- subsample(o2)    # Performs subsampling to compute VIMP confidence intervals.
plot.subsample(oo2)    # Plots subsampling-based VIMP confidence intervals.

# Slide 18: General call to impute
# impute(formula, data = data, ...)    # Supervised imputation syntax.
# impute(data = data, ...)             # Unsupervised imputation syntax.

# Slide 19: OTFI — On-the-fly imputation
data(pbc, package = "randomForestSRC")    # Loads PBC dataset for imputation examples.
pbc.impute <- impute(Surv(days, status) ~ ., data = pbc)    # Performs supervised on-the-fly imputation for survival data.
pbc.impute <- impute(data = pbc)    # Performs unsupervised on-the-fly imputation.

# Slide 20: missForest and mForest
data(pbc, package = "randomForestSRC")    # Reloads PBC dataset for missForest example.
pbc.impute <- impute(data = pbc, mf.q = 1)    # Performs missForest-style imputation using single regression target.

data(housing, package = "randomForestSRC")    # Loads housing dataset for mForest example.
housing.impute <- impute(data = housing, mf.q = 0.5)    # Performs mForest grouped multivariate imputation.
housing.impute <- impute(data = housing, mf.q = 40)     # Performs mForest with larger group size.

# Slide 21: Test-time imputation using impute.learn
# fit         <- impute.learn(...)    # Trains imputation system.
# newdata.imp <- predict(fit, newdata = ...)    # Applies imputation to new data.
# save.impute.learn(fit, path = ...)    # Saves trained imputation system.
# load.fit    <- load.impute.learn(path = ...)    # Loads saved imputation system.

# Slide 23: Test-time imputation example
aq <- airquality[, c("Ozone", "Solar.R", "Wind", "Temp", "Month")]    # Selects subset of airquality dataset for imputation example.
id <- sample(seq_len(nrow(aq)), 100)    # Randomly selects 100 rows for training.
train <- aq[id, ]    # Creates training subset.
test <- aq[-id, ]    # Creates test subset.

fit <- impute.learn(
  data = train,
  mf.q = 1,
  max.iter = 5,
  full.sweep.options = list(ntree = 25, nsplit = 5),
  target.mode = "all"
)    # Trains imputation system using missForest-style grouped regressions.

test.imp <- predict(fit, test, max.predict.iter = 2)    # Applies trained imputation system to test data.

# Slide 26: OOD scoring example
aq <- airquality[, c("Ozone", "Solar.R", "Wind", "Temp", "Month")]    # Reloads airquality subset for OOD scoring example.
id <- sample(seq_len(nrow(aq)), 100)    # Randomly selects training rows.
train <- aq[id, ]    # Creates training subset.
test <- aq[-id, ]    # Creates test subset.

sup.fit <- impute.learn(
  data = train,
  mf.q = 1,
  supervised.formula = Solar.R ~ .,
  supervised.args = list(ntree = 50, nsplit = 5),
  full.sweep.options = list(ntree = 25, nsplit = 5),
  save.ood = TRUE
)    # Trains supervised imputation system with OOD scoring enabled.

ood <- impute.ood(sup.fit, test)    # Computes OOD scores for test data.
print(head(ood$score))    # Displays first few OOD scores.
print(head(ood$score.percentile))    # Displays percentile-based OOD scores.

# Slide 30: Super Greedy Trees (SGT)
library(randomForestSGT)    # Loads SGT package for advanced splitting and model explainers.

# rfsgt(formula, data, ...)    # Canonical SGT interface.
# hcut = 0 gives CART splits; larger hcut values enlarge geometric dictionary.

# Slide 32: Tuning hcut
n <- 2500    # Sets number of observations for Friedman-1 simulation.
p <- 50      # Sets number of noise variables.
noise <- matrix(runif(n * p), ncol = p)    # Generates uniform noise matrix.
dta <- data.frame(mlbench:::mlbench.friedman1(n, sd = 0), noise = noise)    # Creates dataset combining Friedman-1 signal and noise.

filter <- tune.hcut(y ~ ., data = dta, hcut = 3)    # Tunes hcut parameter up to value 3 for SGT splitting.

# Slide 33: Using tuned hcut
o.sgt <- rfsgt(y ~ ., data = dta, filter = filter)    # Fits SGT forest using tuned hcut and selected basis functions.
print(o.sgt)    # Prints SGT forest summary.

# Slide 34: Specific hcut families
o.hcut0 <- rfsgt(y ~ ., data = dta, filter = use.tune.hcut(filter, hcut = 0))    # Fits SGT forest using CART-style axis-aligned splits.
print(o.hcut0)    # Prints SGT forest summary for hcut=0.

o.hcut1 <- rfsgt(y ~ ., data = dta, filter = use.tune.hcut(filter, hcut = 1))    # Fits SGT forest using hyperplane splits.
print(o.hcut1)    # Prints SGT forest summary for hcut=1.

# Slide 35: SGTs as model explainers
o <- rfsgt(y ~ ., data = dta, pure.lasso = TRUE, filter = filter, treesize = 5)    # Fits shallow SGT forest with lasso-based local models.
bo <- get.beta(o, bag = "oob")    # Extracts OOB predictions, local beta coefficients, and partial contributions.
yhat <- bo$predicted    # Stores OOB predictions from SGT model.
beta <- bo$beta[, -1, drop = FALSE]    # Extracts local beta coefficients excluding intercept.
partial <- bo$partial[, -1, drop = FALSE]    # Extracts partial-effect contributions excluding intercept.

# Slide 36: SGT beta and partial effects
print(head(yhat))    # Displays first few OOB predictions from SGT model.
print(head(o$predicted.oob))    # Displays first few OOB predictions from forest ensemble.
print(head(beta[, 1:6]), digits = 2)    # Displays first few local beta coefficients.
print(head(partial[, 1:6]), digits = 2)    # Displays first few partial-effect contributions.

# Slide 37: Random Hazard Forests (RHF)
library(randomForestRHF)    # Loads RHF package for counting-process survival forests.

# rhf(Surv(id, start, stop, event) ~ ., data = data)    # Canonical RHF syntax.

# Slide 38: RHF data format
# id     subject identifier
# start  interval start time
# stop   interval end time
# event  event indicator at stop

# Slide 39: Time-static RHF example
data(peakVO2, package = "randomForestSRC")    # Loads peakVO2 dataset for RHF example.
d <- convert.counting(Surv(ttodead, died) ~ ., data = peakVO2)    # Converts one-row-per-subject data to counting-process format.

f <- "Surv(id, start, stop, event) ~ ."    # Defines RHF formula using counting-process response.
o <- rhf(f, data = d)    # Fits Random Hazard Forest using counting-process data.

# Slide 40: Time-static RHF results
print(o)    # Prints RHF summary including hazard estimates.

# Slide 41: Time-dependent AUC
fit.n1 <- rhf(f, data = d, nodesize = 1)    # Fits RHF with small terminal-node size.
fit.n15 <- rhf(f, data = d, nodesize = 15)  # Fits RHF with larger terminal-node size.

auc.n1.chf <- auct.rhf(fit.n1)    # Computes AUC-t using cumulative hazard marker for nodesize=1.
auc.n1.haz <- auct.rhf(fit.n1, marker = "haz")    # Computes AUC-t using hazard marker for nodesize=1.
auc.n15.chf <- auct.rhf(fit.n15)    # Computes AUC-t using cumulative hazard marker for nodesize=15.
auc.n15.haz <- auct.rhf(fit.n15, marker = "haz")    # Computes AUC-t using hazard marker for nodesize=15.

ylim <- c(0.6, 0.85)    # Sets y-axis limits for AUC-t plots.
par(mfrow = c(2, 2))    # Sets plotting layout to 2x2 grid.

plot(auc.n1.chf, ylim = ylim, main = "nodesize 1, CHF marker")    # Plots AUC-t curve for nodesize=1 using CHF marker.
plot(auc.n1.haz, ylim = ylim, main = "nodesize 1, hazard marker")    # Plots AUC-t curve for nodesize=1 using hazard marker.
plot(auc.n15.chf, ylim = ylim, main = "nodesize 15, CHF marker")    # Plots AUC-t curve for nodesize=15 using CHF marker.
plot(auc.n15.haz, ylim = ylim, main = "nodesize 15, hazard marker")    # Plots AUC-t curve for nodesize=15 using hazard marker.

# Slide 43: Hazard plots
shaz.n15 <- smoothed.hazard(fit.n15)    # Computes smoothed hazard curves for nodesize=15 RHF model.
id <- fit.n15$ensemble.id[1:3]    # Selects first three subject IDs for hazard plotting.

par(mfrow = c(1, 2))    # Sets plotting layout to two side-by-side panels.
plot(fit.n15, idx = id, main = "OOB Hazard")    # Plots OOB hazard curves for selected subjects.
plot(shaz.n15, idx = id, main = "Smoothed Hazard")    # Plots smoothed hazard curves for selected subjects.

# Slide 45: Time-localized VarPro importance
imp.t <- importance.rhf(fit.n15)    # Computes time-localized variable importance across full time grid.
plot(imp.t, type = "dotmatrix")    # Plots time-localized importance using dot-matrix display.
plot(imp.t, type = "lines")        # Plots time-localized importance using line display.

###############################################################################
# MODULE 5 — Unified Orchestrator for ASA Parts I–IV
###############################################################################

# This orchestrator assumes Modules 0–4 have already been sourced or run.

run_master_asa <- function() {    # Defines a master function that runs selected ASA workflows and returns structured results.
  
  message("Starting ASA Master Workflow...")    # Prints a message indicating the workflow has begun.
  
  # --- Dataset Preparation (from Module 0) ---
  gbcs_data <- gbcs_asa    # Stores the ASA‑ready GBCS dataset for use in all modules.
  
  # --- Part I: Training ---
  asa_part1 <- list()    # Creates a list to store Part I results.
  asa_part1$housing_rf <- rfsrc(SalePrice ~ ., data = housing)    # Trains a random forest on the housing dataset as in ASA Part I.
  asa_part1$gbcs_rf_reg <- rfsrc(y_reg ~ ., data = gbcs_data)    # Trains a regression forest on GBCS using survival time as target.
  asa_part1$gbcs_rf_class <- rfsrc(y_class ~ ., data = gbcs_data)    # Trains a classification forest on GBCS using death indicator.
  
  # --- Part II: Inference & Prediction ---
  asa_part2 <- list()    # Creates a list to store Part II results.
  asa_part2$gbcs_surv <- rfsrc(Surv(time, status) ~ ., data = gbcs_data)    # Fits a survival forest on GBCS for inference tasks.
  asa_part2$gbcs_pred <- predict(asa_part2$gbcs_surv, gbcs_data[1:10, ])    # Predicts survival outcomes for first 10 GBCS cases.
  asa_part2$gbcs_partial_age <- plot.variable(asa_part2$gbcs_surv, xvar.names = "age", partial = TRUE)    # Computes partial effect of age.
  
  # --- Part III: Variable Selection ---
  asa_part3 <- list()    # Creates a list to store Part III results.
  asa_part3$vimp <- vimp(rfsrc(y_reg ~ ., data = gbcs_data), importance = "permute")    # Computes permutation VIMP for GBCS regression forest.
  asa_part3$md <- max.subtree(rfsrc(Surv(time, status) ~ ., data = gbcs_data))$order[,1]    # Computes minimal depth for GBCS survival forest.
  asa_part3$varpro <- varpro(Surv(time, status) ~ ., gbcs_data)    # Runs supervised VarPro on GBCS dataset.
  
  # --- Part IV: Advanced Topics ---
  asa_part4 <- list()    # Creates a list to store Part IV results.
  asa_part4$imbalanced <- imbalanced(y_class ~ ., data = gbcs_data)    # Fits imbalanced classifier on GBCS.
  asa_part4$impute_unsup <- impute(data = gbcs_data)    # Performs unsupervised imputation on GBCS.
  asa_part4$sgt <- rfsgt(y_reg ~ ., data = gbcs_data, filter = tune.hcut(y_reg ~ ., data = gbcs_data, hcut = 2))    # Fits SGT forest on GBCS.
  asa_part4$gbcs_count <- convert.counting(Surv(time, status) ~ ., data = gbcs_data)    # Converts GBCS to counting-process format for RHF.
  asa_part4$rhf <- rhf("Surv(id, start, stop, event) ~ .", data = asa_part4$gbcs_count)    # Fits Random Hazard Forest on counting-process GBCS.
  
  # --- Return Structured Output ---
  return(list(
    dataset = gbcs_data,    # Returns ASA-ready GBCS dataset.
    part1 = asa_part1,      # Returns results from ASA Part I.
    part2 = asa_part2,      # Returns results from ASA Part II.
    part3 = asa_part3,      # Returns results from ASA Part III.
    part4 = asa_part4       # Returns results from ASA Part IV.
  ))    # Returns a unified list containing all module outputs.
}

# Example call (commented out for safety):
# results_master <- run_master_asa()    # Runs the full ASA master workflow and stores results.
