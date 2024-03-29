
<!-- README.md is generated from README.Rmd. Please edit that file -->

## megalearner

<!-- badges: start -->

<!-- badges: end -->

**megalearner** is series of R scripts, found in the folder R, that
implement the illustrative example described in the paper Gonzalez
Ginestet, P. et al. (2021). “Survival stacking with multiple data types
using pseudo-observation-based-AUC loss”.

**megalearner** builds a more accurate risk prediction model than
standard methods such as Cox proportional hazard or single machine
learning algorithms that have been already adapted to handle
right-censored observations. **megalearner** optimally stacks
predictions from any method based on different survival outcome types
(either binary or continuous). In the paper, we show that our proposed
method can improve on single survival based methods such as survival
random forest or on other strategies that use a pre-processing step such
as inverse probability of censoring weighted or pseudo-observations.

We illustrate our methodology in two datasets:

1) We use the  Rotterdam data set as training data set and we validate the
predictions in an external data set. For that, we use the data from the
German Breast Cancer Study Group (GBSG). Both data sets are available in
the `survival` R package \[Therneau, 2000\].

2) We use the METABRIC cohort for two endpoints: overall survival (OS)
and recurrence-free survival (RFS) in months.

## Organization of the folders

The scripts depend on each other in the following way:

03-validate.R -\> depends on -\> 02-opt-stack.R -\> depends on -\>
01-train-models.R -\> depends on 00-model-specs.R

01-train-models.R: saves the fit of each algorithm (“fullfits.rds”), the
matrix Z (“Zout.rds”) where the predictions of each algorithms are
recorded and the folds used in the cross-validation (“split.rds”)

02-opt-stack.R: saves the optimal coefficient used to combine the
algorithms (“opt-coeffs.rds”) and the AUC at 5 years cross-validated of
the stacking (“cv-aucs-stack.rds”)

The folder METABRIC contains the R codes to run the procedure in the METABRIC dataset for the two predictive examples.

The folder figure contains the plots in the paper. 

## Installation

Depends on the following R packages: `broom`, `patchwork`, `ggplot2`,
`xtable`, `class`, `glmnet`, `e1071`, `riskRegression`, `CoxBoost`,
`prodlim`, `Matrix`, `randomForestSRC`, `survival`, `eventglm`,
`pseudoloss`, `splines`, `stats`, `graphics`, `grDevices`, `utils`,
`datasets`, `methods`, `base`, `plotROC`

The following packages `CoxBoost`, `eventglm` and `pseudoloss` can be
installed as follows:

``` r
remotes::install_github("binderh/CoxBoost")
remotes::install_github("sachsmc/eventglm")
remotes::install_github("sachsmc/pseudoloss")
```
