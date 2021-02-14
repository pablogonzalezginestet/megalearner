
<!-- README.md is generated from README.Rmd. Please edit that file -->

## megalearner

<!-- badges: start -->

<!-- badges: end -->

**megalearner** is an `R` package that implements the method in the
paper Gonzalez Ginestet, P. et al. (2021). “Survival stacking with
multiple data types using pseudo-observation-based-AUC loss”.

This series of `R` scripts found in the folder R runs the illustrative
example as described in the paper. The scripts depend on each other in
the following way:

03-validate.R -\> depends on -\> 02-opt-stack.R -\> depends on -\>
01-train-models.R -\> depends on 00-model-specs.R

01-train-models.R: saves the fit of each algorithm (“fullfits.rds”), the
matrix Z (“Zout.rds”) and the folds in the cross-validation
(“split.rds”)

02-opt-stack.R: saves the optimal coefficient used to combine the
algorithms (“opt-coeffs.rds”) and the AUC cross-validates of the
stacking (“cv-aucs-stack.rds”)

## Installation

You can install the development version of `megalearner` from GitHub
with:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("pablogonzalezginestet/megalearner")
```

Among other packages, `megalearner` requires to have install these
packages that can be installed as follow:

``` r
remotes::install_github("binderh/CoxBoost")
remotes::install_github("sachsmc/eventglm")
remotes::install_github("sachsmc/pseudoloss")
```
