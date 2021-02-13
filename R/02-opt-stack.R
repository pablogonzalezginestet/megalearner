library(pseudoloss) ## remotes::install_github("sachsmc/pseudoloss")
#source("00-model-specs.R")

Zout <- readRDS("Zout.rds")
fullfits <- readRDS("fullfits.rds")
folds <- readRDS("split.rds")

data("gbsg")
data("rotterdam")

rotterdam <- subset(rotterdam, nodes > 0)

devel <- data.frame(YY = with(rotterdam, Surv(pmin(rtime, dtime), 1.0 * (recur == 1 | unclass(death) == 1))),
                    age = rotterdam$age, meno = unclass(rotterdam$meno),
                    sizemed = 1.0 * (rotterdam$size == "20-50"),
                    sizebig = 1.0 * (rotterdam$size == ">50"),
                    nodes = unclass(rotterdam$nodes), pgr = unclass(rotterdam$pgr),
                    er = unclass(rotterdam$er), hormon = unclass(rotterdam$hormon))

devel$PO <- cumincglm(YY ~ 1, data = devel, time = 5 * 365.25)$y

YYens <- devel$PO[unlist(folds)]
Zmat <- do.call(rbind, Zout)

lassofit <- cv.glmnet(Zmat, YYens, standardize = FALSE, alpha = 0)
lfit <- glmnet(Zmat, YYens, standardize = FALSE,
               lambda = lassofit$lambda.1se, alpha = 0)

kcand <- exp(seq(-5, 0))

set.seed(1890)
ndex <- 1:nrow(Zmat)
part <- as.factor(sample(1:10, length(ndex), replace = TRUE))
folds <- split(ndex, part)

cv.auc <- matrix(NA, nrow = 10, ncol = length(kcand))

for(k in 1:length(kcand)) {
  for(i in 1:length(folds)){

    Zin <- Zmat[unlist(folds[-i]),]
    Zout <- Zmat[folds[[i]],]
    YYin <- YYens[unlist(folds[-i])]
    YYout <- YYens[folds[[i]]]

    fn.auc <- function(beta) {

      yy <- c(Zin %*% matrix(beta, ncol = 1))
      1 - with(calc_roc(yy, matrix(YYin, ncol = 1), ramp = smoothramp),
               calc_auc(fpf, tpf)) + k * sum(abs(beta)^2)
    }

    start.beta <- as.matrix(lfit$beta)[, 1]
    opt.fit <- optim(par = start.beta, fn = fn.auc, method = "BFGS")

    cv.auc[i, k] <- with(calc_roc((Zout %*% opt.fit$par)[, 1],
                                  matrix(YYout, ncol = 1), smoothramp),
                      calc_auc(fpf, tpf))

  }
}

saveRDS(cv.auc, "cv-aucs-stack.rds")

k <- kcand[which.max(colMeans(cv.auc))]
fn.auc <- function(beta) {

  yy <- c(Zmat %*% matrix(beta, ncol = 1))
  1 - with(calc_roc(yy, matrix(YYens, ncol = 1), ramp = smoothramp),
           calc_auc(fpf, tpf)) + k * sum(abs(beta)^2)
}

start.beta <- opt.fit$par

opt.fit <- optim(par = start.beta, fn = fn.auc, method = "BFGS")

saveRDS(opt.fit, "opt-coeffs.rds")
