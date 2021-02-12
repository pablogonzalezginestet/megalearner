library(eventglm) ## remotes::install_github("sachsmc/eventglm")
library(survival)
library(randomForestSRC)
library(CoxBoost)  ## remotes::install_github("binderh/CoxBoost")
library(riskRegression)
library(e1071)
library(splines)
library(glmnet)
library(class)

## models

## survival

stepwise <- function(dset) {

  cfit <- coxph(YY ~ ., data = dset)
  cfin <- step(cfit, trace = 0)
  function(valid) {
    predict(cfin, newdata = valid,
            type = "lp")
  }

}

random.forest <- function(dset, time = 5 * 365.25) {

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL

  cfit <- rfsrc(Surv(time, status) ~ ., data = dset)
  function(valid) {
    res <- predict(cfit, newdata = valid[, -1])
    1 - res$survival[, max(which(res$time.interest <= time))[1]]
  }

}

coxboost <- function(dset) {

  mm <- model.matrix(YY ~ ., data = dset)[, -1]
  cfit <- CoxBoost(time = dset$YY[, "time"], status = dset$YY[, "status"],
           x = mm)
  function(valid) {
    mmv <- model.matrix(YY ~ ., data = valid)[, -1]
    res <- predict(cfit, newdata = mmv,
                   type = "lp")
    res[1, ]
  }

}

## binary

direct.binomial <- function(dset, time = 5 * 365.25) {
  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, 1 - status) ~ 1, data = dset)
  dset$weights <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$time <- dset$status <- NULL

  cfit <- glm(ybin ~ bs(age) + meno + sizemed + sizebig + bs(nodes) + bs(pgr) + bs(er) +
                hormon,
      data = dset, family = "binomial", weights = weights)
  cfin <- step(cfit, trace = 0)
  function(valid) {
      predict(cfit, newdata = valid, type = "response")
  }

}

## bagging ipcw

svm <- function(dset, time = 5 * 365.25){

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, 1 - status) ~ 1, data = dset)

  dset <- dset[!is.na(dset$ybin),]
  wtmp <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$samp.wts <- wtmp / sum(wtmp)
  dset$time <- dset$status <- NULL

  svmboot <- lapply(1:10, function(i) {

    dboot <- dset[sample(1:nrow(dset), nrow(dset),
                         replace = TRUE, prob = dset$samp.wts),]
    dboot$samp.wts <- NULL
    e1071::svm(ybin ~ ., data = dboot)

  })


  function(valid) {

    vboots <- sapply(svmboot, function(cfit) {
      predict(cfit, newdata = valid)
    })
    rowMeans(vboots)

  }
}



knn <- function(dset, time = 5 * 365.25){

  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA

  sfit <- survfit(Surv(time, 1 - status) ~ 1, data = dset)

  dset <- dset[!is.na(dset$ybin),]
  wtmp <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$samp.wts <- wtmp / sum(wtmp)
  dset$time <- dset$status <- NULL


  function(valid) {

    vboots <- sapply(1:10, function(i) {

      dboot <- dset[sample(1:nrow(dset), nrow(dset),
                           replace = TRUE, prob = dset$samp.wts),]
      dboot$samp.wts <- NULL
      tclass <- as.factor(dboot$ybin)
      dboot$ybin <- NULL

      knboot <- class::knn(dboot, valid[, colnames(dboot)],
          cl = tclass, k = 5, prob = TRUE)

      ifelse(knboot == 0, 1 - attr(knboot, "prob"), attr(knboot, "prob"))

    })

    rowMeans(vboots)

  }
}


## pseudo obs

pseudo.glm <- function(dset, time = 5 * 365.25){

  cfit <- cumincglm(YY ~ bs(age) + meno + sizemed + sizebig + bs(nodes) + bs(pgr) + bs(er) +
                      hormon, data = dset, time = time)
  function(valid) {
    predict(cfit, newdata = valid, type = "response")
  }

}


pseudo.glmnet <- function(dset, time = 5 * 365.25) {

  cfit <- cumincglm(YY ~ bs(age) + meno + sizemed + sizebig + bs(nodes) + bs(pgr) + bs(er) +
                      hormon, data = dset, time = time, x = TRUE)

  cnfit <- cv.glmnet(cfit$x, cfit$y)

  function(valid) {
    vx <- model.matrix(~ bs(age) + meno + sizemed + sizebig + bs(nodes) + bs(pgr) + bs(er) +
                         hormon, data = valid)
    predict(cnfit, newx = vx)[, 1]

  }

}


