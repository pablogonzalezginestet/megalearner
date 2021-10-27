library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/")

## get the breast cancer data:  Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)
mycancerstudy <- getCancerStudies(mycgds)[37, 1]
## get all patient ids
mycaselist <- getCaseLists(mycgds,mycancerstudy)[2, 1]

data.clin <- getClinicalData(mycgds, mycaselist)
head(data.clin)
dim(data.clin)
data.clin = data.frame(data.clin)
colnames(data.clin)

metabric = data.clin[c('RFS_MONTHS','RFS_STATUS','AGE_AT_DIAGNOSIS','TUMOR_SIZE',"GRADE","LYMPH_NODES_EXAMINED_POSITIVE", "INFERRED_MENOPAUSAL_STATE","HORMONE_THERAPY","ER_STATUS","PR_STATUS",'NPI',"RADIO_THERAPY")]

head(metabric)
metabric$status = factor(metabric$RFS_STATUS ,c('0:Not Recurred','1:Recurred'), labels = c(0,1))
metabric$meno = factor(metabric$INFERRED_MENOPAUSAL_STATE,c('Post','Pre'), labels = c(1,0))
metabric$hormon = factor(metabric$HORMONE_THERAPY,c('YES','NO'), labels = c(1,0))
metabric$radio = factor(metabric$RADIO_THERAPY,c('YES','NO'), labels = c(1,0))
metabric$er = factor(metabric$ER_STATUS,c('Positive','Negative'), labels = c(1,0))
metabric$pr = factor(metabric$PR_STATUS,c('Positive','Negative'), labels = c(1,0))

metabric$node = ifelse(metabric$LYMPH_NODES_EXAMINED_POSITIVE==0,0,1)
metabric$grade2 = 1*(metabric$GRADE==2)
metabric$grade3 = 1*(metabric$GRADE==3)

library(survival)

devel <- data.frame(YY = with(metabric, Surv(round(RFS_MONTHS), unclass(status))),
                    age = metabric$AGE_AT_DIAGNOSIS, meno = unclass(metabric$meno),
                    nodes = unclass(metabric$LYMPH_NODES_EXAMINED_POSITIVE), 
                    hormon = unclass(metabric$hormon),
                    radio = unclass(metabric$radio),
                    pr= unclass(metabric$pr), er=unclass(metabric$er),
                    grade2 = (metabric$grade2),
                    grade3 = (metabric$grade3),
                    npi = metabric$NPI,
                    size = metabric$TUMOR_SIZE)


devel = devel[!is.na(devel$grade2),]
devel = devel[!is.na(devel$size),]
devel = devel[!is.na(devel$YY),]
summary(devel)
dim(devel)
head(devel)

#loading package
set.seed(420)
data1 = sort(sample(nrow(devel), nrow(devel)*.7))
#creating training data set by selecting the output row values
train<-devel[data1,]
#creating test data set by not selecting the output row values
test<-devel[-data1,]

#################################################################################

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
############### Natively ############################

stepwise <- function(dset) {
  
  cfit <- coxph(YY ~ ., data = dset)
  cfin <- step(cfit, trace = 0)
  function(valid) {
    predict(cfin, newdata = valid,
            type = "lp")
  }
  
}

random.forest <- function(dset, time = 5 * 12) {
  
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

direct.binomial <- function(dset, time = 5 * 12) {
  dset$time <- dset$YY[, "time"]
  dset$status <- dset$YY[, "status"]
  dset$YY <-  NULL
  dset$ybin <- 1.0 * (dset$time < time & dset$status == 1)
  dset$ybin[dset$time < time & dset$status == 0] <- NA
  
  sfit <- survfit(Surv(time, 1 - status) ~ 1, data = dset)
  dset$weights <- 1 / summary(sfit, times = pmin(dset$time, time))$surv
  dset$time <- dset$status <- NULL
  
  cfit <- glm(ybin ~ bs(age) + meno + nodes + hormon + radio + pr + er + grade2 + grade3 + bs(npi) +  bs(size),
              data = dset, family = "binomial", weights = weights)
  cfin <- step(cfit, trace = 0)
  function(valid) {
    predict(cfit, newdata = valid, type = "response")
  }
  
}

## bagging ipcw

svm <- function(dset, time = 5 * 12){
  
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



knn <- function(dset, time = 5 * 12){
  
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

pseudo.glm <- function(dset, time = 5 * 12){
  
  cfit <- cumincglm(YY ~ bs(age) + meno + nodes + hormon + radio + pr + er + grade2 + grade3 + bs(npi) +  bs(size),
                    data = dset, time = time)
  function(valid) {
    predict(cfit, newdata = valid, type = "response")
  }
  
}


pseudo.glmnet <- function(dset, time = 5 * 12) {
  
  cfit <- cumincglm(YY ~ bs(age) + meno + nodes + hormon + radio + pr + er + grade2 + grade3 + bs(npi) +  bs(size), data = dset, time = time, x = TRUE)
  
  cnfit <- cv.glmnet(cfit$x, cfit$y)
  
  function(valid) {
    vx <- model.matrix(~ bs(age) + meno + nodes + hormon + radio + pr + er + grade2 + grade3 + bs(npi) +  bs(size), data = valid)
    predict(cnfit, newx = vx)[, 1]
    
  }
  
}


#16 models
mymods <- list(stepwise, random.forest, coxboost,
               direct.binomial, svm, knn,
               function(xx) pseudo.glm(xx, time = 3 * 12),
               function(xx) pseudo.glm(xx, time = 4 * 12),
               function(xx) pseudo.glm(xx, time = 5 * 12),
               function(xx) pseudo.glm(xx, time = 6 * 12),
               function(xx) pseudo.glm(xx, time = 7 * 12),
               function(xx) pseudo.glmnet(xx, time = 3 * 12),
               function(xx) pseudo.glmnet(xx, time = 4 * 12),
               function(xx) pseudo.glmnet(xx, time = 5 * 12),
               function(xx) pseudo.glmnet(xx, time = 6 * 12),
               function(xx) pseudo.glmnet(xx, time = 7 * 12))

set.seed(420)
ndex <- 1:nrow(train)
part <- as.factor(sample(1:10, length(ndex), replace = TRUE))
folds <- split(ndex, part)

Zout <- vector(mode = "list", length = 10)
for(j in 1:length(folds)) {
  
  training <- train[unlist(folds[-j]), ]
  validation <- train[folds[[j]], ]
  
  Zout[[j]] <- do.call(cbind, lapply(mymods, function(f){
    
    fhat <- f(training)
    fhat(validation)
    
  }))
  
}


fullfits <- lapply(mymods, function(f) {
  
  f(train)
  
})

# fullfits is list where each contain the models
saveRDS(fullfits, "/fullfits_rfs.rds")
#Zout has for each fold the prediction of each model so it is a list where the columns of each list is the model
saveRDS(Zout, "/Zout_rfs.rds") 

saveRDS(folds, "/split_rfs.rds")


##########################  part 2 ############################

library(pseudoloss) ## remotes::install_github("sachsmc/pseudoloss")


Zout <- readRDS("/Zout_rfs.rds")

fullfits <- readRDS("/fullfits_rfs.rds")
folds <- readRDS("/split_rfs.rds")

train$PO <- cumincglm(YY ~ 1, data = train, time = 5 * 12)$y

YYens <- train$PO[unlist(folds)]
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

saveRDS(cv.auc, "/cv-aucs-stack_rfs.rds")

k <- kcand[which.max(colMeans(cv.auc))]
fn.auc <- function(beta) {
  
  yy <- c(Zmat %*% matrix(beta, ncol = 1))
  1 - with(calc_roc(yy, matrix(YYens, ncol = 1), ramp = smoothramp),
           calc_auc(fpf, tpf)) + k * sum(abs(beta)^2)
}

start.beta <- runif(16)
#the optimal alpha
#opt.fit <- optim(par = start.beta, fn = fn.auc, method = "BFGS")
opt.fit <- optim(par = start.beta, fn = fn.auc, method = "L-BFGS-B",lower = 0.001)
saveRDS(opt.fit, "/opt-coeffs_rfs_with_restriction.rds")

###################################### part 3 ######################################
library(pseudoloss)  ##

library(xtable)
library(ggplot2)
library(patchwork)
library(broom)

fullfits <- readRDS("/fullfits_rfs.rds")
cv.auc <- readRDS("/cv-aucs-stack_rfs.rds")
opt.fit <- readRDS("/opt-coeffs_rfs_with_restriction.rds")


valid = test

Zvalid <- lapply(fullfits, function(f) f(valid))
YYensvalid <- cumincglm(YY ~ 1, data = valid, time = 5 * 12)$y

aucalls <- lapply(Zvalid, function(z) {
  
  with(calc_roc(z, matrix(YYensvalid, ncol = 1), ramp = smoothramp),
       calc_auc(fpf, tpf))
  
})

aucalls1 <- lapply(Zvalid, function(z) {
  
  calc_roc(z, matrix(YYensvalid, ncol = 1), ramp = smoothramp)
  
})

names.models <- c("Cox.stepwise", "SRF", "CoxBoost",
                  "IPCW logistic stepwise", "Bagged SVM", "Bagged KNN",
                  "eventglm t1", "eventglm t2", "eventglm t3", "eventglm t4", "eventglm t5",
                  "LASSO eventglm t1", "LASSO eventglm t2", "LASSO eventglm t3",
                  "LASSO eventglm t4", "LASSO eventglm t5")


Zmat.valid <- do.call(cbind, Zvalid)

Zfin.opt <- (Zmat.valid %*% (opt.fit$par) / sum(opt.fit$par))[, 1]
roc.opt <- calc_roc(Zfin.opt, matrix(YYensvalid, ncol = 1), ramp = smoothramp)

round(with(roc.opt, pseudoloss::calc_auc(fpf,tpf)),3)

sumtable <- data.frame(R.package = c("survival", "randomForestSRC", "CoxBoost", "stats",
                                     "e1071", "class", rep("eventglm", 5),
                                     rep("glmnet", 5)),
                       AUC = unlist(aucalls), coefficient = (opt.fit$par) / sum(opt.fit$par))
rownames(sumtable) <- names.models


print(xtable(sumtable[, 1:3], digits = 3))



## roc and predictiveness curve

roc.opt$cut <- sort(unique(Zfin.opt))

index_model = order(-opt.fit$par)[c(1,3,6)]
temp = roc.opt
temp$method = 'Optimized Stack'
for(i in index_model){
  temp1 = aucalls1[[i]]
  temp1$cut <- sort(unique(Zmat.valid[,i]))
  temp1$method = names.models[i]
  temp = rbind(temp,temp1)
}


points_auc = temp[temp$cut>0.3,] %>% 
  group_by(method) %>% 
  slice(which.min(cut))


names.models[c(index_model)]

p1<- ggplot(temp, aes(x = fpf, y = tpf,color = method)) + geom_line() +
  scale_color_manual(
    values = c("#999999", "#E69F00","#F0E442", "#56B4E9"
    ),
    labels = c(  "Cox.stepwise","eventglm t1" ,'Optimized Stack',"SRF"))+
  geom_point(data = points_auc,
             aes(x = fpf, y = tpf, group = method)) +
  theme_bw() +plotROC::style_roc() + geom_abline(intercept = 0, slope = 1, color = "grey90")



Zopt.cut <- pmin(pmax(0, Zfin.opt), .99)


p2 <- ggplot(data.frame(pred.risk = Zopt.cut, pos = YYensvalid), aes(x = pred.risk, y = pos)) +
  stat_smooth(color = "black", method = "loess", span = 1.25, se = FALSE) +
  scale_x_continuous("Predicted 5 year risk") +
  scale_y_continuous("Estimated 5 year risk") +
  geom_rug(aes(x= pred.risk, y = NULL), sides = "b") +
  geom_hline(yintercept = mean(YYensvalid), linetype = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()
p1 + p2
ggsave("/int-plot_metabric_rfs_with_restriction.pdf", width = 7.25, height = 3.75)


valid$Zopt.cut <- Zopt.cut
survConcordance(YY ~ Zopt.cut, data = valid)


valid$group.Z <- cut(valid$Zopt.cut,
                     breaks = c(-1, quantile(valid$Zopt.cut, c(0.25, .5, .75)), 2),
                     labels = c("Q1", "Q2", "Q3", "Q4"))


sfitg <- survfit(YY ~ group.Z, data = valid)
meplot <- tidy(sfitg)
meplot$strata <- gsub("group.Z=", "", meplot$strata, fixed = TRUE)

ggplot(meplot, aes(x = time, y = estimate, color = strata)) +
  geom_step() +
  geom_step(aes(y = conf.low), linetype = 3) +
  geom_step(aes(y = conf.high), linetype = 3) +
  theme_classic() +
  scale_color_grey("Risk quartile", start = .7, end = .1) +
  xlab("Months since diagnosis") + ylab("Survival probability")
ggsave("/km-fig_metabric_rfs_with_restriction.pdf", width = 5.25, height = 3.75)

rpm <- replicate(200, {
  coxfit <- coxph(YY ~ I(Zopt.cut * 10), data = valid[sample(1:nrow(valid), nrow(valid), replace = TRUE), ])
  royston(coxfit)[3]
})
#explained variation
mean(rpm)
sd(rpm)
#######
summary(coxph(YY ~ I(Zopt.cut * 10), data = valid))
summary(coxph(YY ~ group.Z, data = valid))


aucb <- replicate(200, {
  
  smi <- sample(1:length(Zfin.opt), length(Zfin.opt), replace = TRUE)
  roc.opt <- calc_roc(Zfin.opt[smi], matrix(YYensvalid[smi], ncol = 1), ramp = smoothramp)
  with(roc.opt, pseudoloss::calc_auc(fpf,tpf))
  
})
sd(aucb)


summary(sfitg, times= (0:3)*100)
