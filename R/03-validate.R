library(pseudoloss)  ##
#source("00-model-specs.R")

library(xtable)
library(ggplot2)
library(patchwork)
library(broom)

fullfits <- readRDS("fullfits.rds")
cv.auc <- readRDS("cv-aucs-stack.rds")
opt.fit <- readRDS("opt-coeffs.rds")

data("gbsg")

valid <- data.frame(YY = with(gbsg, Surv(rfstime, status)),
                    age = gbsg$age, meno = unclass(gbsg$meno),
                    sizemed = 1.0 * (gbsg$size > 20 & gbsg$size < 50),
                    sizebig = 1.0 * (gbsg$size >50),
                    nodes = unclass(gbsg$nodes), pgr = unclass(gbsg$pgr),
                    er = unclass(gbsg$er), hormon = unclass(gbsg$hormon))
Zvalid <- lapply(fullfits, function(f) f(valid))
YYensvalid <- cumincglm(YY ~ 1, data = valid, time = 5 * 365.25)$y

aucalls <- lapply(Zvalid, function(z) {

  with(calc_roc(z, matrix(YYensvalid, ncol = 1), ramp = smoothramp),
       calc_auc(fpf, tpf))

})

names.models <- c("Cox.stepwise", "survival random forests", "CoxBoost",
                  "IPCW logistic stepwise", "Bagged SVM", "Bagged KNN",
                  "eventglm 2yr", "eventglm 4yr", "eventglm 5yr", "eventglm 6yr", "eventglm 7yr",
                  "LASSO eventglm 2yr", "LASSO eventglm 4yr", "LASSO eventglm 5yr",
                  "LASSO eventglm 6yr", "LASSO eventglm 7yr")


Zmat.valid <- do.call(cbind, Zvalid)

Zfin.opt <- (Zmat.valid %*% (opt.fit$par) / sum(opt.fit$par))[, 1]
roc.opt <- calc_roc(Zfin.opt, matrix(YYensvalid, ncol = 1), ramp = smoothramp)

with(roc.opt, pseudoloss::calc_auc(fpf,tpf))

sumtable <- data.frame(R.package = c("survival", "randomForestSRC", "CoxBoost", "stats",
                                     "e1071", "class", rep("eventglm", 5),
                                     rep("glmnet", 5)),
                       AUC = unlist(aucalls), coefficient = (opt.fit$par) / sum(opt.fit$par))
rownames(sumtable) <- names.models


print(xtable(sumtable[, 1:2], digits = 3))




## roc and predictiveness curve

roc.opt$cut <- sort(unique(Zfin.opt))

p1<- ggplot(roc.opt, aes(x = fpf, y = tpf)) + geom_line() +
  geom_point(data = roc.opt[seq(5, 675, length.out = 6),],
             aes(x = fpf, y = tpf)) +
  geom_text(data = roc.opt[seq(5, 675, length.out = 6),],
            aes(x = fpf - .05, y = tpf + .05, label = round(cut,2))) +
  theme_bw()+ plotROC::style_roc() + geom_abline(intercept = 0, slope = 1, color = "grey90")


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
ggsave("int-plot.pdf", width = 7.25, height = 3.75)

valid$Zopt.cut <- Zopt.cut
survConcordance(YY ~ Zopt.cut, data = valid)


valid$group.Z <- cut(valid$Zopt.cut,
                     breaks = c(-1, quantile(valid$Zopt.cut, c(0.25, .5, .75)), 2),
                     labels = c("Q1", "Q2", "Q3", "Q4"))


sfitg <- survfit(YY ~ group.Z, data = valid)
meplot <- tidy(sfitg)
meplot$strata <- gsub("group.Z=", "", meplot$strata, fixed = TRUE)

ggplot(meplot, aes(x = time / 365.25, y = estimate, color = strata)) +
  geom_step() +
  geom_step(aes(y = conf.low), linetype = 3) +
  geom_step(aes(y = conf.high), linetype = 3) +
  theme_classic() +
  scale_color_grey("Risk quartile", start = .7, end = .1) +
  xlab("Years since surgery") + ylab("Survival probability")
ggsave("km-fig.pdf", width = 5.25, height = 3.75)

rpm <- replicate(200, {
  coxfit <- coxph(YY ~ I(Zopt.cut * 10), data = valid[sample(1:nrow(valid), nrow(valid), replace = TRUE), ])
  royston(coxfit)[3]
})
sd(rpm)
summary(coxph(YY ~ I(Zopt.cut * 10), data = valid))
summary(coxph(YY ~ group.Z, data = valid))


aucb <- replicate(200, {

  smi <- sample(1:length(Zfin.opt), length(Zfin.opt), replace = TRUE)
  roc.opt <- calc_roc(Zfin.opt[smi], matrix(YYensvalid[smi], ncol = 1), ramp = smoothramp)
  with(roc.opt, pseudoloss::calc_auc(fpf,tpf))

})
sd(aucb)

print(xtable(data.frame(statistic = c("AUC 5-years",
                                      "Survival concordance",
                                      "Explained variation",
                                      "HR per 0.1 units",
                                      "HR Q2 vs Q1",
                                      "HR Q3 vs Q1",
                                      "HR Q4 vs Q1"),
                        value = c(0.718, 0.685, 0.193, 1.44, 1.96, 3.10, 5.84),
                        se = c(0.022, 0.018, 0.032, 0.033, 0.21, 0.20, 0.19))), include.rownames = FALSE)

