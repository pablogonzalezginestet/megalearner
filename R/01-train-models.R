#source("00-model-specs.R")
data("gbsg")
data("rotterdam")

mymods <- list(stepwise, random.forest, coxboost,
               direct.binomial, svm, knn,
               function(xx) pseudo.glm(xx, time = 2 * 365.25),
               function(xx) pseudo.glm(xx, time = 4 * 365.25),
               function(xx) pseudo.glm(xx, time = 5 * 365.25),
               function(xx) pseudo.glm(xx, time = 6 * 365.25),
               function(xx) pseudo.glm(xx, time = 7 * 365.25),
               function(xx) pseudo.glmnet(xx, time = 2 * 365.25),
               function(xx) pseudo.glmnet(xx, time = 4 * 365.25),
               function(xx) pseudo.glmnet(xx, time = 5 * 365.25),
               function(xx) pseudo.glmnet(xx, time = 6 * 365.25),
               function(xx) pseudo.glmnet(xx, time = 7 * 365.25))

rotterdam <- subset(rotterdam, nodes > 0)

devel <- data.frame(YY = with(rotterdam, Surv(pmin(rtime, dtime), 1.0 * (recur == 1 | unclass(death) == 1))),
                    age = rotterdam$age, meno = unclass(rotterdam$meno),
                    sizemed = 1.0 * (rotterdam$size == "20-50"),
                    sizebig = 1.0 * (rotterdam$size == ">50"),
                    nodes = unclass(rotterdam$nodes), pgr = unclass(rotterdam$pgr),
                    er = unclass(rotterdam$er), hormon = unclass(rotterdam$hormon))

set.seed(420)
ndex <- 1:nrow(devel)
part <- as.factor(sample(1:10, length(ndex), replace = TRUE))
folds <- split(ndex, part)

Zout <- vector(mode = "list", length = 10)
for(j in 1:length(folds)) {

  training <- devel[unlist(folds[-j]), ]
  validation <- devel[folds[[j]], ]

  Zout[[j]] <- do.call(cbind, lapply(mymods, function(f){

    fhat <- f(training)
    fhat(validation)

  }))

}


fullfits <- lapply(mymods, function(f) {

  f(devel)

})


saveRDS(fullfits, "fullfits.rds")
saveRDS(Zout, "Zout.rds")
saveRDS(folds, "split.rds")
