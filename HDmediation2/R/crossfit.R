crossfit <- function(train, valid, y, type = c("binomial", "continuous"), id = NULL, learners, bound = FALSE) {
    # if(h2o::h2o.clusterIsUp()) {
    #     h2o::h2o.shutdown(prompt = FALSE)
    #     Sys.sleep(3) # Wait for shutdown to complete
    # }

    nfolds <- mlr3superlearner:::set_folds(nrow(train), match.arg(type), train[[y]], FALSE)

    if (match.arg(type) == "binomial") {
        train[[y]] <- as.factor(train[[y]])
    }

    # if (match.arg(type) == "continuous") browser()

    if (is.null(id)) {
        train$tmp_cv_folds <- origami::folds2foldvec(make_folds(train, nfolds, NULL))
    } else {
        train$tmp_cv_folds <- origami::folds2foldvec(make_folds(train, nfolds, train[[id]]))
    }

    train <- h2o::as.h2o(train)
    valid <- lapply(valid, h2o::as.h2o)

    fit <- h2o::h2o.automl(
        x = setdiff(names(train), c(y, id, "tmp_cv_folds")),
        y = y,
        training_frame = train,
        fold_column = "tmp_cv_folds",
        max_models = 10,
        balance_classes = match.arg(type) == "binomial",
        sort_metric = ifelse(match.arg(type) == "binomial", "logloss", "MSE"),
        distribution = ifelse(match.arg(type) == "binomial", "bernoulli", "gaussian"),
        stopping_metric = ifelse(match.arg(type) == "binomial", "logloss", "MSE")
    )

    if (match.arg(type) == "binomial") {
        preds <- lapply(valid, \(x) as.data.frame(predict(fit, x))$p1)
    } else {
        preds <- lapply(valid, \(x) as.vector(predict(fit, x)))
    }

    gc()
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()

    lapply(preds, function(x) bound(x))
}

bound <- function(x, p = 1e-03) {
    pmax(pmin(x, 1 - p), p)
}
