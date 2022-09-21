crossfit <- function(train, valid, y, x, type = c("binomial", "gaussian"), 
                     id = NULL, learners = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.lightgbm", "SL.earth")) {
    fit <- regress(train, y, x, id, match.arg(type), learners = learners)
    lapply(valid, function(newX) predictt(fit, newX[, x]))
}

regress <- function(train, y, x, id, type, learners) {
    if (!is.null(id)) {
        id <- train[, id]
    }
    
    # fit <- glmnet3(train[, x], train[[y]], id = id, family = type)

    family <- ifelse(type == "binomial", binomial(), gaussian())
    fit <- SuperLearner::SuperLearner(
        train[[y]], train[, x], family = family[[1]], SL.library = learners, id = id,
        method = "method.NNLS", env = environment(SuperLearner::SuperLearner)
    )
    
    fit
}

predictt <- function(fit, newX) {
    # predict.glmnet3(fit, newX)
    predict(fit, newX)$pred[, 1]
}
