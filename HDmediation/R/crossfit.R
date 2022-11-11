crossfit <- function(train, valid, y, x, type = c("binomial", "gaussian"), id = NULL, learners) {
    fit <- regress(train, y, x, id, match.arg(type), learners = learners)
    lapply(valid, function(newX) predictt(fit, newX[, x, drop = FALSE]))
}

regress <- function(train, y, x, id, type, learners) {
    if (!is.null(id)) {
        id <- train[, id]
    }

    family <- ifelse(type == "binomial", binomial(), gaussian())
    fit <- SuperLearner::SuperLearner(
        train[[y]], train[, x, drop = FALSE], family = family[[1]],
        SL.library = learners, id = id,
        method = "method.NNLS", env = environment(SuperLearner::SuperLearner)
    )

    fit
}

predictt <- function(fit, newX) {
    predict(fit, newX)$pred[, 1]
}
