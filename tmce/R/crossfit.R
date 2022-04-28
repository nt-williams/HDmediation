crossfit <- function(train, valid, y, x, type = c("binomial", "continuous"), learners) {
    fit <- regress(train, y, x, match.arg(type), learners)
    lapply(valid, function(x) predictt(fit, x))
}

regress <- function(train, y, x, type, learners) {
    family <- ifelse(type == "binomial", binomial(), gaussian())
    fit <- SuperLearner::SuperLearner(
        train[[y]], train[, x], family = family[[1]], SL.library = learners,
        method = "method.NNLS", env = environment(SuperLearner::SuperLearner)
    )

    fit
}

predictt <- function(fit, data) {
    predict(fit, data[, fit$varNames])$pred[, 1]
}
