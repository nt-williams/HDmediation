crossfit <- function(train, valid, y, x, type = c("binomial", "gaussian"), id = NULL, learners, bound = FALSE) {
    fit <- regress(train, y, x, id, match.arg(type), learners = learners)
    lapply(valid, function(newX) predictt(fit, newX[, x, drop = FALSE], bound))
}

regress <- function(train, y, x, id, type, learners) {
    if (!is.null(id)) {
        id <- train[, id]
    }

    family <- ifelse(type == "binomial", binomial(), gaussian())
    fit <- SuperLearner::SuperLearner(
        train[[y]], train[, x, drop = FALSE], family = family[[1]],
        SL.library = learners, id = id,
        method = "method.NNLS", env = environment(SuperLearner::SuperLearner),
        cvControl = SuperLearner::SuperLearner.CV.control(V = 5L)
    )

    fit
}

predictt <- function(fit, newX, bound = FALSE) {
    # predict(fit, newX)
    if (!bound) {
        return(predict(fit, newX)$pred[, 1])
    }
    bound(predict(fit, newX)$pred[, 1])
}

bound <- function(x, p = 1e-03) {
    pmax(pmin(x, 1 - p), p)
}
