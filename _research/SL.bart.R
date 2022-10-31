SL.bart2 <- function(Y, X, newX, family, obsWeights, id, n.samples = 500L,
                     n.burn = 500L, n.chains = 4L, verbose = FALSE) {
    if (!requireNamespace("dbarts", quietly = FALSE)) {
        stop("loading required package (dbarts) failed", call. = FALSE)
    }

    predict_bart <- function(object, newdata) {
        mat <- predict(object, newdata, type = "response")
        apply(mat, 2, mean)
    }

    bart_data <- cbind(Y = Y, X)
    form <- Y ~ .

    model <- dbarts::bart2(formula = form, data = bart_data,
                           n.samples = n.samples, n.burn = n.burn,
                           n.chains = n.chains, verbose = verbose,
                           keepTrees = TRUE)

    pred <- predict_bart(model, newX)
    fit <- list(object = model)
    class(fit) <- c("SL.bart2")
    out <- list(pred = pred, fit = fit)
    return(out)
}

predict.SL.bart2 <- function(object, newdata, family, ...) {
    if (!requireNamespace("dbarts", quietly = FALSE)) {
        stop("loading required package (dbarts) failed", call. = FALSE)
    }

    predict_bart <- function(object, newdata) {
        mat <- predict(object, newdata, type = "response")
        apply(mat, 2, mean)
    }

    pred <- predict_bart(object$object, newdata)
    return(pred)
}
