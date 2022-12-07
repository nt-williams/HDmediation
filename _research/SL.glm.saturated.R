SL.glm.saturated <- function(Y, X, newX, family, obsWeights, ...) {
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    f <- as.formula(paste0("Y ~ .^", ncol(X)))
    fit.glm <- glm(f, data = X, family = family, weights = obsWeights)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}
