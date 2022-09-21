crossfit <- function(train, valid, y, x, type = c("binomial", "gaussian"), id = NULL) {# , learners) {
    fit <- regress(train, y, x, id, match.arg(type))
    lapply(valid, function(newX) predictt(fit, newX[, x]))
}

regress <- function(train, y, x, id, type) { #, learners) {
    if (!is.null(id)) {
        id <- train[, id]
    }

    fit <- glmnet3(train[, x], train[[y]], id = id, family = type)

    # fit <- sal(train[, x], train[[y]], id = id, family = type, ratio = 0.2)

    # family <- ifelse(type == "binomial", binomial(), gaussian())
    # fit <- SuperLearner::SuperLearner(
    #     train[[y]], train[, x], family = family[[1]], SL.library = learners, id = id,
    #     method = "method.NNLS", env = environment(SuperLearner::SuperLearner)
    # )
    
    # bart_data <- cbind(Y = train[[y]], train[, x])
    # form <- Y ~ .
    # 
    # fit <- dbarts::bart2(formula = form, data = bart_data, keepTrees = TRUE, verbose = FALSE)
    fit
}

predictt <- function(fit, newX) {
    predict.glmnet3(fit, newX)
    # predict.sal(fit, newX)
    # predict(fit, data[, fit$varNames])$pred[, 1]
    # mat <- predict(fit, newX, type = "response")
    # apply(mat, 2, mean)
}
