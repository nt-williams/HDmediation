crossfit <- function(train, valid, y, type = c("binomial", "continuous"), id = NULL, learners, bound = FALSE) {
    preds <- mlr3superlearner(data = train,
                              target = y,
                              library = learners,
                              outcome_type = match.arg(type),
                              folds = NULL,
                              newdata = valid,
                              group = id)$preds
    lapply(preds, function(x) bound(x))
}

bound <- function(x, p = 1e-03) {
    pmax(pmin(x, 1 - p), p)
}
