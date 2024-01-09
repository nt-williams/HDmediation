crossfit <- function(train, valid, y, type = c("binomial", "continuous"), id = NULL, learners, bound = FALSE) {
    # fit <- hal9001::fit_hal(train[, setdiff(names(train), y), 
    #                               drop = F], train[[y]], 
    #                  family = ifelse(match.arg(type) == "binomial", "binomial", "gaussian"))
    # lapply(valid, function(x) predict(fit, x))
    preds <- mlr3superlearner(data = train,
                              target = y,
                              library = learners,
                              outcome_type = match.arg(type),
                              folds = NULL,
                              newdata = valid,
                              group = id)$preds
    preds
    # # lapply(preds, function(x) bound(x))
}

bound <- function(x, p = 1e-03) {
    pmax(pmin(x, 1 - p), p)
}
