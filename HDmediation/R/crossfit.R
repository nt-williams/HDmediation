crossfit <- function(train, valid, y, x, type = c("binomial", "gaussian"),
                     id = NULL, learners = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.lightgbm", "SL.earth")) {
    fit <- regress(train, y, x, id, match.arg(type), learners = learners)
    lapply(valid, function(newX) predictt(fit, newX[, x, drop = FALSE]))
}

regress <- function(train, y, x, id, type, learners) {
    if (!is.null(id)) {
        id <- train[, id]
    }

    fit <- glmnet3(train[, x, drop = FALSE], train[[y]], id = id, family = type)

    # family <- ifelse(type == "binomial", binomial(), gaussian())
    # fit <- SuperLearner::SuperLearner(
    #     train[[y]], train[, x, drop = FALSE], family = family[[1]], SL.library = learners, id = id,
    #     method = "method.NNLS", env = environment(SuperLearner::SuperLearner)
    # )

    fit
}

predictt <- function(fit, newX) {
    predict.glmnet3(fit, newX)
    # predict(fit, newX)$pred[, 1]
}

glmnet3 <- function(X, y, family = c("gaussian", "binomial"), id = NULL) {
    if (!is.null(id)) {
        # need to match index with fold number for cv.glmnet
        folds <- origami::make_folds(
            nrow(X), fold_fun = origami::folds_vfold,
            cluster_ids = id, V = 10
        )
        
        foldid <- vector("numeric", nrow(X))
        for (i in 1:nrow(X)) {
            for (v in 1:10) {
                if (i %in% folds[[v]]$validation_set) {
                    foldid[i] <- v
                    break
                }
            }
        }
        
    } else {
        foldid <- NULL
    }
    
    ans <- list(covars = names(X))
    
    if (ncol(X) == 1) {
        x <- as.matrix(X)
        ans$fit <- glm(y ~ ., data = cbind(y = y, X), family = match.arg(family))
    } else {
        f <- as.formula(paste0("~ .^", ncol(X)))
        x <- model.matrix(f, X)[, -1]
        ans$fit <- glmnet::cv.glmnet(x, y, family = match.arg(family), foldid = foldid)
    }
    
    # x <- try(model.matrix(f, X)[, -1])
    # if (inherits(x, "try-error")) browser()
    # ans$fit <- glmnet::cv.glmnet(x, y, family = match.arg(family), foldid = foldid)
    
    ans
}

predict.glmnet3 <- function(object, newx) {
    if (inherits(object$fit, "glm")) {
        return(predict(object$fit, newx, type = "response"))
    }
    
    X <- newx[, object$covars, drop = TRUE]
    f <- as.formula(paste0("~ .^", ncol(X)))
    X <- model.matrix(f, X)[, -1]
    as.vector(predict(object$fit, X, type = "response")[, 1])
}
