glmnet3 <- function(X, y, family = c("gaussian", "binomial"), id = NULL) {
    if (!is.null(id)) {
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

SL.glmnet3 <- function(Y, X, newX, family, obsWeights, id, nrounds = 1000, verbose = -1,
                       learning_rate = 0.1, min_data_in_leaf = 10, max_depth = -1, ...) {
    if (!requireNamespace("glmnet", quietly = FALSE)) {
        stop("loading required package (glmnet) failed", call. = FALSE)
    }

    model <- glmnet3(X, Y, id = id, family = family$family)

    pred <- predict.glmnet3(model, newX)
    fit <- list(object = model)
    class(fit) <- c("SL.glmnet3")
    out <- list(pred = pred, fit = fit)
    return(out)
}

predict.SL.glmnet3 <- function(object, newdata, ...) {
    if (!requireNamespace("glmnet", quietly = FALSE)) {
        stop("loading required package (glmnet) failed", call. = FALSE)
    }

    pred <- predict.glmnet3(object$object, newdata)
    return(pred)
}
