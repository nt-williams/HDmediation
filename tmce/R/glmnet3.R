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

    f <- as.formula(paste0("~ .^", ncol(X)))
    X <- model.matrix(f, X)[, -1]

    ans$fit <- glmnet::cv.glmnet(X, y, family = match.arg(family), foldid = foldid)#,
                                 #lambda = seq(1 / nrow(X)^2, 1 / sqrt(nrow(X)), length.out = 50))

    ans
}

predict.glmnet3 <- function(object, newx) {
    X <- newx[, object$covars]
    f <- as.formula(paste0("~ .^", ncol(X)))
    X <- model.matrix(f, X)[, -1]
    as.vector(predict(object$fit, X, type = "response", s = "lambda.min")[, 1])
}
