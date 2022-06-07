see <- function(data, npsem, folds, learners, ...) {
    cmat <- matrix(nrow = nrow(data), ncol = 1)
    colnames(cmat) <- "c(a,z,m,w)"

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])

        preds <- crossfit(train, list(valid), npsem$S,
                          c(npsem$A, npsem$Z, npsem$M, npsem$W), "binomial", learners)[[1]]
        cmat[folds[[v]]$validation_set, "c(a,z,m,w)"] <- preds
    }
    cmat
}
