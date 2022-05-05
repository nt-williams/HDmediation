e <- function(data, npsem, folds, learners) {
    e <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])

        e[folds[[v]]$validation_set, 1] <-
            crossfit(train, list(valid), npsem$A, c(npsem$M, npsem$W), "binomial", learners)[[1]]
    }
    e
}
