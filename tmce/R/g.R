g <- function(data, npsem, folds, learners) {
    g <- matrix(nrow = nrow(data), ncol = 1)
    for (t in 1:length(npsem$A)) {
        for (v in seq_along(folds)) {
            train <- origami::training(data, folds[[v]])
            valid <- origami::validation(data, folds[[v]])

            g[folds[[v]]$validation_set, 1] <-
                crossfit(train, list(valid), npsem$A, npsem$W, "binomial", learners)[[1]]
        }
    }
    g
}
