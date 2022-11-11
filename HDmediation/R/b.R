b <- function(data, npsem, family, folds, learners, ...) {
    b <- matrix(nrow = nrow(data), ncol = 2)
    colnames(b) <- c("b(0,Z,M,W)", "b(1,Z,M,W)")
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        try(train <- train[train[[npsem$S]] == 1, ], silent = TRUE)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[npsem$A]] <- 1
        valid_0[[npsem$A]] <- 0

        preds <- crossfit(train, list(valid_0, valid_1), npsem$Y,
                          c(npsem$W, npsem$A, npsem$Z, npsem$M),
                          family, learners = learners)
        b[folds[[v]]$validation_set, 1] <- preds[[1]]
        b[folds[[v]]$validation_set, 2] <- preds[[2]]
    }
    b
}
