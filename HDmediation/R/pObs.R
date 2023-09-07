pObs <- function(data, npsem, folds, learners, ...) {
    probs <- matrix(nrow = nrow(data), ncol = 2)
    colnames(probs) <- c("P(delta=1|A=0,Z,M,W)", "P(delta=1|A=1,Z,M,W)")
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        try(train <- train[train[[npsem$S]] == 1, ], silent = TRUE)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[npsem$A]] <- 1
        valid_0[[npsem$A]] <- 0

        preds <- crossfit(train[, c(npsem$cens, npsem$W, npsem$A, npsem$Z, npsem$M)],
                          list(valid_0, valid_1),
                          npsem$cens,
                          "binomial",
                          learners = learners)

        probs[folds[[v]]$validation_set, 1] <- preds[[1]]
        probs[folds[[v]]$validation_set, 2] <- preds[[2]]
    }
    probs
}
