e <- function(data, npsem, folds, learners) {
    emat <- matrix(nrow = nrow(data), ncol = 2)
    colnames(emat) <- c("e(0|m,w)", "e(1|m,w)")

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        preds <- crossfit(train, list(valid), npsem$A, c(npsem$M, npsem$W, npsem$S),
                          "binomial", learners = learners, bound = TRUE)[[1]]
        emat[folds[[v]]$validation_set, "e(0|m,w)"] <- 1 - preds
        emat[folds[[v]]$validation_set, "e(1|m,w)"] <- preds
    }

    emat
}
