g <- function(data, npsem, folds, learners) {
    gmat <- matrix(nrow = nrow(data), ncol = 2)
    colnames(gmat) <- c("g(0|w)", "g(1|w)")

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        preds <- crossfit(train, list(valid), npsem$A, c(npsem$W, npsem$S),
                          "binomial", learners = learners, bound = TRUE)[[1]]

        gmat[folds[[v]]$validation_set, "g(0|w)"] <- 1 - preds
        gmat[folds[[v]]$validation_set, "g(1|w)"] <- preds
    }
    gmat
}
