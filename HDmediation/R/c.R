see <- function(data, npsem, folds, learners, ...) {
    cmat <- matrix(nrow = nrow(data), ncol = 2)
    colnames(cmat) <- c("c(0,z,m,w)", "c(1,z,m,w)")

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        valid_1[[npsem$A]] <- 1
        valid_0[[npsem$A]] <- 0

        preds <- crossfit(train, list(valid_0, valid_1), npsem$S,
                          c(npsem$A, npsem$Z, npsem$M, npsem$W), "binomial",
                          learners = learners, bound = TRUE)
        cmat[folds[[v]]$validation_set, 1] <- preds[[1]]
        cmat[folds[[v]]$validation_set, 2] <- preds[[2]]
    }
    cmat
}
