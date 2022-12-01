u <- function(data, npsem, bb, hm, aprime, folds, learners, ...) {
    u <- matrix(nrow = nrow(data), ncol = 1)
    colnames(u) <- "u(z,w)"

    data[["b(a',Z,M,W)hm(Z,M,W)"]] <- bb[, gl("b({aprime},Z,M,W)")]*hm

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$A]] <- aprime
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        u[folds[[v]]$validation_set, "u(z,w)"] <-
            crossfit(train, list(valid), "b(a',Z,M,W)hm(Z,M,W)",
                     c(npsem$Z, npsem$A, npsem$W, npsem$S),
                     "gaussian", learners = learners)[[1]]
    }
    u
}

ubar <- function(data, npsem, uu, aprime, folds, learners, ...) {
    ubar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(ubar) <- "ubar(w)"

    data[["u(z,w)"]] <- uu[, "u(z,w)"]

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$A]] <- aprime
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        ubar[folds[[v]]$validation_set, "ubar(w)"] <-
            crossfit(train, list(valid), "u(z,w)",
                     c(npsem$A, npsem$W, npsem$S),
                     "gaussian", learners = learners)[[1]]
        }
    ubar
}
