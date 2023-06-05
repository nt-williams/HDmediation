u <- function(data, npsem, bb, hm, aprime, folds, learners, ...) {
    u <- matrix(nrow = nrow(data), ncol = 1)
    colnames(u) <- "u(z,w)"

    data[["tmp_HDmediation_outcome_u_fit"]] <- bb[, gl("b({aprime},Z,M,W)")]*hm

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$A]] <- aprime
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        u[folds[[v]]$validation_set, "u(z,w)"] <-
            crossfit(train[, c("tmp_HDmediation_outcome_u_fit",
                               npsem$Z, npsem$A, npsem$W, npsem$S)],
                     list(valid),
                     "tmp_HDmediation_outcome_u_fit",
                     "continuous",
                     learners = learners)[[1]]
    }
    u
}

ubar <- function(data, npsem, uu, aprime, folds, learners, ...) {
    ubar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(ubar) <- "ubar(w)"

    data[["tmp_HDmediation_outcome_u_fit"]] <- uu[, "u(z,w)"]

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$A]] <- aprime
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        ubar[folds[[v]]$validation_set, "ubar(w)"] <-
            crossfit(train[, c("tmp_HDmediation_outcome_u_fit",
                               npsem$A, npsem$W, npsem$S)],
                     list(valid),
                     "tmp_HDmediation_outcome_u_fit",
                     "continuous",
                     learners = learners)[[1]]
        }
    ubar
}
