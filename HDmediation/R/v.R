v <- function(data, npsem, bb, hz, aprime, folds, learners, ...) {
    v <- matrix(nrow = nrow(data), ncol = 1)
    colnames(v) <- "v(m,w)"

    data[["tmp_HDmediation_outcome_v_fit"]] <- bb[, gl("b({aprime},Z,M,W)")]*hz[, gl("h_z({aprime})")]

    for (fold in seq_along(folds)) {
        train <- origami::training(data, folds[[fold]])
        valid <- origami::validation(data, folds[[fold]])
        valid[[npsem$A]] <- aprime
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        v[folds[[fold]]$validation_set, "v(m,w)"] <-
            crossfit(train[, c("tmp_HDmediation_outcome_v_fit",
                               npsem$M, npsem$A, npsem$W, npsem$S)],
                     list(valid),
                     "tmp_HDmediation_outcome_v_fit",
                     "continuous",
                     learners = learners)[[1]]
    }
    v
}

vbar <- function(data, npsem, vv, astar, folds, learners, ...) {
    vbar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(vbar) <- "vbar(w)"

    data[["tmp_HDmediation_outcome_v_fit"]] <- vv[, "v(m,w)"]

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$A]] <- astar
        try(valid[[npsem$S]] <- 0, silent = TRUE)

        vbar[folds[[v]]$validation_set, "vbar(w)"] <-
            crossfit(train[, c("tmp_HDmediation_outcome_v_fit",
                               npsem$A, npsem$W, npsem$S)],
                     list(valid),
                     "tmp_HDmediation_outcome_v_fit",
                     "continuous",
                     learners = learners)[[1]]
    }
    vbar
}
