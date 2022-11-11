h_z <- function(data, npsem, folds, learners, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_z) <- c("h_z(0)", "h_z(1)")
    for (v in seq_along(folds)) {
        train <- stack_data(origami::training(data, folds[[v]]), npsem)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        try(valid_1[[npsem$S]] <- 0, silent = TRUE)
        try(valid_0[[npsem$S]] <- 0, silent = TRUE)
        valid_1[[npsem$A]] <- 1
        valid_0[[npsem$A]] <- 0

        p_zmw <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                          c(npsem$Z, npsem$M, npsem$W, npsem$S, npsem$A),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners)

        p_mw <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                         c(npsem$M, npsem$W, npsem$S, npsem$A),
                         id = "tmp_tmce_id",
                         "binomial", learners = learners)

        h_z[folds[[v]]$validation_set, "h_z(0)"] <-
            (p_zmw[[1]] / (1 - p_zmw[[1]])) * ((1 - p_mw[[1]]) / p_mw[[1]])

        h_z[folds[[v]]$validation_set, "h_z(1)"] <-
            (p_zmw[[2]] / (1 - p_zmw[[2]])) * ((1 - p_mw[[2]]) / p_mw[[2]])
    }
    h_z
}

stack_data <- function(data, npsem) {
    delta <- data

    # need to rethink how this would work with countinuous variables
    # sample from a uniform distribution with min and max from the observed data
    # would be easiest just to sample from the empirical distribution...
    vals_M <- unique(data[, npsem$M, drop = FALSE])
    for (i in 1:nrow(data)) {
        delta[i, npsem$M] <- vals_M[sample.int(nrow(vals_M), 1), ]
    }

    # delta[[npsem$M]] <- sample(x = vals_M,
    #                            size = nrow(data),
    #                            replace = TRUE,
    #                            prob = rep(1 / length(vals_M), length(vals_M)))

    # # sampling from the empirical distribution P_n(m)
    # delta[, npsem$M] <- sample(data[, npsem$M, drop = FALSE], replace = TRUE)

    out <- rbind(data, delta)
    out[["tmp_tmce_delta"]] <- rep(c(0, 1), each = nrow(data))
    out[["tmp_tmce_id"]] <- rep(1:nrow(data), times = 2)
    out
}
