h_z <- function(data, npsem, aprime, folds, learners, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        train <- stack_data(origami::training(data, folds[[v]]), npsem)
        valid <- origami::validation(data, folds[[v]])
        try(valid[[npsem$S]] <- 0, silent = TRUE)
        valid[[npsem$A]] <- aprime

        # Should S be included?
        p_zmw <- crossfit(train, list(valid), "tmp_tmce_delta",
                          c(npsem$Z, npsem$M, npsem$W, npsem$S, npsem$A),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners)[[1]]

        p_mw <- crossfit(train, list(valid), "tmp_tmce_delta",
                         c(npsem$M, npsem$W, npsem$S, npsem$A),
                         id = "tmp_tmce_id",
                         "binomial", learners = learners)[[1]]

        h_z[folds[[v]]$validation_set, 1] <- (p_zmw / (1 - p_zmw)) * ((1 - p_mw) / p_mw)
    }
    h_z
}

stack_data <- function(data, npsem) {
    delta <- data

    vals_M <- unique(data[, npsem$M, drop = FALSE])
    for (i in 1:nrow(data)) { 
        delta[i, npsem$M] <- vals_M[sample.int(nrow(vals_M), 1), ]
    }
    
    # delta[[npsem$M]] <- sample(x = vals_M,
    #                            size = nrow(data),
    #                            replace = TRUE,
    #                            prob = rep(1 / length(vals_M), length(vals_M)))

    # # sampling from the empirical distribution P_n(m)
    # delta[[npsem$M]] <- sample(data[[npsem$M]], replace = TRUE)

    out <- rbind(data, delta)
    out[["tmp_tmce_delta"]] <- rep(c(0, 1), each = nrow(data))
    out[["tmp_tmce_id"]] <- rep(1:nrow(data), times = 2)
    out
}
