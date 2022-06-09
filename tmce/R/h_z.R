h_z <- function(data, npsem, folds, learners, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        train <- stack_data(origami::training(data, folds[[v]]), npsem)
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$S]] <- 0

        # Should S be included?
        p_zmw <- crossfit(train, list(valid), "tmp_tmce_delta",
                          c(npsem$Z, npsem$M, npsem$W, npsem$S), "binomial", learners)[[1]]

        p_mw <- crossfit(train, list(valid), "tmp_tmce_delta",
                         c(npsem$M, npsem$W, npsem$S), "binomial", learners)[[1]]

        h_z[folds[[v]]$validation_set, 1] <- (p_zmw / (1 - p_zmw)) * ((1 - p_mw) / p_mw)
    }
    h_z
}

stack_data <- function(data, npsem) {
    delta <- data
    # sampling from the empirical distribution P_n(m)
    delta[[npsem$M]] <- sample(data[[npsem$M]], replace = TRUE)

    out <- rbind(data, delta)
    out[["tmp_tmce_delta"]] <- rep(c(0, 1), each = nrow(data))
    out
}
