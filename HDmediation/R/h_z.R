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
                          c(npsem$Z, npsem$M, npsem$W, npsem$A, npsem$S),
                          id = "tmp_tmce_id",
                          "binomial", learners = learners, bound = TRUE)

        p_mw <- crossfit(train, list(valid_0, valid_1), "tmp_tmce_delta",
                         c(npsem$M, npsem$W, npsem$A, npsem$S),
                         id = "tmp_tmce_id",
                         "binomial", learners = learners, bound = TRUE)

        h_z[folds[[v]]$validation_set, "h_z(0)"] <-
            (p_zmw[[1]] / (1 - p_zmw[[1]])) * ((1 - p_mw[[1]]) / p_mw[[1]])

        h_z[folds[[v]]$validation_set, "h_z(1)"] <-
            (p_zmw[[2]] / (1 - p_zmw[[2]])) * ((1 - p_mw[[2]]) / p_mw[[2]])
    }
    h_z
}

uniformly_sample_M <- function(data, M) {
    out <- foreach(m = M,
                   .combine = cbind,
                   .options.future = list(seed = TRUE)) %dofuture% {
                       sample_M(data[[m]])
                   }
    if (length(M) == 1) out <- as.matrix(out)
    colnames(out) <- M
    as.data.frame(out)
}

sample_M <- function(M) {
    vals_M <- unique(M)
    vals_M[sample.int(length(vals_M), size = length(M), replace = TRUE)]
}

stack_data <- function(data, npsem) {
    delta <- data
    delta[, npsem$M] <- uniformly_sample_M(data, npsem$M)[, npsem$M]
    out <- rbind(data, delta)
    out[["tmp_tmce_delta"]] <- rep(c(0, 1), each = nrow(data))
    out[["tmp_tmce_id"]] <- rep(1:nrow(data), times = 2)
    out
}
