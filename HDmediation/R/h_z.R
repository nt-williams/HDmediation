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
        valid_0$tmp_tmce_id <- valid_1$tmp_tmce_id <- 1:nrow(valid_1)

        p_zmw <- crossfit(train[, c("tmp_tmce_delta", "tmp_tmce_id",
                                    npsem$Z, npsem$M, npsem$W, npsem$A, npsem$S)],
                          list(valid_0, valid_1),
                          "tmp_tmce_delta",
                          "binomial",
                          id = "tmp_tmce_id",
                          learners = learners,
                          bound = T)

        p_mw <- crossfit(train[, c("tmp_tmce_delta", "tmp_tmce_id",
                                   npsem$M, npsem$W, npsem$A, npsem$S)],
                         list(valid_0, valid_1),
                         "tmp_tmce_delta",
                         "binomial",
                         id = "tmp_tmce_id",
                         learners = learners,
                         bound = T)

        h_z[folds[[v]]$validation_set, "h_z(0)"] <-
            (p_zmw[[1]] / (1 - p_zmw[[1]])) * ((1 - p_mw[[1]]) / p_mw[[1]])

        h_z[folds[[v]]$validation_set, "h_z(1)"] <-
            (p_zmw[[2]] / (1 - p_zmw[[2]])) * ((1 - p_mw[[2]]) / p_mw[[2]])
    }
    h_z
}

H_factory <- function(x) {
    apply(x, 2, function(j) {
        is_continuous <- any(schoolmath::is.decimal(j))
        if (is_continuous) {
            return(H_factory_continuous(j))
        }
        H_factory_discrete(j)
    }, simplify = FALSE)
}

uniformly_sample_M <- function(data, M) {
    out <- foreach(m = M,
                   .combine = cbind,
                   .options.future = list(seed = TRUE)) %dofuture% {
                       is_continuous <- any(schoolmath::is.decimal(data[[m]]))
                       if (is_continuous) {
                           return(sample_M_continuous(data[[m]]))
                       }
                       sample_M_discrete(data[[m]])
                   }
    if (length(M) == 1) out <- as.matrix(out)
    colnames(out) <- M
    as.data.frame(out)
}

sample_M_discrete <- function(M) {
    vals_M <- unique(M)
    vals_M[sample.int(length(vals_M), size = length(M), replace = TRUE)]
}

sample_M_continuous <- function(M) {
    minx <- min(M)
    maxx <- max(M)
    runif(length(M), minx, maxx)
}

stack_data <- function(data, npsem) {
    delta <- data
    delta[, npsem$M] <- uniformly_sample_M(data, npsem$M)[, npsem$M]

    # # need to rethink how this would work with continuous variables
    # # sample from a uniform distribution with min and max from the observed data
    # # would be easiest just to sample from the empirical distribution...
    # vals_M <- unique(data[, npsem$M, drop = FALSE])
    # for (i in 1:nrow(data)) {
    #     delta[i, npsem$M] <- vals_M[sample.int(nrow(vals_M), 1), ]
    # }

    out <- rbind(data, delta)
    out[["tmp_tmce_delta"]] <- rep(c(0, 1), each = nrow(data))
    out[["tmp_tmce_id"]] <- rep(1:nrow(data), times = 2)
    out
}
