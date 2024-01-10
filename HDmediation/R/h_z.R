h_z <- function(data, npsem, folds, learners, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_z) <- c("h_z(0)", "h_z(1)")
    for (v in seq_along(folds)) {
        train <- stack_data(origami::training(data, folds[[v]]), npsem, 5)
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        try(valid_1[[npsem$S]] <- 0, silent = TRUE)
        try(valid_0[[npsem$S]] <- 0, silent = TRUE)
        valid_1[[npsem$A]] <- 1
        valid_0[[npsem$A]] <- 0
        valid_0$tmp_tmce_id <- valid_1$tmp_tmce_id <- 1:nrow(valid_1)
        # browser()
        
        # p_zmw <- list()
        # fit_zmw <- hal9001::fit_hal(X = as.matrix(train[, c(npsem$Z, npsem$W, npsem$A, npsem$S)]),
        #                  Y = train[["tmp_tmce_delta"]],
        #                  X_unpenalized = as.matrix(train[, npsem$M]),
        #                  max_degree = 2,
        #                  family = "binomial",
        #                  id = train[["tmp_tmce_id"]])
        # p_zmw[[1]] <- predict(fit_zmw,
        #                       valid_0[, c(npsem$Z, npsem$W, npsem$A, npsem$S)],
        #                       as.matrix(valid_0[, npsem$M]))
        # p_zmw[[2]] <- predict(fit_zmw,
        #                       valid_1[, c(npsem$Z, npsem$W, npsem$A, npsem$S)],
        #                       as.matrix(valid_1[, npsem$M]))
        # 
        # 
        # p_mw <- list()
        # fit_mw <- hal9001::fit_hal(X = as.matrix(train[, c(npsem$W, npsem$A, npsem$S)]),
        #                  Y = train[["tmp_tmce_delta"]],
        #                  X_unpenalized = as.matrix(train[, npsem$M]),
        #                  max_degree = 2,
        #                  family = "binomial",
        #                  id = train[["tmp_tmce_id"]])
        # p_mw[[1]] <- predict(fit_mw,
        #                       valid_0[, c(npsem$W, npsem$A, npsem$S)],
        #                       as.matrix(valid_0[, npsem$M]))
        # p_mw[[2]] <- predict(fit_mw,
        #                       valid_1[, c(npsem$W, npsem$A, npsem$S)],
        #                       as.matrix(valid_1[, npsem$M]))

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
        
        num <- (nrow(data) / (nrow(data)*5))*(p_zmw[[1]] / (1 - p_zmw[[1]])) 
        denom <- (nrow(data) / (nrow(data)*5))*(p_mw[[1]] / (1 - p_mw[[1]]))

        h_z[folds[[v]]$validation_set, "h_z(0)"] <- num / denom
            # (p_zmw[[1]] / (1 - p_zmw[[1]])) * ((1 - p_mw[[1]]) / p_mw[[1]])
        
        
        num <- (nrow(data) / (nrow(data)*5))*(p_zmw[[2]] / (1 - p_zmw[[2]])) 
        denom <- (nrow(data) / (nrow(data)*5))*(p_mw[[2]] / (1 - p_mw[[2]]))

        h_z[folds[[v]]$validation_set, "h_z(1)"] <- num / denom
            # (p_zmw[[2]] / (1 - p_zmw[[2]])) * ((1 - p_mw[[2]]) / p_mw[[2]])
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

uniformly_sample_M <- function(data, M, A, W, ndraws = 1) {
    out <- foreach(m = M,
                   .combine = cbind,
                   .options.future = list(seed = TRUE)) %dofuture% {
                       is_continuous <- any(schoolmath::is.decimal(data[[m]]))
                       if (is_continuous) {
                           return(sample_M_continuous(m, data[, c(m, A, W)], ndraws))
                       }
                       sample_M_discrete(data[[m]], ndraws)
                   }
    if (length(M) == 1) out <- as.matrix(out)
    colnames(out) <- M
    as.data.frame(out)
}

sample_M_discrete <- function(M, ndraws) {
    vals_M <- unique(M)
    vals_M[sample.int(length(vals_M), size = length(M)*ndraws, replace = TRUE)]
}

# sample_M_continuous <- function(M) {
#     minx <- min(M)
#     maxx <- max(M)
#     runif(length(M), minx, maxx)
# }

sample_M_continuous <- function(M, data, ndraws) {
    aw <- paste(setdiff(names(data), M), collapse = "+")
    f <- as.formula(paste0(M, "~", aw))
    mu <- predict(lm(f, data = data), data)
    rnorm(length(data[[M]])*ndraws, mu, sd(data[[M]]))
}

stack_data <- function(data, npsem, ndraws = 1) {
    delta <- purrr::map_dfr(1:ndraws, function(x) data)
    delta[, npsem$M] <- 
        uniformly_sample_M(data, npsem$M, npsem$A, npsem$W, ndraws)[, npsem$M]
    data[["tmp_tmce_delta"]] <- rep(0, nrow(data))
    delta[["tmp_tmce_delta"]] <- rep(1, nrow(data)*ndraws)
    out <- rbind(data, delta)
    out[["tmp_tmce_id"]] <- rep(1:nrow(data), times = ndraws + 1)
    out
}
