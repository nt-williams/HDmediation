h_z <- function(data, npsem, aprime, folds, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        train <- stack_data(origami::training(data, folds[[v]]), npsem)
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$S]] <- 0
        valid[[npsem$A]] <- aprime

        # Should S be included?
        p_zmw <- crossfit(train, list(valid), "tmp_tmce_delta",
                          c(npsem$Z, npsem$M, npsem$W, npsem$S, npsem$A),
                          id = "tmp_tmce_id",
                          "binomial")[[1]]

        p_mw <- crossfit(train, list(valid), "tmp_tmce_delta",
                         c(npsem$M, npsem$W, npsem$S, npsem$A),
                         id = "tmp_tmce_id",
                         "binomial")[[1]]

        h_z[folds[[v]]$validation_set, 1] <- (p_zmw / (1 - p_zmw)) * ((1 - p_mw) / p_mw)
    }
    h_z
}

stack_data <- function(data, npsem) {
    delta <- data

    vals_M <- unique(data[[npsem$M]])
    delta[[npsem$M]] <- sample(x = vals_M, 
                               size = nrow(data), 
                               replace = TRUE, 
                               prob = rep(1 / length(vals_M), length(vals_M)))
    
    # # sampling from the empirical distribution P_n(m)
    # delta[[npsem$M]] <- sample(data[[npsem$M]], replace = TRUE)

    out <- rbind(data, delta)
    out[["tmp_tmce_delta"]] <- rep(c(0, 1), each = nrow(data))
    out[["tmp_tmce_id"]] <- rep(1:nrow(data), times = 2)
    out
}

# h_z <- function(data, npsem, aprime, folds, ...) {
#     h_z <- matrix(nrow = nrow(data), ncol = 1)
#     qq <- q(data, npsem, aprime, folds)
#     rr <- r(data, npsem, aprime, folds)
#     h_z[, 1] <- qq[, 1] / rr[, 1]
#     h_z
# }
# 
# q <- function(data, npsem, aprime, folds) {
#     qmat <- matrix(nrow = nrow(data), ncol = 1)
#     colnames(qmat) <- c("q(z|a',w)")
#     
#     for (v in seq_along(folds)) {
#         train <- origami::training(data, folds[[v]])
#         valid <- origami::validation(data, folds[[v]])
#         valid[[npsem$S]] <- 0
#         valid[[npsem$A]] <- aprime
#         
#         z <- valid[folds[[v]]$validation_set, npsem$Z]
# 
#         preds <- crossfit(train, list(valid), npsem$Z, c(npsem$A, npsem$W, npsem$S), "binomial")[[1]]
#         qmat[folds[[v]]$validation_set, "q(z|a',w)"] <- (z * preds) + (1 - z) * (1 - preds)
#     }
#     qmat
# }
#     
# r <- function(data, npsem, aprime, folds) {
#     rmat <- matrix(nrow = nrow(data), ncol = 1)
#     colnames(rmat) <- c("r(z|a',m,w)")
#     
#     for (v in seq_along(folds)) {
#         train <- origami::training(data, folds[[v]])
#         valid <- origami::validation(data, folds[[v]])
#         valid[[npsem$S]] <- 0
#         valid[[npsem$A]] <- aprime
#         
#         z <- valid[folds[[v]]$validation_set, npsem$Z]
#         
#         preds <- crossfit(train, list(valid), npsem$Z, c(npsem$A, npsem$M, npsem$W, npsem$S), "binomial")[[1]]
#         rmat[folds[[v]]$validation_set, "r(z|a',m,w)"] <- (z * preds) + (1 - z) * (1 - preds)
#     }
#     rmat
# }   
#     
