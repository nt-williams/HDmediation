h_z <- function(data, npsem, folds, learners, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 2)
    colnames(h_z) <- c("h_z(0)", "h_z(1)")
    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid_1 <- valid_0 <- origami::validation(data, folds[[v]])
        try(valid_1[[npsem$S]] <- 0, silent = TRUE)
        try(valid_0[[npsem$S]] <- 0, silent = TRUE)
        valid_1[[npsem$A]] <- 1
        valid_0[[npsem$A]] <- 0

        if (is.null(npsem$Z)) {
            qz <- list(rep(1, nrow(valid_1)), 
                       rep(1, nrow(valid_1)))
            rz <- list(rep(1, nrow(valid_1)), 
                       rep(1, nrow(valid_1)))
            z <- rep(1, nrow(valid_1))
        } else {
            qz <- crossfit(train[, c(npsem$Z, npsem$W, npsem$A, npsem$S)],
                           list(valid_0, valid_1),
                           npsem$Z,
                           "binomial",
                           id = NULL,
                           learners = learners,
                           bound = T)
            
            rz <- crossfit(train[, c(npsem$Z, npsem$M, npsem$W, npsem$A, npsem$S)],
                           list(valid_0, valid_1),
                           npsem$Z,
                           "binomial",
                           id = NULL,
                           learners = learners,
                           bound = T)
            
            z <- valid_1[[npsem$Z]] 
        }

        h_z[folds[[v]]$validation_set, "h_z(0)"] <- 
            (z*qz[[1]] + (1 - z)*(1 - qz[[1]])) / 
            (z*rz[[1]] + (1 - z)*(1 - rz[[1]]))
        h_z[folds[[v]]$validation_set, "h_z(1)"] <- 
            (z*qz[[2]] + (1 - z)*(1 - qz[[2]])) / 
            (z*rz[[2]] + (1 - z)*(1 - rz[[2]]))
    }
    h_z
}
