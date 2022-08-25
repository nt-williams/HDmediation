v <- function(data, npsem, bb, hz, aprime, folds, ...) {
    v <- matrix(nrow = nrow(data), ncol = 1)
    colnames(v) <- "v(m,w)"

    data[["b(a',Z,M,W)hz(Z,M,W)"]] <- bb[, gl("b({aprime},Z,M,W)")]*hz[, 1]

    for (fold in seq_along(folds)) {
        train <- origami::training(data, folds[[fold]])
        valid <- origami::validation(data, folds[[fold]])
        valid[[npsem$A]] <- aprime
        valid[[npsem$S]] <- 0

        v[folds[[fold]]$validation_set, "v(m,w)"] <-
            crossfit(train, list(valid), "b(a',Z,M,W)hz(Z,M,W)",
                     c(npsem$M, npsem$A, npsem$W, npsem$S), "gaussian")[[1]]
    }
    v
}

vbar <- function(data, npsem, vv, astar, folds, ...) {
    vbar <- matrix(nrow = nrow(data), ncol = 1)
    colnames(vbar) <- "vbar(w)"

    data[["v(m,w)"]] <- vv[, "v(m,w)"]

    for (v in seq_along(folds)) {
        train <- origami::training(data, folds[[v]])
        valid <- origami::validation(data, folds[[v]])
        valid[[npsem$A]] <- astar
        valid[[npsem$S]] <- 0

        vbar[folds[[v]]$validation_set, "vbar(w)"] <-
            crossfit(train, list(valid), "v(m,w)",
                     c(npsem$A, npsem$W, npsem$S), "gaussian")[[1]]
    }
    vbar
}

true_vbar <- function(data, aprime, astar) {
    intv <- function(m, w, aprime) {
        my(m, 1, aprime, w) * pz(1, aprime, w, 0) + my(m, 0, aprime, w) * pz(0, aprime, w, 0)
    }
    
    my <- function(m, z, a, w) {
        plogis(-log(5) + log(8) * z  + log(4) * m - log(1.2) * w[, "W1"] + log(1.2) * w[, "W1"] * z)
    }
    
    pmaw <- function(m, a, w, s) {
        pm(m, 1, a, w, s) * pz(1, a, w, s) + pm(m, 0, a, w, s ) * pz(0, a, w, s)
    }
    
    pz <- function(z, a, w, s) {
        prob1 <- plogis(-log(2) + (log(4) * a) - log(2) * w[, "W1"] + log(1.4) * s + log(1.43) * s * a)
        z * prob1 + (1 - z) * (1 - prob1)
    }
    
    pm <- function(m, z, a, w, s) {
        prob1 <- plogis(-log(2) + log(4) * z - log(1.4) * w[, "W1"] + log(1.4) * s)
        m * prob1 + (1 - m) * (1 - prob1)
    }
    
    w <- data[, c("W0", "W1")]
    intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)
}

true_v <- function(data, aprime) {
    intv <- function(m, w, aprime) {
        my(m, 1, aprime, w) * pz(1, aprime, w, 0) + my(m, 0, aprime, w) * pz(0, aprime, w, 0)
    }

    my <- function(m, z, a, w) {
        plogis(-log(5) + log(8) * z  + log(4) * m - log(1.2) * w[, "W1"] + log(1.2) * w[, "W1"] * z)
    }

    pz <- function(z, a, w, s) {
        prob1 <- plogis(-log(2) + (log(4) * a) - log(2) * w[, "W1"] + log(1.4) * s + log(1.43) * s * a)
        z * prob1 + (1 - z) * (1 - prob1)
    }

    m <- data$M
    w <- data[, c("W0", "W1")]

    vv <- matrix(intv(m, w, aprime), ncol = 1)
    colnames(vv) <- "v(m,w)"
    vv
}
