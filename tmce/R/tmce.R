tmce <- function(data, A, S, W, Z, M, Y, learners, family, folds = 1) {
    checkmate::assertDataFrame(data[, c(A, S, W, Z, M, Y)])
    checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)

    npsem <- Npsem$new(A = A, S = S, W = W, Z = Z, M = M, Y = Y)
    folds <- make_folds(data, folds)

    gg <- g(data, npsem, folds, learners)
    ee <- e(data, npsem, folds, learners)
    cc <- see(data, npsem, folds, learners)
    bb <- b(data, npsem, "binomial", folds, learners)
    hz <- h_z(data, npsem, folds, learners)
    t <- 1 - mean(data[[npsem$S]])

    eifs <- list()
    for (param in list(c(1, 1), c(1, 0), c(0, 0))) {
        aprime <- param[1]
        astar <- param[2]

        hm <- h_m(hz, gg, ee, aprime, astar)
        uu <- u(data, npsem, bb, hm, aprime, folds, learners)
        uubar <- ubar(data, npsem, uu, aprime, folds, learners)
        vv <- v(data, npsem, bb, hz, aprime, folds, learners)
        vvbar <- vbar(data, npsem, vv, astar, folds, learners)

        # EIF calculation
        D_PY <-
            ((data[[npsem$S]] == 1) & (data[[npsem$A]] == aprime)) /
            (t * gg[, gl("g({aprime}|w)")]) *
            (1 - cc[, 1]) /
            cc[, 1] *
            hm * (data[[npsem$Y]] - bb[, gl("b({aprime},Z,M,W)")])

        D_PZ <-
            ((data[[npsem$S]] == 0) & (data[[npsem$A]] == aprime)) /
            (t * gg[, gl("g({aprime}|w)")]) *
            (uu[, 1] - uubar[, 1])

        D_PM <-
            ((data[[npsem$S]] == 0) & (data[[npsem$A]] == astar)) /
            (t * gg[, gl("g({astar}|w)")]) *
            (vv[, 1] - vvbar[, 1])

        eifs <- c(eifs, list(D_PY + D_PZ + D_PM + (((data[[npsem$S]] == 0) / t) * vvbar[, 1])))
    }

    names(eifs) <- c("11", "10", "00")
    eifs
}
