transported <- function(data, A, S, W, Z, M, Y, family, folds = 1) {
    npsem <- Npsem$new(A = A, S = S, W = W, Z = Z, M = M, Y = Y)
    folds <- make_folds(data, folds)

    gg <- g(data, npsem, folds)
    ee <- e(data, npsem, folds)
    cc <- see(data, npsem, folds)
    bb <- b(data, npsem, family, folds)
    t <- 1 - mean(data[[npsem$S]])

    thetas <- eifs <- eif_comps <- list()
    for (param in list(c(1, 1), c(1, 0), c(0, 0))) {
        aprime <- param[1]
        astar <- param[2]

        hz <- h_z(data, npsem, aprime, folds)
        hm <- h_m(hz, gg, ee, aprime, astar)
        uu <- u(data, npsem, bb, hm, aprime, folds)
        uubar <- ubar(data, npsem, uu, aprime, folds)
        vv <- v(data, npsem, bb, hz, aprime, folds)
        vvbar <- vbar(data, npsem, vv, astar, folds)

        # EIF calculation
        S <- data[[npsem$S]]
        A <- data[[npsem$A]]
        Y <- data[[npsem$Y]]

        eify <-
            ((S == 1) & (A == aprime)) /
            (t * gg[, gl("g({aprime}|w)")]) *
            (1 - cc[, gl("c({aprime},z,m,w)")]) / cc[, gl("c({aprime},z,m,w)")] *
            hm * (Y - bb[, gl("b({aprime},Z,M,W)")])

        eifz <-
            ((S == 0) & (A == aprime)) /
            (t * gg[, gl("g({aprime}|w)")]) *
            (uu[, 1] - uubar[, 1])

        eifm <-
            ((S == 0) & (A == astar)) /
            (t * gg[, gl("g({astar}|w)")]) *
            (vv[, 1] - vvbar[, 1])

        eifw <- (S == 0) / t * (vvbar[, 1] - mean(vvbar[S == 0, 1]))

        theta <- mean(eify + eifz + eifm + eifw) + mean(vvbar[S == 0, 1])
        eif <- eify + eifz + eifm + eifw
        eif_comp <- list(y = eify,
                         z = eifz,
                         m = eifm,
                         w = eifw)

        thetas <- c(thetas, list(theta))
        eifs <- c(eifs, list(eif))
        eif_comps <- c(eif_comps, list(eif_comp))
    }

    names(eifs) <- c("11", "10", "00")
    names(thetas) <- c("11", "10", "00")
    names(eif_comps) <- c("11", "10", "00")

    ans <- list(indirect = thetas$`11` - thetas$`10`,
                direct = thetas$`10` - thetas$`00`)

    ans$var_indirect <- var(eifs$`11` - eifs$`10`)
    ans$var_direct <- var(eifs$`10` - eifs$`00`)
    # ans$eif_components <- eif_comps
    ans
}
