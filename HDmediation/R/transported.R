transported <- function(data, A, S, W, Z, M, Y, family, folds = 1,
                        learners_g = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_e = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_c = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_b = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_hz = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_u = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_ubar = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_v = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                        learners_vbar = c("SL.glm", "SL.glm.interaction", "SL.mean")) {
    npsem <- Npsem$new(A = A, S = S, W = W, Z = Z, M = M, Y = Y)
    folds <- make_folds(data, folds)

    gg <- g(data, npsem, folds, learners_g)
    ee <- e(data, npsem, folds, learners_e)
    cc <- see(data, npsem, folds, learners_c)
    bb <- b(data, npsem, family, folds, learners_b)
    t <- 1 - mean(data[[npsem$S]])
    hz <- h_z(data, npsem, folds, learners_hz)

    thetas <- eifs <- vector("list", 3)
    vvbar <- matrix(nrow = nrow(data), ncol = 3)
    colnames(vvbar) <- c("00", "10", "11")
    names(eifs) <- c("11", "10", "00")
    names(thetas) <- c("11", "10", "00")
    for (param in list(c(1, 1), c(1, 0), c(0, 0))) {
        aprime <- param[1]
        astar <- param[2]

        hm <- h_m(hz, gg, ee, aprime, astar)
        uu <- u(data, npsem, bb, hm, aprime, folds, learners_u)
        uubar <- ubar(data, npsem, uu, aprime, folds, learners_ubar)
        vv <- v(data, npsem, bb, hz, aprime, folds, learners_v)
        vvbar[, paste(param, collapse = "")] <-
            vbar(data, npsem, vv, astar, folds, learners_vbar)

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
            (vv[, 1] - vvbar[, paste(param, collapse = "")])

        eifw <- (S == 0) / t * (vvbar[, paste(param, collapse = "")] -
                                    mean(vvbar[S == 0, paste(param, collapse = "")]))

        theta <- mean(eify + eifz + eifm + eifw) + mean(vvbar[S == 0, paste(param, collapse = "")])
        eif <- eify + eifz + eifm + eifw

        thetas[[paste(param, collapse = "")]] <- theta
        eifs[[paste(param, collapse = "")]] <- eif
    }

    ans <- data.frame(indirect = thetas$`11` - thetas$`10`,
                      direct = thetas$`10` - thetas$`00`,
                      gcomp_indirect = mean(vvbar[, "11"] - vvbar[, "10"]),
                      gcomp_direct = mean(vvbar[, "10"] - vvbar[, "00"]))

    ans$var_indirect <- var(eifs$`11` - eifs$`10`)
    ans$var_direct <- var(eifs$`10` - eifs$`00`)

    ci_indirect <- ans$indirect + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_indirect / nrow(data))
    ci_direct <- ans$direct + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_direct / nrow(data))

    ans$ci_indirect_low <- ci_indirect[1]
    ans$ci_indirect_high <- ci_indirect[2]
    ans$ci_direct_low <- ci_direct[1]
    ans$ci_direct_high <- ci_direct[2]

    ans
}
