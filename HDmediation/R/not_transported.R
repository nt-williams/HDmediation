not_transported <- function(data, A, W, Z, M, Y, family, folds = 1, partial_tmle, bounds,
                            learners_g = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_e = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_b = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_hz = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_u = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_ubar = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_v = c("SL.glm", "SL.glm.interaction", "SL.mean"),
                            learners_vbar = c("SL.glm", "SL.glm.interaction", "SL.mean")) {
    npsem <- Npsem$new(A = A, W = W, Z = Z, M = M, Y = Y)
    folds <- make_folds(data, folds)

    bounds <- scale_y(data[[npsem$Y]], family, bounds)
    data[[npsem$Y]] <- Y <- bounds$y
    A <- data[[npsem$A]]

    gg <- g(data, npsem, folds, learners_g)
    ee <- e(data, npsem, folds, learners_e)
    bb <- b(data, npsem, family, folds, learners_b)
    hz <- h_z(data, npsem, folds, learners_hz)

    thetas <- eifs <- list()
    vvbar <- matrix(nrow = nrow(data), ncol = 3)
    colnames(vvbar) <- c("00", "10", "11")
    for (param in list(c(1, 1), c(1, 0), c(0, 0))) {
        aprime <- param[1]
        astar <- param[2]

        hm <- h_m(hz, gg, ee, aprime, astar)

        ipwy <- (A == aprime) / gg[, gl("g({aprime}|w)")]
        if (partial_tmle) {
            fit <- glm(Y ~ 1, offset = qlogis(bb[, gl("b({aprime},Z,M,W)")]), family = "binomial",
                       subset = A == aprime, weights = ipwy * hm / mean(ipwy * hm))
            bb[, gl("b({aprime},Z,M,W)")] <- plogis(coef(fit) + qlogis(bb[, gl("b({aprime},Z,M,W)")]))
        }

        uu <- u(data, npsem, bb, hm, aprime, folds, learners_u)
        uubar <- ubar(data, npsem, uu, aprime, folds, learners_ubar)
        vv <- v(data, npsem, bb, hz, aprime, folds, learners_v)
        vvbar[, paste(param, collapse = "")] <- vbar(data, npsem, vv, astar, folds, learners_vbar)

        # EIF calculation
        eify <- ipwy * hm / mean(ipwy * hm) * (Y - bb[, gl("b({aprime},Z,M,W)")])

        ipwz <- (A == aprime) / gg[, gl("g({aprime}|w)")]
        eifz <- ipwz / mean(ipwz) * (uu[, 1] - uubar[, 1])

        ipwm <- (A == astar) / gg[, gl("g({astar}|w)")]
        eifm <- ipwm / mean(ipwm) * (vv[, 1] - vvbar[, paste(param, collapse = "")])

        eif <- rescale_y(eify + eifz + eifm + vvbar[, paste(param, collapse = "")], bounds$bounds)
        theta <- mean(eif)

        thetas <- c(thetas, list(theta))
        eifs <- c(eifs, list(eif))
    }

    names(eifs) <- c("11", "10", "00")
    names(thetas) <- c("11", "10", "00")

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
